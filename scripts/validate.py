"""Validate energies, forces, charges, and properties by chemical system."""

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import gzip
import json
import shutil
import tarfile
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for path in ("src", "examples", "tests"):
    sys.path.insert(0, str(ROOT / path))

from ase import Atoms
from cases import CASES
from periodic_water import periodic_cases
from small_cells import small_cell_cases
from water_cluster import comparison, serialize, water_cases
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps


CUO_FORCE_FIELD = "ffield.reax.CuOHCl.2010"


def cuo_surface_case():
    """One slightly perturbed, unrelaxed 96-atom tenorite (010) slab.

    Conventional C2/c cell: Yang et al., Phys. Rev. B 39, 4343 (1989),
    doi:10.1103/PhysRevB.39.4343. ASE places c then a in the surface plane.
    Two layers and a 2x3 surface repeat keep both periodic heights above 10 A.
    """
    from ase.build import surface
    from ase.spacegroup import crystal

    bulk = crystal(
        ["Cu", "O"],
        basis=[(0.25, 0.25, 0), (0, 0.416, 0.25)],
        spacegroup=15,
        cellpar=[4.6837, 3.4226, 5.1288, 90, 99.54, 90],
    )
    slab = surface(bulk, (0, 1, 0), layers=2, vacuum=12).repeat((2, 3, 1))
    # Break exact geometric degeneracies without introducing another fixture.
    slab.positions += np.random.default_rng(260923).normal(0, 0.01, slab.positions.shape)
    slab.wrap()
    return CUO_FORCE_FIELD, slab.get_chemical_symbols(), slab.positions, slab.cell.array, slab.pbc.tolist()


def validation_cases(include_bulk=False, *, include_cuo=False):
    cases = {
        "zno_" + name: ("ffield.reax.ZnOH.2010", s, np.asarray(x, dtype=float), None, False)
        for name, (s, x) in CASES.items()
    }
    cases.update(
        {"cluster_water_" + name: ("ffield.reax.HO.2015", s, x, None, False) for name, (s, x) in water_cases().items()}
    )
    carbon = {
        "methane": (
            ["C"] + ["H"] * 4,
            np.vstack(([0, 0, 0], 0.63 * np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]))),
        ),
        "carbon_monoxide": (["C", "O"], [[0, 0, 0], [1.15, 0, 0]]),
        "carbon_dimer": (["C", "C"], [[0, 0, 0], [1.2, 0, 0]]),
    }
    cases.update(
        {name: ("ffield.reax.CHO.2008", s, np.asarray(x, dtype=float), None, False) for name, (s, x) in carbon.items()}
    )
    cases.update(
        {"periodic_" + name: ("ffield.reax.HO.2015", *values) for name, values in periodic_cases(include_bulk).items()}
    )
    cases.update({"small_" + name: values for name, values in small_cell_cases().items()})
    cases["periodic_carbon_chain"] = (
        "ffield.reax.CHO.2008",
        ["C"] * 3,
        np.array([[0.0, 0, 0], [1.5, 0.4, 0.2], [3.0, -0.2, 0.8]]),
        np.diag([4.5, 12.0, 12.0]),
        [True, False, False],
    )
    if include_cuo:
        cases["surface_cuo_010"] = cuo_surface_case()
    return cases


def backend_differences(actual, expected, atoms):
    """Strict floating-point equivalence, independently of LAMMPS tolerances."""
    errors = {
        "energy_per_atom": abs(actual.energy - expected.energy) / atoms,
        "components_per_atom": max(abs(actual.components[k] - expected.components[k]) for k in actual.components)
        / atoms,
    }
    for key in ("forces", "charges", "bond_orders", "total_bond_orders", "lone_pairs", "bond_counts", "dipole"):
        errors[key] = float(np.max(np.abs(getattr(actual, key) - getattr(expected, key))))
    tolerances = {
        "energy_per_atom": 1e-9,
        "components_per_atom": 1e-9,
        "forces": 2e-8,
        "charges": 1e-10,
        "bond_orders": 1e-10,
        "total_bond_orders": 1e-10,
        "lone_pairs": 1e-10,
        "bond_counts": 0,
        "dipole": 1e-9,
    }
    return {
        "max_absolute_differences": errors,
        "tolerances": tolerances,
        "passed": all(errors[k] <= tolerance for k, tolerance in tolerances.items()),
    }


SYSTEMS = {
    "ffield.reax.HO.2015": "water",
    "ffield.reax.ZnOH.2010": "zno",
    "ffield.reax.CHO.2008": "cho",
    CUO_FORCE_FIELD: "cuo",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Also run fresh LAMMPS single points")
    parser.add_argument("--executable", help="Override lmp_mpi executable")
    parser.add_argument("--system", choices=tuple(SYSTEMS.values()), help="Validate one chemical system")
    parser.add_argument("--include-bulk", action="store_true", help="Also check the 192-atom water box")
    parser.add_argument("--cuo-force-field", type=Path, help="Include the CuO slab using a custom Cu/O parameter file")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "runs" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
    )
    args = parser.parse_args()
    include_cuo = args.system == "cuo" or args.cuo_force_field is not None
    if include_cuo and not args.verify:
        parser.error("CuO validation requires --verify for a fresh LAMMPS comparison")
    args.output.mkdir(parents=True, exist_ok=False)
    all_passed = True
    cases = validation_cases(args.include_bulk, include_cuo=include_cuo)
    systems = tuple(dict.fromkeys(SYSTEMS[case[0]] for case in cases.values()))
    for system in (args.system,) if args.system else systems:
        destination = args.output / system
        destination.mkdir()
        summary = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "system": system,
            "full_derivative": False,
            "force_convention": "fixed_charge",
            "reference_convention": "LAMMPS supercells normalized to input cell for small inputs",
            "units": {"energy": "kcal/mol/input-cell", "forces": "kcal/mol/Angstrom", "charges": "e"},
            "environment": {key: version(key) for key in ("numpy", "autograd", "ase")},
            "python": sys.version,
            "reference_verified": args.verify,
            "cases": {},
        }
        structures, numerical = {}, {}
        for name, (filename, symbols, x, cell, pbc) in cases.items():
            if SYSTEMS[filename] != system:
                continue
            ff = (
                ForceField.from_file(args.cuo_force_field)
                if filename == CUO_FORCE_FIELD and args.cuo_force_field is not None
                else ForceField.bundled(filename)
            )
            structures[name] = dict(
                symbols=symbols, positions=x.tolist(), cell=None if cell is None else cell.tolist(), pbc=pbc
            )
            results, timings = {}, {}
            atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
            for backend in ("replicated", "ase"):
                atoms.calc = ReaxFFCalculator(ff, neighbor_backend=backend)
                start = perf_counter()
                atoms.get_forces()
                timings[backend] = perf_counter() - start
                results[backend] = atoms.calc.evaluation
            numerical[name] = {backend: serialize(result) for backend, result in results.items()}
            row = {
                "atoms": len(x),
                "force_field": filename,
                "force_field_sha256": ff.checksum,
                "optional": name in ("periodic_bulk_water_192", "surface_cuo_010"),
                "energy": results["ase"].energy,
                "evaluation_seconds": timings,
                "backend_comparison": backend_differences(results["ase"], results["replicated"], len(x)),
            }
            row["passed"] = row["backend_comparison"]["passed"]
            if args.verify:
                ref = evaluate_lammps(
                    ff,
                    symbols,
                    x,
                    cell=cell,
                    pbc=pbc,
                    directory=destination / "reference" / name,
                    executable=args.executable,
                )
                numerical[name]["reference"] = serialize(ref)
                row["lammps_comparison"] = comparison(results["ase"], ref, len(x))
                row["reference_version"] = ref.version
                row["passed"] &= row["lammps_comparison"]["passed"]
            summary["cases"][name] = row
            print(f"{system}/{name}: {'PASS' if row['passed'] else 'FAIL'}", flush=True)
        summary["passed"] = all(row["passed"] for row in summary["cases"].values())
        all_passed &= summary["passed"]
        for filename, data in (("summary.json", summary), ("structures.json", structures)):
            (destination / filename).write_text(json.dumps(data, indent=2) + "\n")
        (destination / "results.json.gz").write_bytes(gzip.compress(json.dumps(numerical).encode(), mtime=0))
        if args.verify:
            with tarfile.open(destination / "reference.tar.gz", "w:gz") as archive:
                archive.add(destination / "reference", arcname="reference")
            shutil.rmtree(destination / "reference")
    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
