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
from ase.neighborlist import neighbor_list
from cases import CASES
from periodic_water import periodic_cases
from small_cells import small_cell_cases
from water_cluster import comparison, serialize, water_cases
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps


CUO_FORCE_FIELD = "ffield.reax.CuOHCl.2010"
CHOCL_FORCE_FIELD = "ffield.reax.CHOCl.2021"
PILOT_CASES = ("periodic_bulk_water_192", "periodic_bulk_zno_128", "surface_cuo_010")


def zno_bulk_case():
    """One 128-atom wurtzite benchmark cell; dimensions all exceed 10 Angstrom.

    Representative fixture parameters, not a relaxed equilibrium prediction:
    a=3.25 Angstrom, c=5.21 Angstrom, u=0.382, with seeded 0.01 Angstrom noise.
    """
    from ase.build import bulk

    atoms = bulk("ZnO", "wurtzite", a=3.25, c=5.21, u=0.382, orthorhombic=True).repeat((4, 2, 2))
    atoms.positions += np.random.default_rng(260924).normal(0, 0.01, atoms.positions.shape)
    atoms.wrap()
    return "ffield.reax.ZnOH.2010", atoms.get_chemical_symbols(), atoms.positions, atoms.cell.array, atoms.pbc.tolist()


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
    if include_bulk:
        cases["periodic_bulk_zno_128"] = zno_bulk_case()
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


def charged_validation_cases():
    """Shared ionic examples and charged variants of existing periodic fixtures.

    Values extend the neutral case tuple with a net charge in e.
    """
    ff = "ffield.reax.HO.2015"
    symbols, water = water_cases()["monomer"]
    cases = {
        "hydroxide": (ff, symbols[:2], water[:2].copy(), None, False, -1.0),
        "hydronium": (
            ff,
            symbols + ["H"],
            np.vstack((water, [-0.3, -0.45, 0.85])),
            None,
            False,
            1.0,
        ),
    }
    for name in ("water_4A", "partial_pbc_water"):
        cases["charged_" + name] = (*small_cell_cases()[name], -1.0)
    return cases


def independent_charge_check(ff, symbols, x, cell, pbc, total_charge, charges):
    """Assemble QEq independently with ASE images; eliminate the last charge.

    This checks the QEq quadratic, not stationarity of the reported ReaxFF
    energy, which uses slightly different conversion constants.
    """
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
    cutoff = ff.general[12]
    i, j, distances = neighbor_list("ijd", atoms, cutoff)
    scaled = distances / cutoff
    taper = 1 - 35 * scaled**4 + 84 * scaled**5 - 70 * scaled**6 + 20 * scaled**7
    gamma = np.array([ff.atoms[s]["gamma"] for s in symbols])
    shield = taper / (distances**3 + (gamma[i] * gamma[j]) ** -1.5) ** (1 / 3)
    hessian = np.diag([ff.atoms[s]["eta"] for s in symbols])
    np.add.at(hessian, (i, j), 14.4 * shield)
    chi = np.array([ff.atoms[s]["chi"] for s in symbols])
    n = len(symbols)
    particular = np.zeros(n)
    particular[-1] = total_charge
    # q = particular + basis @ z enforces the constraint exactly.
    basis = np.vstack((np.eye(n - 1), -np.ones(n - 1)))
    reduced = basis.T @ hessian @ basis
    expected = particular + basis @ np.linalg.solve(reduced, -basis.T @ (hessian @ particular + chi))
    potential = hessian @ charges + chi
    errors = {
        "total_charge": float(abs(np.sum(charges) - total_charge)),
        "charges": float(np.max(abs(charges - expected))),
        "chemical_potential_spread": float(np.ptp(potential)),
    }
    limits = {"total_charge": 1e-8, "charges": 1e-10, "chemical_potential_spread": 1e-9}
    minimum_curvature = float(np.linalg.eigvalsh(reduced).min()) if n > 1 else None
    return {
        "method": "Independent ASE image assembly and charge-constraint elimination",
        "maximum_absolute_errors": errors,
        "limits": limits,
        "minimum_reduced_curvature": minimum_curvature,
        "passed": all(errors[k] <= limits[k] for k in errors) and (minimum_curvature is None or minimum_curvature > 0),
    }


def chocl_validation_cases():
    """Hur parameters: representative geometries, not optimized barriers.

    SN2 distances follow Bucko (2008) Table 1; hydrogen orientations are idealized.
    All values extend the neutral validation tuple with net charge in e.
    """
    phi = np.arange(3) * 2 * np.pi / 3
    h = np.column_stack((np.cos(phi), np.sin(phi), np.zeros(3)))
    tetrahedral_h = h * (1.10 * np.sqrt(8 / 9)) + [0, 0, -1.10 / 3]
    methyl = np.vstack(([0, 0, 0], tetrahedral_h, [0, 0, 1.85]))
    return {
        "chocl_chloromethane": (CHOCL_FORCE_FIELD, ["C", "H", "H", "H", "Cl"], methyl, None, False, 0.0),
        "chocl_hydrogen_chloride": (
            CHOCL_FORCE_FIELD,
            ["H", "Cl"],
            np.array([[0.0, 0, 0], [1.28, 0, 0]]),
            None,
            False,
            0.0,
        ),
        "chocl_formyl_chloride": (
            CHOCL_FORCE_FIELD,
            ["C", "H", "O", "Cl"],
            np.array([[0.0, 0, 0], [-0.55, -0.953, 0], [1.2, 0, 0], [-0.9, 1.559, 0]]),
            None,
            False,
            0.0,
        ),
        "chocl_sn2_reactant": (
            CHOCL_FORCE_FIELD,
            ["C", "H", "H", "H", "Cl", "Cl"],
            np.vstack((methyl, [0, 0, -3.15])),
            None,
            False,
            -1.0,
        ),
        "chocl_sn2_symmetric": (
            CHOCL_FORCE_FIELD,
            ["C", "H", "H", "H", "Cl", "Cl"],
            np.vstack(([0, 0, 0], 1.08 * h, [0, 0, 2.35], [0, 0, -2.35])),
            None,
            False,
            -1.0,
        ),
    }


SYSTEMS = {
    "ffield.reax.HO.2015": "water",
    "ffield.reax.ZnOH.2010": "zno",
    "ffield.reax.CHO.2008": "cho",
    CUO_FORCE_FIELD: "cuo",
    CHOCL_FORCE_FIELD: "cho",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Also run fresh LAMMPS single points")
    parser.add_argument("--include-charged", action="store_true", help="Also validate net charge with independent QEq")
    parser.add_argument("--include-chocl", action="store_true", help="Include the bundled Hur C/H/O/Cl parameters")
    parser.add_argument("--executable", help="Override lmp_mpi executable")
    parser.add_argument("--system", choices=tuple(SYSTEMS.values()), help="Validate one chemical system")
    parser.add_argument(
        "--include-bulk", action="store_true", help="Also check 192-atom water and 128-atom ZnO bulk cells"
    )
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
    charged = charged_validation_cases() if args.include_charged else {}
    if args.include_chocl:
        charged.update(chocl_validation_cases())
    cases.update({name: values[:5] for name, values in charged.items()})
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
            total_charge = charged[name][5] if name in charged else 0.0
            structures[name] = dict(
                symbols=symbols,
                positions=x.tolist(),
                cell=None if cell is None else cell.tolist(),
                pbc=pbc,
                total_charge=total_charge,
            )
            results, timings = {}, {}
            atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
            for backend in ("replicated", "ase"):
                atoms.calc = ReaxFFCalculator(ff, neighbor_backend=backend, total_charge=total_charge)
                start = perf_counter()
                atoms.get_forces()
                timings[backend] = perf_counter() - start
                results[backend] = atoms.calc.evaluation
            numerical[name] = {backend: serialize(result) for backend, result in results.items()}
            row = {
                "atoms": len(x),
                "force_field": filename,
                "force_field_sha256": ff.checksum,
                "optional": name in PILOT_CASES,
                "energy": results["ase"].energy,
                "evaluation_seconds": timings,
                "backend_comparison": backend_differences(results["ase"], results["replicated"], len(x)),
            }
            row["passed"] = row["backend_comparison"]["passed"]
            row["total_charge"] = total_charge
            row["reference_verified"] = args.verify
            row["reference_charge_mode"] = "supplied" if total_charge != 0 else "equilibrated"
            if total_charge != 0:
                row["independent_qeq"] = independent_charge_check(
                    ff, symbols, x, cell, pbc, total_charge, results["ase"].charges
                )
                row["passed"] &= row["independent_qeq"]["passed"]
                row["reference_note"] = "LAMMPS uses supplied charges for energy/forces; QEq validated independently"
            if args.verify:
                ref = evaluate_lammps(
                    ff,
                    symbols,
                    x,
                    cell=cell,
                    pbc=pbc,
                    directory=destination / "reference" / name,
                    executable=args.executable,
                    supplied_charges=results["ase"].charges if total_charge != 0 else None,
                )
                numerical[name]["reference"] = serialize(ref)
                row["lammps_comparison"] = comparison(results["ase"], ref, len(x))
                row["reference_version"] = ref.version
                row["passed"] &= row["lammps_comparison"]["passed"]
            summary["cases"][name] = row
            print(f"{system}/{name}: {'PASS' if row['passed'] else 'FAIL'}", flush=True)
        summary["passed"] = all(row["passed"] for row in summary["cases"].values())
        summary["reference_verified"] = all(row["reference_verified"] for row in summary["cases"].values())
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
