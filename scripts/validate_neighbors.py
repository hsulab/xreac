"""Compare ASE neighbors with native replication, optionally against lmp_mpi."""
import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for path in ("src", "examples", "tests"):
    sys.path.insert(0, str(ROOT/path))

from ase import Atoms
from cases import CASES
from periodic_water import periodic_cases
from small_cells import small_cell_cases
from water_cluster import comparison, serialize, water_cases
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps


def neighbor_cases():
    cases = {"zno_"+name: ("ffield.reax.ZnOH", s, np.asarray(x, dtype=float), None, False)
             for name, (s, x) in CASES.items()}
    cases.update({"cluster_water_"+name: ("qeq_ff.water", s, x, None, False)
                  for name, (s, x) in water_cases().items()})
    carbon = {
        "methane": (["C"]+["H"]*4, np.vstack(([0, 0, 0], .63*np.array(
            [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])))),
        "carbon_monoxide": (["C", "O"], [[0, 0, 0], [1.15, 0, 0]]),
        "carbon_dimer": (["C", "C"], [[0, 0, 0], [1.2, 0, 0]]),
        "carbon_torsion": (["C"]*4, [[0, 0, 0], [1.5, .1, 0], [2.3, 1.2, .2], [3.7, .9, .8]]),
    }
    cases.update({name: ("ffield.reax.cho", s, np.asarray(x, dtype=float), None, False)
                  for name, (s, x) in carbon.items()})
    cases.update({"periodic_"+name: ("qeq_ff.water", *values) for name, values in periodic_cases().items()})
    cases.update({"small_"+name: values for name, values in small_cell_cases().items()})
    cases["periodic_carbon_chain"] = ("ffield.reax.cho", ["C"]*3,
        np.array([[0., 0, 0], [1.5, .4, .2], [3., -.2, .8]]), np.diag([4.5, 12., 12.]), [True, False, False])
    for distance in (4.9999, 5., 5.0001, 9.9999, 10., 10.0001):
        cases[f"cutoff_{distance:g}"] = ("ffield.reax.ZnOH", ["Zn", "O"],
            np.array([[0., 0, 0], [distance, 0, 0]]), None, False)
    return cases


def backend_differences(actual, expected, atoms):
    """Strict floating-point equivalence, independently of LAMMPS tolerances."""
    errors = {"energy_per_atom": abs(actual.energy-expected.energy)/atoms,
              "components_per_atom": max(abs(actual.components[k]-expected.components[k])
                                         for k in actual.components)/atoms}
    for key in ("forces", "charges", "bond_orders", "total_bond_orders", "lone_pairs", "bond_counts", "dipole"):
        errors[key] = float(np.max(np.abs(getattr(actual, key)-getattr(expected, key))))
    tolerances = {"energy_per_atom": 1e-9, "components_per_atom": 1e-9, "forces": 2e-8,
                  "charges": 1e-10, "bond_orders": 1e-10, "total_bond_orders": 1e-10,
                  "lone_pairs": 1e-10, "bond_counts": 0, "dipole": 1e-9}
    return {"max_absolute_differences": errors, "tolerances": tolerances,
            "passed": all(errors[k] <= tolerance for k, tolerance in tolerances.items())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Also run fresh LAMMPS single points")
    parser.add_argument("--executable", help="Override lmp_mpi executable")
    parser.add_argument("--output", type=Path, default=ROOT/"validation"/"runs"/
                        ("neighbors-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "full_derivative": False, "force_convention": "fixed_charge",
               "reference_convention": "LAMMPS supercells normalized to input cell for small inputs",
               "units": {"energy": "kcal/mol/input-cell", "forces": "kcal/mol/Angstrom", "charges": "e"},
               "environment": {"python": sys.version, "numpy": version("numpy"),
                               "autograd": version("autograd"), "ase": version("ase")},
               "reference_verified": args.verify, "cases": {}}
    for name, (filename, symbols, x, cell, pbc) in neighbor_cases().items():
        ff = ForceField.bundled(filename)
        work = args.output/name
        work.mkdir()
        structure = dict(symbols=symbols, positions=x.tolist(),
                         cell=None if cell is None else cell.tolist(), pbc=pbc)
        (work/"structure.json").write_text(json.dumps(structure, indent=2)+"\n")
        results, timings = {}, {}
        atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
        for backend in ("replicated", "ase"):
            atoms.calc = ReaxFFCalculator(ff, neighbor_backend=backend)
            start = perf_counter()
            atoms.get_forces()
            timings[backend] = perf_counter()-start
            results[backend] = atoms.calc.evaluation
            data = serialize(results[backend])
            data["neighbor_backend"] = backend
            (work/(backend+".json")).write_text(json.dumps(data, indent=2)+"\n")
        row = {"atoms": len(x), "force_field": filename, "force_field_sha256": ff.checksum,
               "energy": results["ase"].energy, "evaluation_seconds": timings,
               "replicated_cell_repetitions": results["replicated"].cell_repetitions,
               "ase_cell_repetitions": results["ase"].cell_repetitions,
               "backend_comparison": backend_differences(results["ase"], results["replicated"], len(x))}
        row["passed"] = row["backend_comparison"]["passed"]
        if args.verify:
            ref = evaluate_lammps(ff, symbols, x, cell=cell, pbc=pbc,
                                  directory=work/"lammps", executable=args.executable)
            (work/"reference.json").write_text(json.dumps(serialize(ref), indent=2)+"\n")
            row["lammps_comparison"] = comparison(results["ase"], ref, len(x))
            row["reference_version"] = ref.version
            row["passed"] &= row["lammps_comparison"]["passed"]
        summary["cases"][name] = row
        force_error = row["backend_comparison"]["max_absolute_differences"]["forces"]
        print(f"{name}: {'PASS' if row['passed'] else 'FAIL'}, backend force difference {force_error:.3g}", flush=True)
    summary["passed"] = all(row["passed"] for row in summary["cases"].values())
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    if not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
