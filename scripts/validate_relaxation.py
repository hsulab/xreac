"""Retain fixed-charge FIRE relaxations and verify final forces with lmp_mpi."""

import argparse
import json
import os
from pathlib import Path
import sys
from importlib.metadata import version

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "tests"))

from cases import CASES
from water_cluster import comparison, serialize, water_cases
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("ase", "native"), default="ase")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        args.output = ROOT / "validation" / f"relaxation-{args.backend}-catorch3"
    args.output.mkdir(parents=True, exist_ok=False)
    cases = [
        ("zno", ForceField.bundled("ffield.reax.ZnOH"), (["Zn", "O"], [[0, 0, 0], [2.3, 0.1, 0.2]]), 1e-5),
        ("cluster20", ForceField.bundled("ffield.reax.ZnOH"), CASES["cluster20"], 1e-4),
        ("water_monomer", ForceField.bundled("qeq_ff.water"), water_cases()["monomer"], 1e-5),
    ]
    report = {
        "optimizer": "FIRE",
        "backend": args.backend,
        "environment": os.environ.get("CONDA_DEFAULT_ENV"),
        "python": sys.version,
        "numpy": np.__version__,
        "autograd": version("autograd"),
        "ase": version("ase") if args.backend == "ase" else None,
        "full_derivative": False,
        "force_convention": "fixed_charge",
        "energy_units": "kcal/mol",
        "force_units": "kcal/mol/Angstrom",
        "cases": {},
    }
    for name, ff, (symbols, x), tolerance in cases:
        directory = args.output / name
        directory.mkdir()
        (directory / "initial.json").write_text(
            json.dumps({"symbols": symbols, "positions": np.asarray(x).tolist()}, indent=2) + "\n"
        )
        calculator = Calculator(ff)
        initial = calculator.evaluate(symbols, x)
        relaxed = calculator.relax(symbols, x, force_tolerance=tolerance, backend=args.backend)
        ref = evaluate_lammps(ff, symbols, relaxed.positions, directory=directory / "lammps")
        (directory / "python.json").write_text(json.dumps(serialize(relaxed.evaluation), indent=2) + "\n")
        (directory / "reference.json").write_text(json.dumps(serialize(ref), indent=2) + "\n")
        xyz = [str(len(symbols)), "FIRE relaxation with fixed-charge forces; Angstrom"]
        xyz += [s + " " + " ".join(f"{v:.17g}" for v in row) for s, row in zip(symbols, relaxed.positions)]
        (directory / "structure.xyz").write_text("\n".join(xyz) + "\n")
        row = comparison(relaxed.evaluation, ref, len(symbols))
        row.update(
            {
                "iterations": relaxed.iterations,
                "converged": relaxed.converged,
                "message": relaxed.message,
                "force_tolerance": tolerance,
                "max_iterations": 500,
                "initial_energy": initial.energy,
                "final_energy": relaxed.evaluation.energy,
                "max_force": float(np.max(abs(relaxed.evaluation.forces))),
                "reference_max_force": float(np.max(abs(ref.forces))),
                "reference_version": ref.version,
                "force_field_sha256": ff.checksum,
            }
        )
        row["passed"] = bool(
            row["passed"]
            and row["converged"]
            and not relaxed.evaluation.full_derivative
            and row["max_force"] <= tolerance
            and row["reference_max_force"] <= tolerance + 1e-7
        )
        report["cases"][name] = row
        print(
            f"{name}: {'PASS' if row['passed'] else 'FAIL'}; {relaxed.iterations} iterations; max force {row['max_force']:.3e}; LAMMPS {row['reference_max_force']:.3e}",
            flush=True,
        )
    report["passed"] = all(row["passed"] for row in report["cases"].values())
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
