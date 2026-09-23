"""Retain lmp_mpi comparisons and fail if the declared tolerances are exceeded."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import numpy as np
from cases import CASES
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "runs" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
    )
    parser.add_argument("--executable", default=None)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    ff = ForceField.bundled("ffield.reax.ZnOH.2010")
    calc = Calculator(ff)
    cases = dict(CASES)
    for r in (1.5, 2.5, 3.0, 3.5, 4.0, 4.9999, 5.0001, 9.999, 10.0, 10.001):
        cases[f"stretch_{r}"] = (["Zn", "O"], [[0, 0, 0], [r, 0, 0]])
    relaxed = calc.relax(["Zn", "O"], [[0, 0, 0], [2.3, 0.1, 0.2]], force_tolerance=1e-5)
    cases["relaxed_zno"] = (["Zn", "O"], relaxed.positions)
    symbols20, x20 = CASES["cluster20"]
    relaxed20 = calc.relax(symbols20, x20)
    cases["relaxed_cluster20"] = (symbols20, relaxed20.positions)
    report = {
        "force_field_sha256": ff.checksum,
        "reference_version": None,
        "energy_tolerance_per_atom": 1e-5,
        "charge_tolerance": 1e-6,
        "force_tolerance": 1e-4,
        "relaxation_converged": relaxed.converged and relaxed20.converged,
        "cluster20_relaxation_iterations": relaxed20.iterations,
        "cluster20_relaxation_max_force": float(np.max(abs(relaxed20.evaluation.forces))),
        "force_comparison": "fixed_charge",
        "full_derivative": False,
        "cases": {},
    }
    for name, (symbols, positions) in cases.items():
        actual = calc.evaluate(symbols, positions)
        ref = evaluate_lammps(ff, symbols, positions, directory=args.output / name, executable=args.executable)
        report["reference_version"] = ref.version
        errors = {k: abs(actual.components[k] - ref.components[k]) for k in actual.components}
        row = {
            "atoms": len(symbols),
            "energy_error_per_atom": abs(actual.energy - ref.energy) / len(symbols),
            "component_errors": errors,
            "max_charge_error": float(np.max(abs(actual.charges - ref.charges))),
            "max_force_error": float(np.max(abs(actual.forces - ref.forces))),
        }
        row["passed"] = bool(
            row["energy_error_per_atom"] <= 1e-5
            and max(errors.values()) / len(symbols) <= 1e-5
            and row["max_charge_error"] <= 1e-6
            and row["max_force_error"] <= 1e-4
        )
        report["cases"][name] = row
        (args.output / name / "python.json").write_text(
            json.dumps(
                {
                    "energy": actual.energy,
                    "components": actual.components,
                    "charges": actual.charges.tolist(),
                    "forces": actual.forces.tolist(),
                    "full_derivative": actual.full_derivative,
                    "force_convention": actual.force_convention,
                },
                indent=2,
            )
            + "\n"
        )
        print(f"{name}: {'PASS' if row['passed'] else 'FAIL'}, force error {row['max_force_error']:.3g}", flush=True)
    report["passed"] = bool(report["relaxation_converged"] and all(row["passed"] for row in report["cases"].values()))
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Report: {args.output / 'summary.json'}")
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
