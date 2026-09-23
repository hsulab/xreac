"""Calculate primitive periodic cells and verify against larger LAMMPS supercells.

Small LAMMPS primitive cells are not used as ground truth: their hydrogen-bond
exclusions can depend on the original atom IDs. Results use input-cell units.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from water_cluster import comparison, serialize, water_cases
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def small_cell_cases():
    water = "ffield.reax.HO.2015"
    symbols, x = water_cases()["monomer"]
    cases = {"water_4A": (water, symbols, x, np.diag([4.0] * 3), [True] * 3)}
    rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
    tilted = np.array([[4.5, 0, 0], [1.0, 5.0, 0], [0.6, 0.4, 5.5]]) @ rotation
    cases["partial_pbc_water"] = (water, symbols, x @ rotation, tilted, [True, True, False])
    cases["zno_4A"] = (
        "ffield.reax.ZnOH.2010",
        ["Zn", "O"],
        np.array([[0.0, 0, 0], [1.9, 0.1, 0.2]]),
        np.diag([4.0] * 3),
        [True] * 3,
    )
    cases["zinc_chain"] = (
        "ffield.reax.ZnOH.2010",
        ["Zn"],
        np.zeros((1, 3)),
        np.diag([2.5, 12.0, 12.0]),
        [True, False, False],
    )
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "validation"
        / "runs"
        / ("small-cell-support-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")),
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {
        "full_derivative": False,
        "force_convention": "fixed_charge",
        "reference": "LAMMPS on replicated supercells; results normalized to input cell",
        "units": {"energy": "kcal/mol/input-cell", "forces": "kcal/mol/Angstrom", "charges": "e"},
        "python": sys.version,
        "cases": {},
    }
    for name, (filename, symbols, x, cell, pbc) in small_cell_cases().items():
        work = args.output / name
        work.mkdir()
        ff = ForceField.bundled(filename)
        (work / "structure.json").write_text(
            json.dumps(dict(symbols=symbols, positions=x.tolist(), cell=cell.tolist(), pbc=pbc), indent=2) + "\n"
        )
        start = perf_counter()
        result = Calculator(ff).evaluate(symbols, x, cell=cell, pbc=pbc)
        elapsed = perf_counter() - start
        (work / "python.json").write_text(json.dumps(serialize(result), indent=2) + "\n")
        row = {
            "atoms": len(x),
            "internal_atoms": int(len(x) * np.prod(result.cell_repetitions)),
            "cell_repetitions": result.cell_repetitions,
            "force_field": filename,
            "force_field_sha256": ff.checksum,
            "energy": result.energy,
            "hydrogen_bond_energy": result.components["hydrogen_bond"],
            "evaluation_seconds": elapsed,
        }
        if args.verify:
            ref = evaluate_lammps(ff, symbols, x, cell=cell, pbc=pbc, directory=work / "lammps")
            (work / "reference.json").write_text(json.dumps(serialize(ref), indent=2) + "\n")
            row.update(comparison(result, ref, len(x)))
            row["reference_version"] = ref.version
        summary["cases"][name] = row
        print(
            f"{name}: E={result.energy:.9f}, {len(x)} atoms / {row['internal_atoms']} internal, "
            f"{elapsed:.3f}s" + (f", {'PASS' if row['passed'] else 'FAIL'}" if args.verify else ""),
            flush=True,
        )
    if args.verify:
        summary["passed"] = all(row["passed"] for row in summary["cases"].values())
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if args.verify and not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
