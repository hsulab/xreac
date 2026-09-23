"""Fixed-cell periodic water examples, optionally verified against lmp_mpi.

The bulk case has 64 molecules (192 atoms) at approximately 1 g/cm^3.
These deterministic test structures are not equilibrated liquid snapshots.
"""

import argparse
from datetime import datetime, timezone
from itertools import product
import json
from importlib.metadata import version
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from water_cluster import comparison, serialize, water_cases
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def water_box(counts=(4, 4, 4), cell=None, pbc=(True, True, True)):
    """Place randomly oriented waters on a jittered grid, with a fixed seed."""
    cell = np.diag([12.48] * 3) if cell is None else np.asarray(cell, dtype=float)
    rng = np.random.default_rng(731)
    monomer = water_cases()["monomer"][1]
    positions = []
    for index in product(*(range(n) for n in counts)):
        center = ((np.asarray(index) + 0.12) / counts) @ cell + rng.normal(0, 0.04, 3)
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        positions.extend(monomer @ rotation + center)
    positions = np.asarray(positions)
    # Deliberately split some molecules across the primary-cell faces.
    fractional = positions @ np.linalg.inv(cell)
    fractional[:, pbc] %= 1
    positions = fractional @ cell
    return ["O", "H", "H"] * int(np.prod(counts)), positions


def periodic_cases():
    cell = np.diag([12.0] * 3)
    symbols, dimer = water_cases()["distorted_dimer"]
    boundary = (dimer + [11.5, 11.6, 11.8]) % 12
    separated = dimer.copy()
    separated[3:] += [3.15, 0.2, 0.3]
    tilted = np.array([[12.0, 0, 0], [3.0, 12.0, 0], [1.0, 2.0, 12.0]])
    rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
    tilted = tilted @ rotation
    sparse_symbols, sparse = water_box((2, 2, 2), tilted)
    _, slab = water_box((2, 2, 2), tilted, pbc=(True, True, False))
    bulk_symbols, bulk = water_box()
    return {
        "boundary_dimer": (symbols, boundary, cell, [True] * 3),
        "multiple_images": (symbols, separated, cell, [True] * 3),
        "triclinic_water": (sparse_symbols, sparse, tilted, [True] * 3),
        "partial_pbc_water": (sparse_symbols, slab, tilted, [True, True, False]),
        "bulk_water_192": (bulk_symbols, bulk, np.diag([12.48] * 3), [True] * 3),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--executable", help="Override lmp_mpi executable")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "runs" / ("periodic-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")),
    )
    args = parser.parse_args()
    ff = ForceField.bundled("qeq_ff.water")
    calc = Calculator(ff)
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {
        "force_field": ff.path.name,
        "sha256": ff.checksum,
        "full_derivative": False,
        "force_convention": "fixed_charge",
        "units": {"energy": "kcal/mol", "forces": "kcal/mol/Angstrom", "charges": "e"},
        "structures": "Deterministic test geometries; not equilibrated liquid snapshots",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": version("numpy"),
            "autograd": version("autograd"),
        },
        "reference_verified": args.verify,
        "cases": {},
    }
    for name, (symbols, positions, cell, pbc) in periodic_cases().items():
        directory = args.output / name
        directory.mkdir()
        structure = dict(symbols=symbols, positions=positions.tolist(), cell=cell.tolist(), pbc=pbc)
        (directory / "structure.json").write_text(json.dumps(structure, indent=2) + "\n")
        lattice = " ".join(f"{v:.17g}" for v in cell.ravel())
        flags = " ".join("T" if flag else "F" for flag in pbc)
        xyz = [str(len(symbols)), f'Lattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="{flags}"']
        xyz += [s + " " + " ".join(f"{v:.17g}" for v in x) for s, x in zip(symbols, positions)]
        (directory / "structure.extxyz").write_text("\n".join(xyz) + "\n")
        start = perf_counter()
        result = calc.evaluate(symbols, positions, cell=cell, pbc=pbc)
        elapsed = perf_counter() - start
        (directory / "python.json").write_text(json.dumps(serialize(result), indent=2) + "\n")
        row = {
            "atoms": len(symbols),
            "energy": result.energy,
            "evaluation_seconds": elapsed,
            "cell": cell.tolist(),
            "pbc": pbc,
            "total_charge": float(result.charges.sum()),
        }
        if args.verify:
            ref = evaluate_lammps(
                ff, symbols, positions, cell=cell, pbc=pbc, directory=directory / "lammps", executable=args.executable
            )
            (directory / "reference.json").write_text(json.dumps(serialize(ref), indent=2) + "\n")
            row.update(comparison(result, ref, len(symbols)))
            row["reference_version"] = ref.version
            row["energy_difference"] = result.energy - ref.energy
        summary["cases"][name] = row
        print(
            f"{name}: E={result.energy:.9f} kcal/mol, {elapsed:.3f}s"
            + (f", reference {'PASS' if row['passed'] else 'FAIL'}" if args.verify else ""),
            flush=True,
        )
    if args.verify:
        summary["passed"] = all(row["passed"] for row in summary["cases"].values())
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Results: {args.output}")
    if args.verify and not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
