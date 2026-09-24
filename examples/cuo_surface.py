"""Time one CuO(010) slab and always verify its results against lmp_mpi.

The bundled Cu/O/H/Cl parameters have a separate CC BY-NC 4.0 license; see data/README.md.
Use --source-root to measure an extracted historical src/ tree on the same slab.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
THREAD_ENV = {
    name: "1"
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    )
}
os.environ.update(THREAD_ENV)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ffield",
        type=Path,
        default=ROOT / "data/ffield.reax.CuOHCl.2010",
        help="Cu/O parameter file (default: bundled)",
    )
    parser.add_argument("--executable", default=os.environ.get("XREAC_LAMMPS", "lmp_mpi"))
    parser.add_argument(
        "--source-root", type=Path, default=ROOT, help="Repository or extracted tree containing src/xreac"
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--batch-seconds", type=float, default=0.3)
    parser.add_argument(
        "--lammps-calls", type=int, default=100, help="Fresh and reusable evaluations per timing batch"
    )
    parser.add_argument(
        "--skip-lammps-timing", action="store_true", help="Still perform LAMMPS numerical verification"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/runs" / datetime.now(timezone.utc).strftime("cuo-%Y%m%dT%H%M%SZ"),
    )
    args = parser.parse_args()
    if min(args.repeats, args.batch_seconds, args.lammps_calls) <= 0:
        parser.error("Repeat counts and batch duration must be positive")
    source = args.source_root.resolve() / "src"
    if not (source / "xreac/__init__.py").is_file():
        parser.error("--source-root must contain src/xreac")
    executable = shutil.which(args.executable)
    if executable is None:
        parser.error("LAMMPS is required for this example; use --executable or XREAC_LAMMPS")
    # Import the selected implementation before the shared repository helpers,
    # so comparisons against a historical revision use that revision's code.
    sys.path.insert(0, str(source))
    import xreac
    from xreac import Calculator, ForceField
    from xreac.reference import evaluate_lammps

    if Path(xreac.__file__).resolve().parent != source / "xreac":
        raise RuntimeError("Imported an unexpected xreac source tree")
    for path in ("scripts", "examples", "tests"):
        sys.path.insert(0, str(ROOT / path))
    import numpy as np
    from ase import Atoms
    from ase.io import write
    from benchmark_lammps import lammps_timing, native_timing
    from validate import backend_differences, validation_cases
    from water_cluster import comparison, serialize
    from xreac.ase import ReaxFFCalculator

    ff = ForceField.from_file(args.ffield)
    _, symbols, x, cell, pbc = validation_cases(include_cuo=True)["surface_cuo_010"]
    args.output.mkdir(parents=True, exist_ok=False)
    structure = dict(symbols=symbols, positions=x.tolist(), cell=cell.tolist(), pbc=pbc)
    (args.output / "structure.json").write_text(json.dumps(structure, indent=2) + "\n")
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
    write(args.output / "cuo_010.extxyz", atoms)
    result, timing = native_timing(Calculator(ff), symbols, x, cell, pbc, args.repeats, args.batch_seconds)
    atoms.calc = ReaxFFCalculator(ff, neighbor_backend="ase")
    atoms.get_forces()
    ase_result = atoms.calc.evaluation
    reference = evaluate_lammps(
        ff, symbols, x, cell=cell, pbc=pbc, executable=executable, directory=args.output / "reference"
    )
    checks = {"native": comparison(result, reference, len(x)), "ase": comparison(ase_result, reference, len(x))}
    neighbors = backend_differences(result, ase_result, len(x))
    passed = all(row["passed"] for row in checks.values()) and neighbors["passed"]
    passed &= reference.cell_repetitions == (1, 1, 1)
    numerical = dict(native=serialize(result), ase=serialize(ase_result), reference=serialize(reference))
    (args.output / "results.json").write_text(json.dumps(numerical, indent=2) + "\n")
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        case="surface_cuo_010",
        atoms=len(x),
        description="Unrelaxed, stoichiometric 96-atom CuO(010) slab; 12 A vacuum; 0.01 A seeded displacement",
        structure_source="https://doi.org/10.1103/PhysRevB.39.4343",
        force_field=str(args.ffield.resolve()),
        force_field_sha256=ff.checksum,
        source_root=str(args.source_root.resolve()),
        source_sha256={
            str(p.relative_to(source)): sha256(p.read_bytes()).hexdigest()
            for p in sorted((source / "xreac").glob("*.py"))
        },
        driver_revision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        platform=platform.platform(),
        python=sys.version,
        environment={name: version(name) for name in ("numpy", "autograd", "ase")},
        thread_environment=THREAD_ENV,
        full_derivative=False,
        force_convention="fixed_charge",
        xreac=timing,
        lammps_comparison=checks,
        neighbor_comparison=neighbors,
        reference_version=reference.version,
        reference_cell_repetitions=reference.cell_repetitions,
        passed=bool(passed),
    )
    if not passed:
        (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        raise SystemExit(f"CuO comparison FAILED; see {args.output}")
    print(
        f"CuO(010), {len(x)} atoms: LAMMPS and neighbor checks PASS; xreac {1000 * timing['median_seconds']:.3f} ms",
        flush=True,
    )
    for mode in () if args.skip_lammps_timing else ("fresh", "steady"):
        elapsed, energy, state = lammps_timing(reference.directory, executable, args.repeats, args.lammps_calls, mode)
        errors = dict(
            energy_per_atom=abs(energy - result.energy) / len(x),
            forces=float(np.max(abs(state[:, 2:] - result.forces))),
            charges=float(np.max(abs(state[:, 1] - result.charges))),
        )
        if any(errors[key] > checks["native"]["limits"][key] for key in errors):
            raise RuntimeError(f"Timed LAMMPS {mode} results differ: {errors}")
        elapsed["final_state_errors"] = errors
        report[f"lammps_{mode}"] = elapsed
        report[f"slowdown_vs_{mode}"] = timing["median_seconds"] / elapsed["median_seconds"]
        print(
            f"LAMMPS {mode}: {1000 * elapsed['median_seconds']:.3f} ms; xreac ratio {report[f'slowdown_vs_{mode}']:.2f}x",
            flush=True,
        )
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
