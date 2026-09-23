"""Compare single-CPU xreac evaluations with in-process lmp_mpi timings."""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time

# Set before importing numerical libraries; inherited by the LAMMPS process.
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

ROOT = Path(__file__).resolve().parents[1]
for path in ("src", "scripts", "tests", "examples"):
    sys.path.insert(0, str(ROOT / path))

import numpy as np
from cases import cluster
from validate import SYSTEMS, validation_cases
from water_cluster import comparison
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def distribution(samples, count):
    per_call = [value / count for value in samples]
    return {
        "calls_per_batch": count,
        "batch_seconds": samples,
        "median_seconds": statistics.median(per_call),
        "min_seconds": min(per_call),
        "max_seconds": max(per_call),
    }


def native_timing(calc, symbols, x, cell, pbc, repeats, target_seconds):
    def evaluate():
        return calc.evaluate(symbols, x, cell=cell, pbc=pbc, full_derivative=False)

    start = time.perf_counter()
    result = evaluate()
    first = time.perf_counter() - start
    count = max(1, math.ceil(target_seconds / first))
    wall, cpu = [], []
    for _ in range(repeats):
        cpu_start, start = time.process_time(), time.perf_counter()
        for _ in range(count):
            evaluate()
        wall.append(time.perf_counter() - start)
        cpu.append(time.process_time() - cpu_start)
    timing = distribution(wall, count)
    timing.update(first_seconds=first, cpu_to_wall_ratio=sum(cpu) / sum(wall))
    return result, timing


def lammps_timing(work, executable, repeats, count, mode):
    # Reuse the exact validated system, force field, and QEq tolerance.
    original = (work / "in.lammps").read_text()
    lines = original.splitlines()
    qeq = next(line for line in lines if line.startswith("fix charges "))
    setup = original.split("compute reax all pair reaxff")[0]
    script = (
        "echo none\n"
        + setup
        + """compute reax all pair reaxff
compute bondinfo all reaxff/atom
compute properties all reduce sum c_bondinfo[1] c_bondinfo[2] c_bondinfo[3]
compute dipole all dipole
thermo_style custom step pe c_reax[*] c_properties[*] c_dipole[*]
thermo_modify format float %.17g
thermo 1000000000
run 0
"""
    )
    if mode == "fresh":
        # New QEq fix clears extrapolation history. run 0 pre yes also rebuilds
        # neighbors, matching the stateless Calculator.evaluate API.
        operation = f"unfix charges\nset group all charge 0\n{qeq}\nrun 0 pre yes post no\n"
        script += "log none\n"
    else:
        # No integrator: same fixed geometry, QEq every step, neighbor reuse.
        # This is an optimistic throughput baseline, not an MD benchmark.
        operation = f"run {count} pre yes post yes\n"
    for repeat in range(repeats):
        script += f"variable start{repeat} timer\n"
        script += operation * count if mode == "fresh" else operation
        script += f"variable stop{repeat} timer\n"
        script += f'print "$(v_stop{repeat}-v_start{repeat}:%.17g)" append {mode}.times screen no\n'
    # Preserve the final computed state independently of the initial validation.
    script += f'print "$(pe:%.17g)" file {mode}.energy screen no\n'
    script += f"write_dump all custom {mode}.dump id q fx fy fz modify sort id format float %.17g\n"
    (work / f"{mode}.in").write_text(script)
    command = [executable, "-in", f"{mode}.in", "-log", f"{mode}.log", "-screen", "none", "-nocite"]
    start = time.perf_counter()
    proc = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=300, check=True)
    elapsed = time.perf_counter() - start
    (work / f"{mode}.stderr").write_text(proc.stderr)
    log = (work / f"{mode}.log").read_text()
    if "convergence failed" in log.lower():
        raise RuntimeError(f"QEq did not converge: {work}")
    if not re.search(r"1 MPI tasks x (?:no|1) OpenMP threads", log):
        raise RuntimeError(f"Could not confirm single-rank, single-thread LAMMPS: {work}")
    samples = np.loadtxt(work / f"{mode}.times", ndmin=1).tolist()
    timing = distribution(samples, count)
    timing["process_wall_seconds"] = elapsed
    timing["mpi_tasks"] = 1
    timing["openmp_threads"] = 0 if "x no OpenMP threads" in log else 1
    if mode == "steady":
        loops = re.findall(r"Loop time of ([\d.eE+-]+) on 1 procs for (\d+) steps", log)
        timing["loop_seconds_per_step"] = [float(t) / int(n) for t, n in loops if int(n)]
    rows = (work / f"{mode}.dump").read_text().splitlines()
    offset = next(i for i, line in enumerate(rows) if line.startswith("ITEM: ATOMS")) + 1
    state = np.loadtxt(rows[offset:], ndmin=2)
    return timing, float((work / f"{mode}.energy").read_text()), state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", default=os.environ.get("XREAC_LAMMPS", "lmp_mpi"))
    parser.add_argument("--include-large", action="store_true", help="Include 100/200-atom Zn/O and 192-atom water")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--batch-seconds", type=float, default=0.3)
    parser.add_argument("--lammps-calls", type=int, default=1000)
    parser.add_argument(
        "--skip-lammps-timing",
        action="store_true",
        help="Time only xreac; still verify each case with a fresh LAMMPS single point",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation" / "runs" / datetime.now(timezone.utc).strftime("cpu-benchmark-%Y%m%dT%H%M%SZ"),
    )
    args = parser.parse_args()
    if min(args.repeats, args.batch_seconds, args.lammps_calls) <= 0:
        parser.error("Repeat counts and batch duration must be positive")
    executable = shutil.which(args.executable)
    if executable is None:
        parser.error(f"Executable not found: {args.executable}")
    args.output.mkdir(parents=True, exist_ok=False)
    cases = validation_cases(args.include_large)
    names = ["cluster_water_monomer", "cluster_water_distorted_dimer", "methane", "zno_cluster20"]
    if args.include_large:
        for n in (100, 200):
            symbols, x = cluster(n)
            cases[f"zno_cluster{n}"] = ("ffield.reax.ZnOH.2010", symbols, x, None, False)
            names.append(f"zno_cluster{n}")
        names.append("periodic_bulk_water_192")
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
            for path in sorted((ROOT / "src" / "xreac").glob("*.py"))
        },
        "platform": platform.platform(),
        "cpu": (
            subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"], text=True).strip()
            if sys.platform == "darwin"
            else platform.processor()
        ),
        "python": sys.version,
        "environment": {name: version(name) for name in ("numpy", "autograd", "ase")},
        "thread_environment": THREAD_ENV,
        "executable": executable,
        "repeats": args.repeats,
        "lammps_timing": not args.skip_lammps_timing,
        "method": {
            "xreac": "Calculator.evaluate, native neighbors, fresh QEq, fixed-charge forces, all properties",
            "lammps_fresh": "Repeated run 0; rebuild neighbors and recreate QEq fix with zero charges each call; all properties",
            "lammps_steady": "Fixed geometry, no integrator, QEq each step with history and neighbor reuse; properties at batch endpoints",
            "timing": "Median batch wall time per call; force-field loading, process startup and final dumps excluded",
            "cpu_limit": "One MPI rank and numerical-library thread; no OS core affinity is imposed",
            "caveat": "Fresh LAMMPS includes command parsing/setup overhead; steady is an optimistic fixed-geometry baseline, not MD",
        },
        "cases": {},
    }
    for name in names:
        filename, symbols, x, cell, pbc = cases[name]
        ff = ForceField.bundled(filename)
        result, native = native_timing(Calculator(ff), symbols, x, cell, pbc, args.repeats, args.batch_seconds)
        work = args.output / name
        ref = evaluate_lammps(ff, symbols, x, cell=cell, pbc=pbc, executable=executable, directory=work)
        checks = comparison(result, ref, len(x))
        if not checks["passed"] or ref.cell_repetitions != (1, 1, 1):
            raise RuntimeError(f"Reference mismatch or unequal atom counts: {name}")
        row = {
            "system": SYSTEMS[filename],
            "atoms": len(x),
            "force_field": filename,
            "force_field_sha256": ff.checksum,
            "reference_version": ref.version,
            "comparison": checks,
            "xreac": native,
        }
        for mode in () if args.skip_lammps_timing else ("fresh", "steady"):
            timing, energy, state = lammps_timing(work, executable, args.repeats, args.lammps_calls, mode)
            errors = {
                "energy_per_atom": abs(energy - result.energy) / len(x),
                "charges": float(np.max(abs(state[:, 1] - result.charges))),
                "forces": float(np.max(abs(state[:, 2:] - result.forces))),
            }
            if any(errors[key] > checks["limits"][key] for key in errors):
                raise RuntimeError(f"Timed {mode} calculation mismatch: {name}: {errors}")
            timing["final_state_errors"] = errors
            row[f"lammps_{mode}"] = timing
            row[f"slowdown_vs_{mode}"] = native["median_seconds"] / timing["median_seconds"]
        report["cases"][name] = row
        (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        message = f"{name} ({len(x)} atoms): xreac {1000 * native['median_seconds']:.3f} ms"
        if args.skip_lammps_timing:
            message += "; LAMMPS verification PASS"
        else:
            message += (
                f"; LAMMPS fresh {1000 * row['lammps_fresh']['median_seconds']:.3f} ms "
                f"({row['slowdown_vs_fresh']:.1f}x); steady {1000 * row['lammps_steady']['median_seconds']:.3f} ms "
                f"({row['slowdown_vs_steady']:.1f}x)"
            )
        print(message, flush=True)
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
