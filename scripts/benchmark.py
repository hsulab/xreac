"""Benchmark isolated clusters in separate processes (including peak process RSS)."""

import argparse
from datetime import datetime, timezone
import json
import platform
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))

import numpy as np
from importlib.metadata import version
from cases import cluster
from xreac import Calculator, ForceField
from xreac.energy import EnergyModel
from xreac.neighbors import replicated_neighbors


def worker(n, relax_iterations=0):
    ff = ForceField.bundled("ffield.reax.ZnOH.2010")
    calc = Calculator(ff)
    symbols, x = cluster(n)
    neighbors, _ = replicated_neighbors(x, ff.general[12])
    model = EnergyModel(ff, symbols, neighbors)
    # Benchmark only fixed-charge forces, matching the LAMMPS convention.
    start = time.perf_counter()
    calc.evaluate(symbols, x, full_derivative=False)
    first = time.perf_counter() - start
    energies, evaluations = [], []
    for _ in range(3):
        start = time.perf_counter()
        model.components(x)
        energies.append(time.perf_counter() - start)
        start = time.perf_counter()
        calc.evaluate(symbols, x, full_derivative=False)
        evaluations.append(time.perf_counter() - start)
    result = {
        "atoms": n,
        "first_evaluation_seconds": first,
        "median_energy_seconds": statistics.median(energies),
        "median_evaluation_seconds": statistics.median(evaluations),
        "full_derivative": False,
        "force_convention": "fixed_charge",
    }
    if relax_iterations:
        start = time.perf_counter()
        relaxed = calc.relax(symbols, x, max_iterations=relax_iterations, backend="native")
        result.update(
            relaxation_backend="native",
            relaxation_seconds=time.perf_counter() - start,
            relaxation_iterations=relaxed.iterations,
            relaxation_converged=relaxed.converged,
            relaxation_max_force=float(np.max(abs(relaxed.evaluation.forces))),
        )
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["peak_process_rss_mib"] = rss / (1024**2 if sys.platform == "darwin" else 1024)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=int)
    parser.add_argument(
        "--relax-iterations",
        type=int,
        default=0,
        help="Optionally time native FIRE for this many steps (fixed-charge forces only)",
    )
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "runs" / "benchmark.json")
    args = parser.parse_args()
    if args.relax_iterations < 0:
        parser.error("--relax-iterations must be nonnegative")
    if args.worker:
        print(json.dumps(worker(args.worker, args.relax_iterations)))
        return
    report = {
        "date": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "autograd": version("autograd"),
        "force_field_sha256": ForceField.bundled("ffield.reax.ZnOH.2010").checksum,
        "full_derivative": False,
        "force_convention": "fixed_charge",
        "relaxation_backend": "native" if args.relax_iterations else None,
        "relaxation_max_iterations": args.relax_iterations,
        "note": "Three timed repeats; fixed-charge forces only. Optional relaxation uses native FIRE. RSS includes interpreter, numerical libraries, and any requested relaxation.",
        "results": [],
    }
    for n in (20, 100, 200):
        child = subprocess.run(
            [sys.executable, __file__, "--worker", str(n), "--relax-iterations", str(args.relax_iterations)],
            capture_output=True,
            text=True,
            check=True,
        )
        row = json.loads(child.stdout)
        report["results"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
