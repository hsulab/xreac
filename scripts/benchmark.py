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


def worker(n):
    ff = ForceField.zno()
    calc = Calculator(ff)
    symbols, x = cluster(n)
    model = EnergyModel(ff, symbols)
    # Benchmark only fixed-charge forces, matching the LAMMPS convention.
    start = time.perf_counter()
    calc.evaluate(symbols, x, full_derivative=False)
    first = time.perf_counter()-start
    energies, evaluations = [], []
    for _ in range(3):
        start = time.perf_counter()
        model.components(x)
        energies.append(time.perf_counter()-start)
        start = time.perf_counter()
        calc.evaluate(symbols, x, full_derivative=False)
        evaluations.append(time.perf_counter()-start)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {"atoms": n, "first_evaluation_seconds": first,
            "median_energy_seconds": statistics.median(energies),
            "median_evaluation_seconds": statistics.median(evaluations),
            "full_derivative": False, "force_convention": "fixed_charge",
            "peak_process_rss_mib": rss/(1024**2 if sys.platform == "darwin" else 1024)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=int)
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "benchmark-fixed-charge.json")
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker)))
        return
    report = {"date": datetime.now(timezone.utc).isoformat(), "python": sys.version,
              "platform": platform.platform(), "numpy": np.__version__,
              "autograd": version("autograd"), "force_field_sha256": ForceField.zno().checksum,
              "full_derivative": False, "force_convention": "fixed_charge",
              "note": "Three timed repeats; evaluations compute only fixed-charge forces. No charge-response forces or relaxation; RSS includes interpreter and numerical libraries.",
              "results": []}
    for n in (20, 100, 200):
        child = subprocess.run([sys.executable, __file__, "--worker", str(n)], capture_output=True, text=True, check=True)
        row = json.loads(child.stdout)
        report["results"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
