"""Experiment with QEq factorization/history reuse during bulk-water Berendsen MD.

This standalone benchmark temporarily intercepts the existing fixed-charge QEq
solve. It does not change the public calculator or its default direct solver.
Every candidate is compared to a fresh direct solve on the same matrix.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter
from unittest.mock import patch
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / name) for name in ("examples", "scripts", "src", "tests")]
# Also establishes single-thread numerical-library environment before NumPy.
from water_md import initial_water, run_ase, verify_frames, THREAD_ENV

import numpy as np
from scipy.linalg import lu_factor, lu_solve, LinAlgWarning
from ase.io import write
from autograd.tracer import Box
from xreac import energy


class CachedQEq:
    """Reuse LU as an approximate inverse, with checked iterative refinement.

    The current KKT matrix is always rebuilt. Cached factors are only a
    preconditioner: corrections solve old_KKT * dx = rhs - current_KKT * x.
    A new LU is built if six corrections fail to converge. Previous solutions
    supply a linear predictor; neither charges nor the current matrix are frozen.
    This experimental solver only accepts ordinary NumPy arrays, not autodiff.
    """

    def __init__(self, tolerance=1e-12, max_corrections=6):
        self.tolerance = tolerance
        self.max_corrections = max_corrections
        self.factor = None
        self.history = []
        self.rebuilds = 0
        self.corrections = 0
        self.max_residual = 0.0

    def solve(self, matrix, rhs):
        if isinstance(matrix, Box) or isinstance(rhs, Box):
            raise TypeError("Experimental QEq cache supports fixed-charge forces only")
        if self.history and self.history[-1].shape != rhs.shape:
            self.factor, self.history = None, []
        solution = None
        if self.factor is not None:
            guess = self.history[-1].copy()
            if len(self.history) == 2:
                guess = 2 * guess - self.history[-2]
            for _ in range(self.max_corrections):
                residual = rhs - matrix @ guess
                if np.max(abs(residual)) <= self.tolerance:
                    solution = guess
                    break
                guess += lu_solve(self.factor, residual, check_finite=False)
                self.corrections += 1
        if solution is None:
            with warnings.catch_warnings():
                warnings.simplefilter("error", LinAlgWarning)
                self.factor = lu_factor(matrix, check_finite=True)
            solution = lu_solve(self.factor, rhs, check_finite=True)
            self.rebuilds += 1
        residual = float(np.max(abs(rhs - matrix @ solution)))
        if not np.isfinite(solution).all() or not np.isfinite(residual) or residual > self.tolerance:
            raise RuntimeError(f"QEq residual {residual} exceeds {self.tolerance}")
        self.max_residual = max(self.max_residual, residual)
        self.history = (self.history + [solution.copy()])[-2:]
        return solution


class PairedQEqProbe:
    """Compare cached and fresh solvers at every geometry; use cached charges."""

    def __init__(self, direct):
        self.direct = direct
        self.cache = CachedQEq()
        self.samples = []
        self.max_charge_difference = 0.0

    def __call__(self, matrix, rhs):
        # Only the QEq KKT solve is intercepted. Assert its neutral constraint.
        if isinstance(matrix, Box) or matrix.shape != (len(rhs), len(rhs)) or matrix[-1, -1] != 0:
            raise AssertionError("Expected a numerical QEq KKT matrix")
        np.testing.assert_array_equal(matrix[-1, :-1], np.ones(len(rhs) - 1))
        times, results = {}, {}
        operations = [("direct", self.direct), ("cached", self.cache.solve)]
        if len(self.samples) % 2:
            operations.reverse()
        for label, solve in operations:
            start = perf_counter()
            results[label] = solve(matrix, rhs)
            times[label] = perf_counter() - start
        error = float(np.max(abs(results["direct"][:-1] - results["cached"][:-1])))
        self.max_charge_difference = max(self.max_charge_difference, error)
        if error > 1e-10 or abs(results["cached"][:-1].sum()) > 1e-10:
            raise RuntimeError(f"Cached QEq differs from direct solve: {error}")
        self.samples.append(times)
        return results["cached"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--warmup", type=int, default=100)
    parser.add_argument("--timestep-fs", type=float, default=0.25)
    parser.add_argument("--temperature", type=float, default=300)
    parser.add_argument("--taut-fs", type=float, default=100)
    parser.add_argument("--seed", type=int, default=260924)
    parser.add_argument("--executable", default="lmp_mpi")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/runs" / datetime.now(timezone.utc).strftime("water-qeq-%Y%m%dT%H%M%SZ"),
    )
    args = parser.parse_args()
    if (
        args.steps < 2
        or args.warmup < 0
        or any(not np.isfinite(v) or v <= 0 for v in (args.timestep_fs, args.temperature, args.taut_fs))
    ):
        parser.error("Need steps >= 2, warmup >= 0, and positive finite timestep, temperature and taut")
    executable = shutil.which(args.executable)
    if executable is None:
        parser.error(f"Executable not found: {args.executable}")
    args.output.mkdir(parents=True, exist_ok=False)
    initial, ff = initial_water(args.temperature, args.seed)
    probe = PairedQEqProbe(np.linalg.solve)
    qeq_seconds = []
    electrostatics = energy.EnergyModel.electrostatics

    def timed_electrostatics(model, distances, fixed_charges=None):
        if fixed_charges is not None:
            return electrostatics(model, distances, fixed_charges)
        start = perf_counter()
        result = electrostatics(model, distances)
        qeq_seconds.append(perf_counter() - start)
        return result

    # Scoped to this isolated experimental process, restored before references.
    # The sole core linalg.solve call is the fixed-charge QEq KKT solve.
    with (
        patch.object(energy.np.linalg, "solve", probe),
        patch.object(energy.EnergyModel, "electrostatics", timed_electrostatics),
    ):
        frames, states, timing = run_ase(
            initial, ff, args.steps, args.warmup, args.timestep_fs, args.temperature, args.taut_fs
        )
    expected_calls = 1 + args.warmup + args.steps
    if len(probe.samples) != expected_calls or len(qeq_seconds) != expected_calls:
        raise RuntimeError(f"Expected {expected_calls} QEq solves, got {len(probe.samples)}")
    samples = probe.samples[1 + args.warmup :]
    direct = sum(row["direct"] for row in samples)
    cached = sum(row["cached"] for row in samples)
    qeq_total = sum(qeq_seconds[1 + args.warmup :])
    write(args.output / "samples.extxyz", frames)
    checks = verify_frames(frames, states, ff, executable, args.output / "references")
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        git_revision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        source_sha256={
            str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
            for p in [
                *sorted((ROOT / "src/xreac").glob("*.py")),
                ROOT / "examples/water_md.py",
                Path(__file__).resolve(),
            ]
        },
        case="periodic_bulk_water_192",
        atoms=len(initial),
        force_field_sha256=ff.checksum,
        settings=dict(
            steps=args.steps,
            warmup=args.warmup,
            timestep_fs=args.timestep_fs,
            temperature_K=args.temperature,
            taut_fs=args.taut_fs,
            seed=args.seed,
        ),
        thread_environment=THREAD_ENV,
        environment={name: version(name) for name in ("numpy", "autograd", "ase", "scipy")},
        method="Every moving-MD QEq matrix solved with direct NumPy and cached SciPy LU/refinement; alternate solver order. Cached charges drive MD. Rebuild matrix and converge every step. Full MD timer includes both solvers and verification overhead; savings below isolate solver time, not a measured production speedup.",
        md=timing,
        qeq=dict(
            calls=len(probe.samples),
            production_calls=len(samples),
            direct_seconds=direct,
            cached_seconds=cached,
            solver_speedup=direct / cached,
            estimated_md_saving_fraction=(direct - cached) / (timing["seconds"] - cached),
            ideal_solve_elimination_ceiling_fraction=direct / (timing["seconds"] - cached),
            estimated_direct_qeq_assembly_and_solve_seconds=qeq_total - cached,
            estimated_direct_qeq_fraction=(qeq_total - cached) / (timing["seconds"] - cached),
            lu_builds=probe.cache.rebuilds,
            refinement_corrections=probe.cache.corrections,
            max_absolute_kkt_residual=probe.cache.max_residual,
            residual_tolerance=probe.cache.tolerance,
            max_charge_difference=probe.max_charge_difference,
            direct_median_seconds=float(np.median([s["direct"] for s in samples])),
            cached_median_seconds=float(np.median([s["cached"] for s in samples])),
        ),
        states=states,
        checks=checks,
    )
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.output / "solver-times.json").write_text(json.dumps(probe.samples) + "\n")
    print(json.dumps(report["qeq"], indent=2), flush=True)
    print(f"All direct-solve and LAMMPS checks PASS. Results: {args.output}")


if __name__ == "__main__":
    main()
