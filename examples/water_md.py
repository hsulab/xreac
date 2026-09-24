"""Single-CPU bulk-water Berendsen MD with mandatory LAMMPS verification."""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from time import perf_counter, process_time

# Set before importing NumPy/ASE; inherited by the single-rank LAMMPS process.
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
for directory in ("src", "scripts", "tests", "examples"):
    sys.path.insert(0, str(ROOT / directory))

import numpy as np
from ase import Atoms, units
from ase.io import write
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from validate import validation_cases
from water_cluster import comparison, serialize
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps, REFERENCE_VERSION


def initial_water(temperature, seed):
    filename, symbols, x, cell, pbc = validation_cases(include_bulk=True)["periodic_bulk_water_192"]
    ff = ForceField.bundled(filename)
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
    # ASE's default elemental masses differ slightly from this ReaxFF file.
    atoms.set_masses([ff.atoms[s]["mass"] for s in symbols])
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature, force_temp=True, rng=np.random.default_rng(seed))
    Stationary(atoms, preserve_temperature=True)
    return atoms, ff


def run_ase(initial, ff, steps, warmup, timestep_fs, temperature, taut_fs):
    atoms = initial.copy()
    atoms.calc = ReaxFFCalculator(ff, neighbor_backend="ase", neighbor_skin=0.3)
    frames, states, blocks = [], [], []
    with NVTBerendsen(
        atoms,
        timestep=timestep_fs * units.fs,
        temperature_K=temperature,
        taut=taut_fs * units.fs,
        fixcm=False,
    ) as md:
        atoms.get_forces()
        md.run(warmup)
        warmup_builds = atoms.calc.neighbor_list_builds
        for count in (0, steps // 2, steps - steps // 2):
            if count:
                cpu, start = process_time(), perf_counter()
                md.run(count)
                elapsed = perf_counter() - start
                blocks.append(dict(steps=count, seconds=elapsed, cpu_seconds=process_time() - cpu))
            atoms.get_forces()
            if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.get_temperature()):
                raise RuntimeError("Non-finite MD state")
            frames.append(atoms.copy())
            states.append(
                dict(
                    step=md.nsteps - warmup,
                    temperature_K=atoms.get_temperature(),
                    result=serialize(atoms.calc.evaluation),
                )
            )
            print(f"ASE step {md.nsteps - warmup}/{steps}: T={atoms.get_temperature():.2f} K", flush=True)
    elapsed = sum(block["seconds"] for block in blocks)
    return (
        frames,
        states,
        dict(
            seconds=elapsed,
            seconds_per_step=elapsed / steps,
            ns_per_day=steps * timestep_fs / 1e6 / elapsed * 86400,
            cpu_to_wall_ratio=sum(block["cpu_seconds"] for block in blocks) / elapsed,
            blocks=blocks,
            neighbor_builds_during_production=atoms.calc.neighbor_list_builds - warmup_builds,
            neighbor_builds_total=atoms.calc.neighbor_list_builds,
        ),
    )


def read_dump(path):
    lines = path.read_text().splitlines()
    frames = []
    for index, line in enumerate(lines):
        if line == "ITEM: TIMESTEP":
            step = int(lines[index + 1])
            count = int(lines[index + 3])
        elif line.startswith("ITEM: ATOMS"):
            columns = line.split()[2:]
            rows = np.asarray([[float(x) for x in row.split()] for row in lines[index + 1 : index + 1 + count]])
            rows = rows[np.argsort(rows[:, columns.index("id")])]
            frames.append((step, {name: rows[:, i] for i, name in enumerate(columns)}))
    return frames


def run_lammps(initial, initial_reference, steps, warmup, timestep_fs, temperature, taut_fs, executable, work):
    work.mkdir()
    for name in ("atoms.data", "ffield"):
        shutil.copyfile(initial_reference / name, work / name)
    # LAMMPS real units use Angstrom/fs; ASE velocities use internal time units.
    velocities = initial.get_velocities() * units.fs
    with (work / "atoms.data").open("a") as stream:
        stream.write("\nVelocities\n\n")
        for i, velocity in enumerate(velocities):
            stream.write(f"{i + 1} " + " ".join(f"{v:.17g}" for v in velocity) + "\n")
    setup = (initial_reference / "in.lammps").read_text().split("compute reax all pair reaxff")[0]
    setup = setup.replace("neighbor 2.0 bin", "neighbor 0.6 bin")
    interval = steps // 2
    script = (
        setup
        + f"""timestep {timestep_fs:.17g}
fix integrator all nve
fix thermostat all temp/berendsen {temperature:.17g} {temperature:.17g} {taut_fs:.17g}
compute_modify thermostat_temp extra/dof 0
thermo_style custom step temp pe
compute_modify thermo_temp extra/dof 0
thermo_modify format float %.17g
thermo {interval}
run {warmup}
reset_timestep 0
dump samples all custom {interval} md.dump id xu yu zu q fx fy fz vx vy vz
dump_modify samples sort id format float %.17g
run {steps}
write_dump all custom final.dump id xu yu zu q fx fy fz vx vy vz modify sort id format float %.17g
"""
    )
    (work / "in.md").write_text(script)
    command = [executable, "-in", "in.md", "-log", "log.md", "-screen", "none", "-nocite"]
    proc = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=max(300, steps), check=False)
    (work / "stderr.txt").write_text(proc.stderr)
    log = (work / "log.md").read_text()
    if proc.returncode or "convergence failed" in log.lower() or "ERROR" in log:
        raise RuntimeError(f"LAMMPS MD failed; see {work}: {log[-2000:]}")
    if REFERENCE_VERSION not in log or not re.search(r"1 MPI tasks x (?:no|1) OpenMP threads", log):
        raise RuntimeError("Unexpected LAMMPS version or thread/rank count")
    loops = re.findall(r"Loop time of ([\d.eE+-]+) on \d+ procs for (\d+) steps", log)
    elapsed, count = float(loops[-1][0]), int(loops[-1][1])
    if count != steps:
        raise RuntimeError("Unexpected number of timed LAMMPS steps")
    frames = read_dump(work / "md.dump")
    final = read_dump(work / "final.dump")[0]
    if frames[-1][0] != steps:
        frames.append(final)
    # Keep only start, middle and end, including for odd step counts.
    frames = [frame for frame in frames if frame[0] in (0, steps // 2, steps)]
    return frames, dict(
        seconds=elapsed, seconds_per_step=elapsed / steps, ns_per_day=steps * timestep_fs / 1e6 / elapsed * 86400
    )


def verify_frames(frames, states, ff, executable, work):
    checks = []
    for atoms, state in zip(frames, states):
        step = state["step"]
        ref = evaluate_lammps(
            ff,
            atoms.get_chemical_symbols(),
            atoms.positions,
            cell=atoms.cell.array,
            pbc=atoms.pbc,
            executable=executable,
            directory=work / f"ase-{step}",
        )
        # Re-evaluate with a fresh calculator; compare saved MD state as well.
        atoms.calc = ReaxFFCalculator(ff)
        atoms.get_forces()
        result = atoms.calc.evaluation
        check = comparison(result, ref, len(atoms))
        saved = state["result"]
        errors = {
            key: float(np.max(np.abs(np.asarray(saved[key]) - getattr(result, key))))
            for key in ("energy", "forces", "charges")
        }
        if not check["passed"] or any(value > 2e-8 for value in errors.values()):
            raise RuntimeError(f"ASE MD reference mismatch at step {step}: {check}, {errors}")
        (work / f"ase-{step}" / "results.json").write_text(
            json.dumps(dict(ase=saved, reference=serialize(ref))) + "\n"
        )
        checks.append(dict(step=step, comparison=check, md_vs_fresh_errors=errors))
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--warmup", type=int, default=100)
    parser.add_argument("--timestep-fs", type=float, default=0.25)
    parser.add_argument("--temperature", type=float, default=300)
    parser.add_argument("--taut-fs", type=float, default=100)
    parser.add_argument("--seed", type=int, default=260924)
    parser.add_argument("--executable", default=os.environ.get("XREAC_LAMMPS", "lmp_mpi"))
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/runs" / datetime.now(timezone.utc).strftime("water-md-%Y%m%dT%H%M%SZ"),
    )
    args = parser.parse_args()
    if (
        args.steps < 2
        or args.warmup < 0
        or any(not np.isfinite(v) or v <= 0 for v in (args.timestep_fs, args.temperature, args.taut_fs))
    ):
        parser.error("Need steps >= 2, warmup >= 0, and finite positive timestep, temperature and taut")
    executable = shutil.which(args.executable)
    if executable is None:
        parser.error(f"Executable not found: {args.executable}")
    args.output.mkdir(parents=True, exist_ok=False)
    initial, ff = initial_water(args.temperature, args.seed)
    write(args.output / "initial.extxyz", initial)
    ref = evaluate_lammps(
        ff,
        initial.get_chemical_symbols(),
        initial.positions,
        cell=initial.cell.array,
        pbc=initial.pbc,
        executable=executable,
        directory=args.output / "initial-reference",
    )
    settings = (args.steps, args.warmup, args.timestep_fs, args.temperature, args.taut_fs)
    frames, states, ase_timing = run_ase(initial, ff, *settings)
    write(args.output / "ase-samples.extxyz", frames)
    checks = verify_frames(frames, states, ff, executable, args.output / "references")
    lmp_frames, lmp_timing = run_lammps(initial, ref.directory, *settings, executable, args.output / "lammps-md")
    lmp_checks = []
    for step, data in lmp_frames:
        atoms = initial.copy()
        atoms.positions[:] = np.column_stack([data[key] for key in ("xu", "yu", "zu")])
        atoms.calc = ReaxFFCalculator(ff)
        atoms.get_forces()
        result = atoms.calc.evaluation
        reference = evaluate_lammps(
            ff,
            atoms.get_chemical_symbols(),
            atoms.positions,
            cell=atoms.cell.array,
            pbc=atoms.pbc,
            executable=executable,
            directory=args.output / "references" / f"lammps-{step}",
        )
        check = comparison(result, reference, len(atoms))
        errors = dict(
            charges=float(np.max(abs(data["q"] - result.charges))),
            forces=float(np.max(abs(np.column_stack([data[key] for key in ("fx", "fy", "fz")]) - result.forces))),
        )
        if not check["passed"] or any(errors[key] > check["limits"][key] for key in errors):
            raise RuntimeError(f"LAMMPS MD state mismatch at {step}: {check}, {errors}")
        lmp_checks.append(dict(step=step, comparison=check, md_vs_xreac_errors=errors))
        (args.output / "references" / f"lammps-{step}" / "results.json").write_text(
            json.dumps(dict(ase=serialize(result), reference=serialize(reference))) + "\n"
        )
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        git_revision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        source_sha256={
            str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
            for p in [*sorted((ROOT / "src/xreac").glob("*.py")), Path(__file__).resolve()]
        },
        platform=platform.platform(),
        environment={name: version(name) for name in ("numpy", "autograd", "ase")},
        thread_environment=THREAD_ENV,
        case="periodic_bulk_water_192",
        atoms=len(initial),
        force_field=ff.path.name,
        force_field_sha256=ff.checksum,
        settings=dict(
            steps=args.steps,
            warmup=args.warmup,
            timestep_fs=args.timestep_fs,
            temperature_K=args.temperature,
            taut_fs=args.taut_fs,
            seed=args.seed,
        ),
        method=dict(
            neighbors="ASE skin=0.3 A per atom, LAMMPS pair skin=0.6 A; displacement-based rebuilds",
            qeq="Fresh direct xreac solve each step; LAMMPS qeq/reaxff each step with history, tolerance 1e-12",
            timing="Moving MD, warmup and reference checks excluded; ASE wall time and LAMMPS internal loop time; one numerical thread/rank, no affinity",
            thermostat="Fixed volume Berendsen, no constraints, matched masses and initial velocities; 3N temperature DOF and no ongoing COM removal",
            caveat="ASE scales velocities before Verlet; LAMMPS scales after. Validate sampled identical geometries, not trajectory identity. Unrelaxed initial fixture; this is not an equilibrated-liquid study.",
        ),
        ase=ase_timing,
        lammps=lmp_timing,
        slowdown=ase_timing["seconds_per_step"] / lmp_timing["seconds_per_step"],
        ase_states=states,
        ase_checks=checks,
        lammps_checks=lmp_checks,
    )
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"ASE {ase_timing['seconds_per_step'] * 1000:.3f} ms/step; LAMMPS {lmp_timing['seconds_per_step'] * 1000:.3f} ms/step; ratio {report['slowdown']:.2f}x; all checks PASS",
        flush=True,
    )
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
