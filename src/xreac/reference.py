"""Development-only single-rank lmp_mpi reference harness."""
from dataclasses import dataclass
from pathlib import Path
import json
import os
import shutil
import subprocess
import tempfile

import numpy as np

from .calculator import validate_input
from .energy import COMPONENTS

REFERENCE_VERSION = "LAMMPS (22 Jul 2025 - Update 4)"


@dataclass(frozen=True)
class ReferenceResult:
    energy: float
    forces: np.ndarray
    charges: np.ndarray
    components: dict[str, float]
    version: str
    directory: Path


def evaluate_lammps(force_field, symbols, positions, *, executable=None,
                    directory=None, timeout=120, expected_version=REFERENCE_VERSION):
    """Run and retain a reproducible reference calculation; never silently skip."""
    symbols, x = validate_input(symbols, positions)
    executable = executable or os.environ.get("XREAC_LAMMPS", "lmp_mpi")
    executable = shutil.which(str(executable))
    if executable is None:
        raise FileNotFoundError("lmp_mpi unavailable; set XREAC_LAMMPS to its executable path")
    work = Path(directory) if directory is not None else Path(tempfile.mkdtemp(prefix="xreac-lmp-"))
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        raise ValueError(f"Reference directory must be empty: {work}")
    work = work.resolve()
    # Keep the exact input file alongside outputs, including its citation header.
    raw = force_field.path.read_bytes()
    from hashlib import sha256
    if sha256(raw).hexdigest() != force_field.checksum:
        raise ValueError("Force-field file changed after loading")
    (work / "ffield").write_bytes(raw)
    types = list(dict.fromkeys(symbols))
    lines = ["xreac reference", "", f"{len(x)} atoms", f"{len(types)} atom types", ""]
    for axis, dim in enumerate("xyz"):
        lines.append(f"{x[:, axis].min()-15:.17g} {x[:, axis].max()+15:.17g} {dim}lo {dim}hi")
    lines += ["", "Masses", ""]
    lines += [f"{i+1} {force_field.atoms[s]['mass']:.17g}" for i, s in enumerate(types)]
    lines += ["", "Atoms # charge", ""]
    lines += [f"{i+1} {types.index(s)+1} 0 " + " ".join(f"{v:.17g}" for v in pos)
              for i, (s, pos) in enumerate(zip(symbols, x))]
    (work / "atoms.data").write_text("\n".join(lines)+"\n")
    terms = " ".join(f"$(c_reax[{i}]:%.17g)" for i in range(1, 15))
    script = f"""units real
atom_style charge
boundary f f f
read_data atoms.data
pair_style reaxff NULL tabulate 0 enobonds yes
pair_coeff * * ffield {' '.join(types)}
fix charges all qeq/reaxff 1 0 {force_field.general[12]:.17g} 1e-12 reaxff maxiter 2000
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
compute reax all pair reaxff
thermo_style custom step pe c_reax[*]
thermo_modify format float %.17g
run 0
print "$(pe:%.17g) {terms}" file energy.txt screen no
write_dump all custom atoms.dump id type q x y z fx fy fz modify sort id format float %.17g
"""
    (work / "in.lammps").write_text(script)
    command = [executable, "-in", "in.lammps", "-log", "log.lammps", "-nocite"]
    try:
        proc = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        for name, content in (("stdout.txt", exc.stdout), ("stderr.txt", exc.stderr)):
            if isinstance(content, bytes):
                content = content.decode(errors="replace")
            (work / name).write_text(content or "")
        raise RuntimeError(f"lmp_mpi exceeded {timeout}s; see {work}") from exc
    (work / "stdout.txt").write_text(proc.stdout)
    (work / "stderr.txt").write_text(proc.stderr)
    if proc.returncode:
        raise RuntimeError(f"lmp_mpi failed ({proc.returncode}); see {work}\n{proc.stdout[-2000:]}\n{proc.stderr[-1000:]}")
    if "convergence failed" in proc.stdout.lower():
        raise RuntimeError(f"LAMMPS QEq did not converge; see {work}")
    values = np.loadtxt(work / "energy.txt", ndmin=1)
    rows = (work / "atoms.dump").read_text().splitlines()
    start = next(i for i, line in enumerate(rows) if line.startswith("ITEM: ATOMS"))+1
    atoms = np.loadtxt(rows[start:], ndmin=2)
    version = next((line for line in proc.stdout.splitlines() if line.startswith("LAMMPS (")), "unknown")
    metadata = {"executable": executable, "version": version, "command": command,
                "force_field_sha256": force_field.checksum, "qeq_tolerance": 1e-12}
    (work / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    if version != expected_version:
        raise RuntimeError(f"Reference version mismatch: expected {expected_version!r}, got {version!r}; see {work}")
    if values.shape != (15,) or atoms.shape != (len(symbols), 9) or not np.isfinite(values).all() or not np.isfinite(atoms).all():
        raise RuntimeError(f"Invalid or non-finite reference output; see {work}")
    return ReferenceResult(float(values[0]), atoms[:, 6:9], atoms[:, 2],
                           dict(zip(COMPONENTS, map(float, values[1:]))), version, work)
