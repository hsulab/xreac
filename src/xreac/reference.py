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
from .geometry import Boundary

REFERENCE_VERSION = "LAMMPS (22 Jul 2025 - Update 4)"


@dataclass(frozen=True)
class ReferenceResult:
    energy: float
    forces: np.ndarray
    charges: np.ndarray
    components: dict[str, float]
    version: str
    directory: Path
    total_bond_orders: np.ndarray
    lone_pairs: np.ndarray
    bond_counts: np.ndarray
    dipole: np.ndarray


def evaluate_lammps(force_field, symbols, positions, *, executable=None,
                    directory=None, timeout=120, expected_version=REFERENCE_VERSION,
                    cell=None, pbc=None):
    """Retain a single-point reference with fresh QEq and optional fixed-cell PBC.

    cell/pbc follow Calculator.evaluate(). A rotated triclinic cell is mapped
    to LAMMPS coordinates and vector results are rotated back. Image flags
    preserve the supplied coordinate branch for dipole comparisons.
    """
    symbols, x = validate_input(symbols, positions, cell=cell, pbc=pbc)
    force_field.validate_model(symbols)
    boundary = Boundary(cell, pbc)
    boundary.validate_cutoff(force_field.general[12], min(5.0, force_field.general[12]))
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
    rotation = np.eye(3)
    origin = np.zeros(3)
    images = np.zeros((len(x), 3), dtype=int)
    reference_x = x
    if boundary.periodic:
        # Pad open directions, preserving the periodic lattice vectors. Convert
        # the cell to LAMMPS restricted triclinic coordinates by a rigid rotation.
        reference_cell = boundary.cell.copy()
        fractional = x @ boundary.inverse
        for axis in np.flatnonzero(~boundary.pbc):
            margin = 15 / boundary.heights[axis]
            lo = min(0., fractional[:, axis].min()) - margin
            hi = max(1., fractional[:, axis].max()) + margin
            origin += lo * boundary.cell[axis]
            reference_cell[axis] *= hi - lo
        q, r = np.linalg.qr(reference_cell.T)
        rotation = q @ np.diag(np.sign(np.diag(r)))
        restricted = reference_cell @ rotation
        fractional = (x - origin) @ np.linalg.inv(reference_cell)
        images[:, boundary.pbc] = np.floor(fractional[:, boundary.pbc]).astype(int)
        reference_x = (fractional - images) @ restricted
        for axis, dim in enumerate("xyz"):
            lines.append(f"0 {restricted[axis, axis]:.17g} {dim}lo {dim}hi")
        lines.append(f"{restricted[1, 0]:.17g} {restricted[2, 0]:.17g} {restricted[2, 1]:.17g} xy xz yz")
    else:
        for axis, dim in enumerate("xyz"):
            lines.append(f"{x[:, axis].min()-15:.17g} {x[:, axis].max()+15:.17g} {dim}lo {dim}hi")
    lines += ["", "Masses", ""]
    lines += [f"{i+1} {force_field.atoms[s]['mass']:.17g}" for i, s in enumerate(types)]
    lines += ["", "Atoms # charge", ""]
    lines += [f"{i+1} {types.index(s)+1} 0 " + " ".join(f"{v:.17g}" for v in pos)
              + " " + " ".join(str(v) for v in images[i])
              for i, (s, pos) in enumerate(zip(symbols, reference_x))]
    (work / "atoms.data").write_text("\n".join(lines)+"\n")
    terms = " ".join(f"$(c_reax[{i}]:%.17g)" for i in range(1, 15))
    script = f"""units real
atom_style charge
boundary {' '.join('p' if flag else 'f' for flag in boundary.pbc)}
read_data atoms.data
pair_style reaxff NULL tabulate 0 enobonds yes
pair_coeff * * ffield {' '.join(types)}
fix charges all qeq/reaxff 1 0 {force_field.general[12]:.17g} 1e-12 reaxff maxiter 2000
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
compute reax all pair reaxff
compute bondinfo all reaxff/atom
compute dipole all dipole
dump atoms all custom 1 atoms.dump id type q x y z fx fy fz c_bondinfo[1] c_bondinfo[2] c_bondinfo[3]
dump_modify atoms sort id format float %.17g
thermo_style custom step pe c_reax[*] c_dipole[1] c_dipole[2] c_dipole[3]
thermo_modify format float %.17g
run 0
print "$(pe:%.17g) {terms}" file energy.txt screen no
print "$(c_dipole[1]:%.17g) $(c_dipole[2]:%.17g) $(c_dipole[3]:%.17g)" file dipole.txt screen no
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
    dipole = np.loadtxt(work / "dipole.txt", ndmin=1)
    version = next((line for line in proc.stdout.splitlines() if line.startswith("LAMMPS (")), "unknown")
    metadata = {"executable": executable, "version": version, "command": command,
                "force_field_sha256": force_field.checksum, "qeq_tolerance": 1e-12,
                "cell": boundary.cell.tolist() if boundary.periodic else None,
                "pbc": boundary.pbc.tolist(), "rotation": rotation.tolist(),
                "origin": origin.tolist()}
    (work / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    if version != expected_version:
        raise RuntimeError(f"Reference version mismatch: expected {expected_version!r}, got {version!r}; see {work}")
    if values.shape != (15,) or atoms.shape != (len(symbols), 12) or dipole.shape != (3,) or not all(np.isfinite(v).all() for v in (values, atoms, dipole)):
        raise RuntimeError(f"Invalid or non-finite reference output; see {work}")
    return ReferenceResult(float(values[0]), atoms[:, 6:9] @ rotation.T, atoms[:, 2],
                           dict(zip(COMPONENTS, map(float, values[1:]))), version, work,
                           atoms[:, 9], atoms[:, 10], atoms[:, 11].astype(int), dipole @ rotation.T)
