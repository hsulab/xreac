"""Development-only single-rank lmp_mpi reference harness."""

from dataclasses import dataclass, replace
from pathlib import Path
import json
import os
import shutil
import subprocess
import tempfile

import numpy as np

from .calculator import validate_input
from .energy import COMPONENTS
from .geometry import Boundary, validate_expansion_limit

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
    cell_repetitions: tuple[int, int, int] = (1, 1, 1)


def evaluate_lammps(
    force_field,
    symbols,
    positions,
    *,
    executable=None,
    directory=None,
    timeout=120,
    expected_version=REFERENCE_VERSION,
    cell=None,
    pbc=None,
    allow_small_cell=False,
    max_expanded_atoms=512,
    supplied_charges=None,
):
    """Retain a single-point reference with fresh QEq and optional fixed-cell PBC.

    Small cells are replicated and results normalized to the input cell.
    cell/pbc follow Calculator.evaluate(). A rotated triclinic cell is mapped
    to LAMMPS coordinates and vector results are rotated back. Image flags
    preserve the supplied coordinate branch for dipole comparisons.
    allow_small_cell=True is for diagnostics outside LAMMPS's documented QEq
    cell-size range, not an accepted reference for validating small-cell support.
    supplied_charges is an optional finite (N,) array in e. When provided,
    disable LAMMPS QEq and use these charges unchanged, including nonzero net
    charge. This validates energies and fixed-charge forces, not the QEq solve.
    """
    symbols, x = validate_input(symbols, positions, cell=cell, pbc=pbc)
    force_field.validate_model(symbols)
    if supplied_charges is not None:
        supplied_charges = np.array(supplied_charges, dtype=float, copy=True)
        if supplied_charges.shape != (len(symbols),) or not np.isfinite(supplied_charges).all():
            raise ValueError("supplied_charges must be a finite (N,) array")
    boundary = Boundary(cell, pbc)
    if not isinstance(allow_small_cell, bool):
        raise ValueError("allow_small_cell must be a boolean")
    validate_expansion_limit(max_expanded_atoms)
    within_cell_limits = True
    try:
        boundary.validate_cutoff(force_field.general[12], min(5.0, force_field.general[12]))
    except ValueError:
        within_cell_limits = False
        if not allow_small_cell:
            repetitions, shifts, expanded_cell = boundary.supercell(
                force_field.general[12], min(5.0, force_field.general[12]), len(x), max_expanded_atoms
            )
            copies, n = len(shifts), len(x)
            ref = evaluate_lammps(
                force_field,
                symbols * copies,
                np.reshape(x[None, :, :] + shifts[:, None, :], (-1, 3)),
                cell=expanded_cell,
                pbc=boundary.pbc,
                executable=executable,
                directory=directory,
                timeout=timeout,
                expected_version=expected_version,
                max_expanded_atoms=max_expanded_atoms,
                supplied_charges=None if supplied_charges is None else np.tile(supplied_charges, copies),
            )
            folded = {}
            for key, tolerance in (
                ("forces", 1e-4),
                ("charges", 1e-6),
                ("total_bond_orders", 1e-8),
                ("lone_pairs", 1e-8),
                ("bond_counts", 0),
            ):
                values = getattr(ref, key).reshape((copies, n) + getattr(ref, key).shape[1:])
                if np.max(abs(values - values[0])) > tolerance:
                    raise RuntimeError(f"LAMMPS replicated {key} differ between copies; see {ref.directory}")
                folded[key] = values.mean(axis=0) if key != "bond_counts" else values[0].copy()
            metadata_path = ref.directory / "metadata.json"
            metadata = json.loads(metadata_path.read_text())
            metadata.update(
                input_cell=boundary.cell.tolist(),
                input_atoms=n,
                cell_repetitions=repetitions.tolist(),
                energy_divisor=copies,
                returned_results="Per input cell; forces/properties averaged across equivalent copies",
                input_total_charge=0.0 if supplied_charges is None else float(supplied_charges.sum()),
            )
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            (ref.directory / "primitive.json").write_text(
                json.dumps(
                    dict(
                        symbols=symbols,
                        positions=x.tolist(),
                        cell=boundary.cell.tolist(),
                        pbc=boundary.pbc.tolist(),
                        supplied_charges=None if supplied_charges is None else supplied_charges.tolist(),
                    ),
                    indent=2,
                )
                + "\n"
            )
            return replace(
                ref,
                energy=ref.energy / copies,
                components={key: value / copies for key, value in ref.components.items()},
                dipole=np.sum(
                    (x - np.average(x, axis=0, weights=[force_field.atoms[s]["mass"] for s in symbols]))
                    * folded["charges"][:, None],
                    axis=0,
                ),
                cell_repetitions=tuple(map(int, repetitions)),
                **folded,
            )
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
            lo = min(0.0, fractional[:, axis].min()) - margin
            hi = max(1.0, fractional[:, axis].max()) + margin
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
            lines.append(f"{x[:, axis].min() - 15:.17g} {x[:, axis].max() + 15:.17g} {dim}lo {dim}hi")
    lines += ["", "Masses", ""]
    lines += [f"{i + 1} {force_field.atoms[s]['mass']:.17g}" for i, s in enumerate(types)]
    lines += ["", "Atoms # charge", ""]
    initial_charges = np.zeros(len(x)) if supplied_charges is None else supplied_charges
    lines += [
        f"{i + 1} {types.index(s) + 1} {initial_charges[i]:.17g} "
        + " ".join(f"{v:.17g}" for v in pos)
        + " "
        + " ".join(str(v) for v in images[i])
        for i, (s, pos) in enumerate(zip(symbols, reference_x))
    ]
    (work / "atoms.data").write_text("\n".join(lines) + "\n")
    terms = " ".join(f"$(c_reax[{i}]:%.17g)" for i in range(1, 15))
    qeq = (
        f"fix charges all qeq/reaxff 1 0 {force_field.general[12]:.17g} 1e-12 reaxff maxiter 2000"
        if supplied_charges is None
        else "# Supplied charges; no QEq"
    )
    checkqeq = "" if supplied_charges is None else " checkqeq no"
    script = f"""units real
atom_style charge
boundary {" ".join("p" if flag else "f" for flag in boundary.pbc)}
read_data atoms.data
pair_style reaxff NULL tabulate 0 enobonds yes{checkqeq}
pair_coeff * * ffield {" ".join(types)}
{qeq}
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
        raise RuntimeError(
            f"lmp_mpi failed ({proc.returncode}); see {work}\n{proc.stdout[-2000:]}\n{proc.stderr[-1000:]}"
        )
    if "convergence failed" in proc.stdout.lower():
        raise RuntimeError(f"LAMMPS QEq did not converge; see {work}")
    values = np.loadtxt(work / "energy.txt", ndmin=1)
    rows = (work / "atoms.dump").read_text().splitlines()
    start = next(i for i, line in enumerate(rows) if line.startswith("ITEM: ATOMS")) + 1
    atoms = np.loadtxt(rows[start:], ndmin=2)
    dipole = np.loadtxt(work / "dipole.txt", ndmin=1)
    version = next((line for line in proc.stdout.splitlines() if line.startswith("LAMMPS (")), "unknown")
    metadata = {
        "executable": executable,
        "version": version,
        "command": command,
        "force_field_sha256": force_field.checksum,
        "qeq_tolerance": 1e-12 if supplied_charges is None else None,
        "charge_mode": "equilibrated" if supplied_charges is None else "supplied",
        "total_charge": float(initial_charges.sum()),
        "cell": boundary.cell.tolist() if boundary.periodic else None,
        "pbc": boundary.pbc.tolist(),
        "rotation": rotation.tolist(),
        "origin": origin.tolist(),
        "allow_small_cell": allow_small_cell,
        "within_validated_cell_limits": within_cell_limits,
    }
    (work / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if version != expected_version:
        raise RuntimeError(f"Reference version mismatch: expected {expected_version!r}, got {version!r}; see {work}")
    if (
        values.shape != (15,)
        or atoms.shape != (len(symbols), 12)
        or dipole.shape != (3,)
        or not all(np.isfinite(v).all() for v in (values, atoms, dipole))
    ):
        raise RuntimeError(f"Invalid or non-finite reference output; see {work}")
    if supplied_charges is not None and not np.allclose(atoms[:, 2], supplied_charges, atol=1e-12, rtol=0):
        raise RuntimeError(f"LAMMPS changed supplied charges; see {work}")
    return ReferenceResult(
        float(values[0]),
        atoms[:, 6:9] @ rotation.T,
        atoms[:, 2],
        dict(zip(COMPONENTS, map(float, values[1:]))),
        version,
        work,
        atoms[:, 9],
        atoms[:, 10],
        atoms[:, 11].astype(int),
        dipole @ rotation.T,
    )
