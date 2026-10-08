"""Minimise 1, 2, 4, or 8 H2O on four-layer p(1x4) anatase (101)."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.build import surface
from ase.constraints import FixAtoms
from ase.io import write
from ase.optimize import FIRE
from ase.spacegroup import crystal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps
from water_cluster import comparison

FORCE_FIELD = "ffield.reax.TiOH.Monti2012"
TI_O_ADSORPTION_LIMIT = 2.6  # Broad structural check, not a fitted equilibrium distance.


def adsorption_status(atoms, waters):
    """Distinguish optimiser convergence from intact water bound to surface Ti.

    Distances above 2.6 A do not reproduce the Ti-bound molecular adsorption
    state. This broad check never restrains the geometry.
    """
    ti = np.flatnonzero(atoms.numbers[:192] == 22)
    rows = []
    for n in range(waters):
        oxygen = 192 + 3 * n
        distances = atoms.get_distances(oxygen, ti, mic=True)
        nearest = int(np.argmin(distances))
        oh = atoms.get_distances(oxygen, [oxygen + 1, oxygen + 2], mic=True)
        rows.append(
            dict(
                oh_distances_angstrom=oh.tolist(),
                nearest_ti_index=int(ti[nearest]),
                nearest_ti_o_angstrom=float(distances[nearest]),
                molecular=bool(np.all((oh > 0.7) & (oh < 1.25))),
                ti_bound=bool(distances[nearest] <= TI_O_ADSORPTION_LIMIT),
            )
        )
    return dict(
        passed=all(row["molecular"] and row["ti_bound"] for row in rows),
        ti_o_limit_angstrom=TI_O_ADSORPTION_LIMIT,
        water_geometry=rows,
    )


def geometry(waters=1):
    """192 slab atoms: four 12-atom repeats, then 1x4 in-plane replication.

    Rectangular surface axes are [10-1] and [010]. The bottom stoichiometric
    repeat (48 atoms) is fixed; tag zero identifies the adsorbates.
    """
    if waters not in (1, 2, 4, 8):
        raise ValueError("waters must be 1, 2, 4, or 8")
    a, c, u = 3.784, 9.515, 0.208
    bulk = crystal(
        ["Ti", "O"],
        basis=[(0, 0, 0), (0, 0, u)],
        spacegroup=141,
        setting=1,
        cellpar=[a, a, c, 90, 90, 90],
    )
    slab = surface(bulk, (1, 0, 1), layers=4)
    spacing = a * c / np.hypot(a, c)
    lateral_shift = a * a / np.hypot(a, c)
    # Move the terminal O plane to the bottom, including the lateral stacking
    # shift. This exposes Ti5c/O2c on both faces instead of Ti4c/O1c.
    terminal = slab.positions[:, 2] > 3.5 * spacing
    slab.positions[terminal] -= [4 * lateral_shift, 0, 4 * spacing]
    slab.set_tags(np.floor((slab.positions[:, 2] + u * spacing + 1e-6) / spacing).astype(int) + 1)
    slab = slab.repeat((1, 4, 1))
    slab.center(vacuum=12, axis=2)
    slab.wrap()
    ti = np.flatnonzero(slab.numbers == 22)
    sites = ti[np.isclose(slab.positions[ti, 2], slab.positions[ti, 2].max())]
    # Equivalent neighbouring Ti5c sites along [010], 3.784 A apart.
    site = sites[np.lexsort((slab.positions[sites, 1], slab.positions[sites, 0]))[0]]
    origin = slab.positions[site].copy() + [0, a, 2.3]
    other_row = sites[~np.isclose(slab.positions[sites, 0], slab.positions[site, 0])]
    second_site = other_row[np.argmin(slab.positions[other_row, 1])]
    second_origin = slab.positions[second_site].copy() + [0, a, 2.3]
    bond, half_angle = 0.9572, np.deg2rad(104.52 / 2)
    for n in range(waters):
        # Four waters fill one Ti5c row; eight fill both rows, without
        # duplicating sites through the periodic y boundary.
        oxygen = (origin if n < 4 else second_origin) + [0, (n % 4) * a, 0]
        offsets = np.array(
            [
                [0, 0, 0],
                [bond * np.cos(half_angle), bond * np.sin(half_angle), 0],
                [bond * np.cos(half_angle), -bond * np.sin(half_angle), 0],
            ]
        )
        slab += Atoms("OH2", positions=oxygen + offsets)
    slab.set_constraint(FixAtoms(mask=slab.get_tags() == 1))
    slab.wrap()
    return slab


def verify(atoms, ff, directory):
    probe = atoms.copy()
    probe.calc = ReaxFFCalculator(ff)
    probe.get_potential_energy()
    reference = evaluate_lammps(
        ff,
        probe.get_chemical_symbols(),
        probe.positions,
        cell=probe.cell.array,
        pbc=probe.pbc,
        directory=directory,
    )
    result = comparison(probe.calc.evaluation, reference, len(probe))
    if not result["passed"]:
        raise RuntimeError(f"LAMMPS comparison failed: {result}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--waters", type=int, choices=(1, 2, 4, 8), nargs="+", default=[1, 2])
    parser.add_argument("--ffield", type=Path)
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--fmax", type=float, default=0.02, help="eV/Angstrom on mobile atoms")
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--verify", action="store_true", help="Compare initial/final fixed-charge results with LAMMPS")
    parser.add_argument(
        "--full-derivative", action="store_true", help="Diagnostic charge-response forces; default matches LAMMPS"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/runs" / ("anatase-water-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")),
    )
    args = parser.parse_args()
    if args.steps < 1 or not np.isfinite(args.fmax) or args.fmax <= 0:
        parser.error("steps and fmax must be positive")
    ff = ForceField.from_file(args.ffield) if args.ffield else ForceField.bundled(FORCE_FIELD)
    args.output.mkdir(parents=True, exist_ok=False)
    passed = True
    for waters in dict.fromkeys(args.waters):
        destination = args.output / f"{waters}water"
        destination.mkdir()
        atoms = geometry(waters)
        write(destination / "initial.xyz", atoms, format="extxyz")
        if args.build_only:
            continue
        atoms.calc = ReaxFFCalculator(ff, full_derivative=args.full_derivative)
        report = dict(
            waters=waters,
            coverage_monolayers=waters / 8,
            slab_atoms=192,
            fixed_atoms=48,
            layers=4,
            repeat=[1, 4],
            cell_angstrom=atoms.cell.array.tolist(),
            pbc=atoms.pbc.tolist(),
            lattice_angstrom=dict(a=3.784, c=9.515),
            oxygen_internal_coordinate=0.208,
            force_field=ff.path.name,
            force_field_sha256=ff.checksum,
            full_derivative=args.full_derivative,
            force_convention="charge_response" if args.full_derivative else "fixed_charge",
            total_charge=0,
            fmax_ev_angstrom=args.fmax,
            initial_energy_ev=atoms.get_potential_energy(),
        )
        if args.verify:
            report["initial_reference"] = verify(atoms, ff, destination / "initial_lammps")
        with FIRE(
            atoms,
            dt=0.05,
            dtmax=0.3,
            maxstep=0.1,
            logfile=str(destination / "optimize.log"),
            trajectory=str(destination / "relaxation.traj"),
        ) as optimizer:
            converged = optimizer.run(fmax=args.fmax, steps=args.steps)
            report["steps"] = optimizer.nsteps
        report.update(
            converged=bool(converged),
            energy_ev=atoms.get_potential_energy(),
            max_mobile_force_ev_angstrom=float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
            charges_e=atoms.get_charges().tolist(),
        )
        adsorption = adsorption_status(atoms, waters)
        report["water_geometry"] = adsorption.pop("water_geometry")
        report["adsorption"] = adsorption
        report["passed"] = bool(converged and adsorption["passed"])
        write(destination / "relaxed.xyz", atoms, format="extxyz")
        if args.verify:
            report["final_reference"] = verify(atoms, ff, destination / "final_lammps")
        (destination / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        print(
            f"{waters} water: converged={converged}, Ti-bound molecular adsorption={adsorption['passed']}, "
            f"E={report['energy_ev']:.8f} eV",
            flush=True,
        )
        passed &= report["passed"]
    print(f"Results: {args.output}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
