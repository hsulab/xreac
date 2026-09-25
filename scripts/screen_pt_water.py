"""Screen near-flat atop water starting orientations before attempting Pt(111) NEB.

Motif: Michaelides et al., PRL 90, 216102 (2003).
O-Pt starting distance: Arnadottir et al., Surface Science 606, 233 (2012).
These are constructed starting geometries, not coordinates extracted from DFT.
"""

import argparse
import json
from pathlib import Path

from ase.io import write
from ase.neighborlist import neighbor_list
from ase.optimize import FIRE

from neb_pt_water import ForceField, ReaxFFCalculator, geometry, np
from xreac.energy import EnergyModel


def closest_pt_h_cutoff(a, ff, metal="Pt"):
    symbols = a.get_chemical_symbols()
    model = EnergyModel(ff, symbols, neighbor_list("ijS", a, ff.general[12]), a.cell.array, a.pbc)
    _, distances = model.edges.geometry(a.positions)
    bo = model.bond_orders(distances)[0]
    indices = [i for i in range(len(bo)) if symbols[model.i[i]] == metal and symbols[model.j[i]] == "H"]
    edge = min(indices, key=lambda i: abs(bo[i] - 0.001))
    return {
        metal.lower(): int(model.i[edge]),
        "h": int(model.j[edge]),
        "bond_order": float(bo[edge]),
        "distance_from_cutoff": float(abs(bo[edge] - 0.001)),
    }


def describe(a):
    v = a.positions[13:15] - a.positions[12]
    n = np.cross(v[0], v[1])
    return dict(
        o_pt_angstrom=a.get_distance(8, 12, mic=True),
        plane_tilt_degrees=float(np.degrees(np.arccos(np.clip(abs(n[2]) / np.linalg.norm(n), 0, 1)))),
        oh_angstrom=np.linalg.norm(v, axis=1).tolist(),
        hoh_degrees=a.get_angle(13, 12, 14),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--ffield", type=Path, help="Explicit experimental parameter file")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    ff = ForceField.from_file(args.ffield) if args.ffield else ForceField.bundled("ffield.reax.PtNiCHO.2016")
    report = dict(
        force_field=ff.path.name,
        sha256=ff.checksum,
        full_derivative=True,
        fmax_target_ev_angstrom=0.02,
        max_steps=args.steps,
        cases=[],
    )
    for azimuth in (0, 30, 60):
        for tilt in (-15, 0, 15):
            a = geometry()
            a.positions[12] = a.positions[8] + [0, 0, 2.34]
            half = np.radians(104.5 / 2)
            v = 0.97 * np.array([[-np.sin(half), np.cos(half), 0], [np.sin(half), np.cos(half), 0]])
            t, z = np.radians([tilt, azimuth])
            rx = np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])
            rz = np.array([[np.cos(z), -np.sin(z), 0], [np.sin(z), np.cos(z), 0], [0, 0, 1]])
            a.positions[13:15] = a.positions[12] + v @ rx.T @ rz.T
            a.calc = ReaxFFCalculator(ff, full_derivative=True)
            label = f"azimuth-{azimuth}-tilt-{tilt}"
            start = describe(a)
            opt = FIRE(a, dt=0.02, maxstep=0.02, logfile=str(args.output / f"{label}.log"))
            history = []
            for _ in opt.irun(fmax=0.02, steps=args.steps):
                history.append(a.positions.copy())
                if len(history) > 30 and np.max(abs(history[-1] - history[-30])) < 1e-10:
                    break
            fmax = float(np.linalg.norm(a.get_forces(), axis=1).max())
            row = dict(
                label=label,
                start=start,
                final=describe(a),
                energy_ev=a.get_potential_energy(),
                fmax_ev_angstrom=fmax,
                converged=fmax <= 0.02,
                steps=opt.nsteps,
                positions=a.positions.tolist(),
                closest_pt_h_cutoff=closest_pt_h_cutoff(a, ff),
            )
            report["cases"].append(row)
            # JSON retains full precision; extxyz is for visualization.
            write(args.output / f"{label}.extxyz", a)
            (args.output / "screen.json").write_text(json.dumps(report, indent=2) + "\n")
            print(label, row["energy_ev"], fmax, row["converged"], flush=True)


if __name__ == "__main__":
    main()
