"""Quick unconstrained ASE NEB of charged SN2, initialized from scan_sn2.py output."""

import argparse
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write
from ase.mep import NEB
from ase.optimize import BFGS, FIRE

from scan_sn2 import CHOCL_FORCE_FIELD, ROOT, SYMBOLS, ForceField
from xreac.ase import ReaxFFCalculator

EV_TO_KJ_MOL = 96.4853321233


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", type=Path, default=ROOT / "validation/cho/sn2_scan.json")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/runs/sn2-neb")
    parser.add_argument("--images", type=int, default=9)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--perturbation", type=float, default=0.015)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    scan = json.loads(args.scan.read_text())
    ff = ForceField.bundled(CHOCL_FORCE_FIELD)

    def atoms(x):
        a = Atoms(SYMBOLS, positions=x)
        a.calc = ReaxFFCalculator(ff, total_charge=-1, full_derivative=True)
        return a

    original = np.array([r["positions"] for r in scan["continuous"]])
    start = atoms(original[0])
    opt = BFGS(start, logfile=str(args.output / "endpoint.log"), maxstep=0.03)
    endpoint_converged = opt.run(fmax=0.01, steps=250)
    # Symmetry-related product retains atom labels, exchanging the two chlorines.
    finish = start.positions.copy()
    finish[:, 2] *= -1
    finish[[4, 5]] = finish[[5, 4]]
    images = []
    rng = np.random.default_rng(20260925)
    for t in np.linspace(0, 1, args.images):
        j = t * (len(original) - 1)
        lo, hi = int(np.floor(j)), int(np.ceil(j))
        x = (1 - (j - lo)) * original[lo] + (j - lo) * original[hi]
        x += (1 - t) * (start.positions - original[0]) + t * (finish - original[-1])
        # Break exact collinearity and methyl symmetry so NEB can explore bending.
        if 0 < t < 1:
            x += args.perturbation * np.sin(np.pi * t) * rng.normal(size=x.shape)
        images.append(atoms(x))
    neb = NEB(images, method="improvedtangent", k=0.1, remove_rotation_and_translation=True)
    initial_energies = [a.get_potential_energy() for a in images]
    write(args.output / "initial.extxyz", images)
    history = []

    def record():
        energies = np.array([a.get_potential_energy() for a in images])
        forces = neb.get_forces()
        history.append(
            dict(
                climb=neb.climb,
                barrier_kj_mol=float((energies.max() - energies[0]) * EV_TO_KJ_MOL),
                fmax=float(np.linalg.norm(forces, axis=1).max()),
            )
        )
        if len(history) % 25 == 0:
            print(history[-1], flush=True)

    stages = []
    for climb in (False, True):
        neb.climb = climb
        opt = FIRE(
            neb,
            dt=0.03,
            maxstep=0.04,
            logfile=str(args.output / f"neb-{climb}.log"),
            trajectory=str(args.output / f"neb-{climb}.traj"),
        )
        opt.attach(record, interval=1)
        converged = opt.run(fmax=args.fmax, steps=args.steps)
        stages.append(dict(climb=climb, converged=bool(converged), steps=opt.nsteps))
        write(args.output / f"band-{climb}.extxyz", images)
        if not converged:
            # An unresolved regular band is not a sound starting point for
            # climbing-image optimization, especially if all images lie below
            # the endpoints because an intervening barrier was undersampled.
            break
    energies = np.array([a.get_potential_energy() for a in images])
    distances = [[float(a.get_distance(0, i)) for i in (4, 5)] for a in images]
    dense = []
    probe = atoms(images[0].positions)
    for i in range(len(images) - 1):
        for fraction in np.linspace(0, 1, 21)[:-1]:
            probe.positions = (1 - fraction) * images[i].positions + fraction * images[i + 1].positions
            dense.append(
                dict(
                    image_coordinate=float(i + fraction),
                    energy_kj_mol=float((probe.get_potential_energy() - energies[0]) * EV_TO_KJ_MOL),
                )
            )
    report = dict(
        force_field=CHOCL_FORCE_FIELD,
        checksum=ff.checksum,
        total_charge=-1,
        full_derivative=True,
        fmax_tolerance=args.fmax,
        perturbation=args.perturbation,
        converged=bool(stages[-1]["climb"] and stages[-1]["converged"]),
        final_neb_fmax=history[-1]["fmax"],
        dense_path=dense,
        dense_path_barrier_kj_mol=max(r["energy_kj_mol"] for r in dense),
        endpoint_converged=bool(endpoint_converged),
        endpoint_fmax=float(np.linalg.norm(start.get_forces(), axis=1).max()),
        stages=stages,
        history=history,
        energies_ev=energies.tolist(),
        relative_energies_kj_mol=((energies - energies[0]) * EV_TO_KJ_MOL).tolist(),
        barrier_kj_mol=float((energies.max() - energies[0]) * EV_TO_KJ_MOL),
        initial_barrier_kj_mol=float((max(initial_energies) - initial_energies[0]) * EV_TO_KJ_MOL),
        scan_barrier_kj_mol=scan["continuous_barrier"] * 4.184,
        endpoint_energy_ev=float(energies[0]),
        c_cl_distances=distances,
        positions=[a.positions.tolist() for a in images],
        symbols=list(SYMBOLS),
        charges=[a.get_charges().tolist() for a in images],
    )
    (args.output / "neb.json").write_text(json.dumps(report, indent=2) + "\n")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
    ax.plot(range(args.images), report["relative_energies_kj_mol"], "o", label="NEB images")
    ax.plot(
        [r["image_coordinate"] for r in dense], [r["energy_kj_mol"] for r in dense], label="Interpolated path check"
    )
    ax.axhline(report["scan_barrier_kj_mol"], ls="--", color="gray", label="Previous constrained scan barrier")
    ax.set(
        xlabel="Image index",
        ylabel="Potential energy above reactant (kJ/mol)",
        title=f"Charged SN2: barrier {report['barrier_kj_mol']:.1f} kJ/mol",
    )
    ax.legend(fontsize=8)
    fig.savefig(args.output / "neb.png", dpi=180)
    print(
        json.dumps({k: report[k] for k in ("barrier_kj_mol", "endpoint_converged", "endpoint_fmax", "stages")}),
        flush=True,
    )


if __name__ == "__main__":
    main()
