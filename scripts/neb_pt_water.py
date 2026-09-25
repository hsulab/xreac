"""Exploratory H2O -> OH + H on p(2x2) Pt(111) or Ni(111).

The bottom four metal atoms are fixed. This is a neutral, vacuum, fixed-cell
potential-energy calculation, not an electrochemical or free-energy barrier.
"""

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from ase import Atoms
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps
from water_cluster import comparison


def geometry(product=False, lattice=3.95, layers=3, vacuum=24.0, metal="Pt"):
    slab = fcc111(metal, size=(2, 2, layers), a=lattice, vacuum=vacuum / 2)
    origin = slab.positions[4 * (layers - 1)].copy()
    if product:
        # OH near an atop site; departing H starts at a separate fcc hollow.
        xyz = [origin + [0, 0, 2.0], origin + [-0.6, 0, 2.8]]
        hollow = np.array([2 / 3, 1 / 3]) @ slab.info["adsorbate_info"]["cell"]
        xyz.append(origin + [hollow[0], hollow[1], 1.0])
    else:
        xyz = [origin + [0, 0, 2.3], origin + [-0.76, 0.59, 2.3], origin + [0.76, 0.59, 2.3]]
    slab += Atoms("OH2", positions=xyz)
    slab.set_constraint(FixAtoms(indices=np.flatnonzero(slab.get_tags() == layers)))
    return slab


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "validation/runs/pt-water-gai2016")
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--images", type=int, default=7)
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--metal", choices=("Pt", "Ni"), default="Pt")
    parser.add_argument("--lattice", type=float, help="Defaults: Pt 3.95, Ni 3.52 Angstrom")
    parser.add_argument("--layers", type=int, default=3)
    parser.add_argument(
        "--vacuum", type=float, default=24.0, help="Total padding outside metal slab in Angstrom (half per face)"
    )
    parser.add_argument("--single-point", action="store_true")
    parser.add_argument("--bottom-height", type=float, help="Place the bottom layer at this z height in Angstrom")
    parser.add_argument("--ffield", type=Path, help="Explicit experimental parameter file")
    parser.add_argument("--restart-endpoints", type=Path)
    parser.add_argument("--restart-band", type=Path)
    parser.add_argument("--neb-dtmax", type=float, default=1.0)
    parser.add_argument("--perturbation", type=float, default=0.0)
    parser.add_argument(
        "--allow-unconverged-endpoints",
        action="store_true",
        help="Run an exploratory band, retaining endpoint convergence warnings",
    )
    args = parser.parse_args()
    if args.lattice is None:
        args.lattice = {"Pt": 3.95, "Ni": 3.52}[args.metal]
    if args.layers < 2 or args.vacuum <= 0 or args.images < 3:
        parser.error("Require at least two layers, positive vacuum, and at least three images")
    args.output.mkdir(parents=True, exist_ok=False)
    ff = ForceField.from_file(args.ffield) if args.ffield else ForceField.bundled("ffield.reax.PtNiCHO.2016")
    report = dict(
        force_field=ff.path.name,
        sha256=ff.checksum,
        lattice_angstrom=args.lattice,
        metal=args.metal,
        slab=f"p(2x2) {args.metal}(111), {args.layers} layers, bottom layer fixed",
        layers=args.layers,
        vacuum_padding_angstrom=args.vacuum,
        pbc=[True, True, False],
        image_count=args.images,
        bottom_height_angstrom=args.bottom_height,
        total_charge=0,
        full_derivative=True,
        comparison_units=dict(energy="kcal/mol", forces="kcal/mol/Angstrom", charges="e", dipole="e Angstrom"),
        endpoint_fmax_tolerance_ev_angstrom=0.02,
        max_steps_per_stage=args.steps,
        restart_endpoints=None if args.restart_endpoints is None else str(args.restart_endpoints),
        restart_band=None if args.restart_band is None else str(args.restart_band),
        neb_dtmax=args.neb_dtmax,
        fmax_tolerance_ev_angstrom=args.fmax,
        references={},
        endpoints={},
    )

    def save():
        (args.output / "results.json").write_text(json.dumps(report, indent=2) + "\n")

    def attach(a, full=True):
        a.calc = ReaxFFCalculator(ff, full_derivative=full)
        return a

    def place_slab(a):
        if args.bottom_height is not None:
            a.translate([0, 0, args.bottom_height - a.positions[:4, 2].mean()])
            if a.positions[:, 2].min() < 0 or a.positions[:, 2].max() >= a.cell[2, 2]:
                raise ValueError("Requested slab placement puts atoms outside the nonperiodic cell")
        return a

    def reference(a, label, required=True):
        probe = attach(a.copy(), full=False)
        probe.get_potential_energy()
        actual = probe.calc.evaluation
        try:
            ref = evaluate_lammps(
                ff, a.get_chemical_symbols(), a.positions, cell=a.cell.array, pbc=a.pbc, directory=args.output / label
            )
        except RuntimeError as exc:
            report["references"][label] = dict(passed=False, error=str(exc))
            save()
            print(label, report["references"][label], flush=True)
            if required:
                raise
            return
        result = comparison(actual, ref, len(a))
        result["cell_repetitions"] = ref.cell_repetitions
        report["references"][label] = result
        print(label, result, flush=True)
        save()
        if not result["passed"] and required:
            raise RuntimeError("LAMMPS comparison failed; inspect saved results")

    endpoints = [
        attach(geometry(product=p, lattice=args.lattice, layers=args.layers, vacuum=args.vacuum, metal=args.metal))
        for p in (False, True)
    ]
    if args.restart_endpoints:
        endpoints = [attach(read(args.restart_endpoints / f"{name}.extxyz")) for name in ("reactant", "product")]
        for a in endpoints:
            if a.get_chemical_symbols() != [args.metal] * (4 * args.layers) + ["O", "H", "H"]:
                raise ValueError("Restart endpoints must match the requested layer count and atom ordering")
            expected = geometry(lattice=args.lattice, layers=args.layers, vacuum=args.vacuum, metal=args.metal)
            if not np.allclose(a.cell.array, expected.cell.array) or not np.array_equal(a.pbc, expected.pbc):
                raise ValueError("Restart cell must match --layers, --lattice, --vacuum and slab periodicity")
    if args.perturbation:
        rng = np.random.default_rng(20260925)
        for a in endpoints:
            x = a.positions.copy()
            x[4:] += rng.normal(scale=args.perturbation, size=x[4:].shape)
            a.set_positions(x)
    for a in endpoints:
        place_slab(a)
    report["perturbation_angstrom"] = args.perturbation
    for name, a in zip(("reactant", "product"), endpoints):
        t0 = time.monotonic()
        e = a.get_potential_energy()
        print(name, "initial energy", e, "seconds", time.monotonic() - t0, flush=True)
        write(args.output / f"{name}-initial.extxyz", a)
        reference(a, f"reference-{name}-initial")
        if args.single_point:
            continue
        opt = FIRE(
            a,
            dt=0.02,
            maxstep=0.02,
            logfile=str(args.output / f"{name}.log"),
            trajectory=str(args.output / f"{name}.traj"),
        )
        ok = opt.run(fmax=0.02, steps=args.steps)
        report["endpoints"][name] = dict(
            converged=bool(ok),
            steps=opt.nsteps,
            energy_ev=a.get_potential_energy(),
            fmax_ev_angstrom=float(np.linalg.norm(a.get_forces(), axis=1).max()),
            oh_distances_angstrom=[a.get_distance(len(a) - 3, i, mic=True) for i in (len(a) - 2, len(a) - 1)],
        )
        write(args.output / f"{name}.extxyz", a)
        print(name, report["endpoints"][name], flush=True)
        save()
    if args.single_point:
        return
    endpoints_converged = all(r["converged"] for r in report["endpoints"].values())
    report["endpoints_converged"] = endpoints_converged
    if not endpoints_converged and not args.allow_unconverged_endpoints:
        report["status"] = "Endpoint relaxation incomplete; NEB not started"
        save()
        return
    if max(report["endpoints"]["product"]["oh_distances_angstrom"]) < 1.4:
        report["status"] = "Product recombined to water; no distinct dissociated endpoint"
        save()
        return
    images = [endpoints[0]] + [attach(endpoints[0].copy()) for _ in range(args.images - 2)] + [endpoints[1]]
    if args.restart_band:
        restored = read(args.restart_band, index=":")
        if len(restored) != args.images:
            raise ValueError("Restart band length must match --images")
        for a in restored:
            if (
                a.get_chemical_symbols() != endpoints[0].get_chemical_symbols()
                or not np.allclose(a.cell.array, endpoints[0].cell.array)
                or not np.array_equal(a.pbc, endpoints[0].pbc)
            ):
                raise ValueError("Restart band atoms, cell and periodicity must match the endpoints")
        images = [endpoints[0]] + [attach(place_slab(a)) for a in restored[1:-1]] + [endpoints[1]]
    neb = NEB(images, k=0.1, method="improvedtangent")
    if not args.restart_band:
        neb.interpolate(method="idpp", mic=True)
    report["stages"] = []
    for climb in (False, True):
        neb.climb = climb
        opt = FIRE(
            neb,
            dt=0.03,
            dtmax=args.neb_dtmax,
            maxstep=0.04,
            logfile=str(args.output / f"neb-{climb}.log"),
            trajectory=str(args.output / f"neb-{climb}.traj"),
        )
        started = time.monotonic()
        ok = opt.run(fmax=args.fmax, steps=args.steps)
        report["stages"].append(
            dict(climb=climb, converged=bool(ok), steps=opt.nsteps, wall_seconds=time.monotonic() - started)
        )
        energies = np.array([a.get_potential_energy() for a in images])
        report.update(
            energies_ev=energies.tolist(),
            relative_energies_ev=(energies - energies[0]).tolist(),
            barrier_ev=float(energies.max() - energies[0]),
            reaction_energy_ev=float(energies[-1] - energies[0]),
            neb_fmax_ev_angstrom=float(np.linalg.norm(neb.get_forces(), axis=1).max()),
            status="Converged CI-NEB" if climb and ok and endpoints_converged else "Unconverged reaction path",
        )
        write(args.output / "band.extxyz", images)
        write(args.output / "band.traj", images)
        save()
        print(report["stages"][-1], "barrier eV", report["barrier_ev"], flush=True)
        if not ok:
            break
    for i in sorted({0, int(np.argmax(energies)), len(images) - 1}):
        reference(images[i], f"reference-image-{i}", required=False)


if __name__ == "__main__":
    main()
