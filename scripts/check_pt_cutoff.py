"""Probe the stalled Pt/water endpoint on either side of the ReaxFF angle cutoff."""

import argparse
import json
from pathlib import Path

from ase.neighborlist import neighbor_list

from neb_pt_water import Atoms, ForceField, ReaxFFCalculator, comparison, evaluate_lammps, np
from xreac.energy import EnergyModel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primitive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    p = json.loads(args.primitive.read_text())
    ff = ForceField.bundled("ffield.reax.PtNiCHO.2016")
    rows = []
    for dz in (-1e-6, 0, 1e-6):
        a = Atoms(p["symbols"], positions=p["positions"], cell=p["cell"], pbc=p["pbc"])
        a.positions[13, 2] += dz
        a.calc = ReaxFFCalculator(ff)
        energy = a.get_potential_energy()
        actual = a.calc.evaluation
        model = EnergyModel(ff, p["symbols"], neighbor_list("ijS", a, 10), a.cell.array, a.pbc)
        _, distances = model.edges.geometry(a.positions)
        bo = model.bond_orders(distances)[0]
        ix = np.flatnonzero((model.i == 8) & (model.j == 13))
        ix = ix[np.argmin(distances[ix])]
        ref = evaluate_lammps(
            ff, p["symbols"], a.positions, cell=a.cell.array, pbc=a.pbc, directory=args.output / f"offset-{dz}"
        )
        rows.append(
            dict(
                h13_z_offset_angstrom=dz,
                pt8_h13_bond_order=float(bo[ix]),
                energy_ev=energy,
                h13_force_ev_angstrom=a.get_forces()[13].tolist(),
                comparison=comparison(actual, ref, len(a)),
            )
        )
    report = dict(
        source=str(args.primitive), cutoff=0.001, units="comparison errors in kcal/mol, Angstrom, e", rows=rows
    )
    (args.output / "cutoff.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
