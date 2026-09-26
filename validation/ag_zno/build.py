"""Build one nonpolar ZnO(10-10) slab and one supported tetrahedral Ag4."""

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, surface
from ase.constraints import FixAtoms
from ase.io import write


def build():
    here = Path(__file__).resolve().parent
    # Three-index (100) in the hexagonal primitive cell is (10-10).
    cell = bulk("ZnO", "wurtzite", a=3.25, c=5.21, u=0.382)
    slab = surface(cell, (1, 0, 0), layers=4).repeat((4, 3, 1))
    z = slab.positions[:, 2]
    slab.set_tags(np.floor((z - z.min() + 1e-6) / (3.25 * np.sqrt(3) / 2)).astype(int) + 1)
    slab.set_constraint(FixAtoms(mask=slab.get_tags() == 1))
    slab.center(vacuum=15, axis=2)
    write(here / "clean.xyz", slab)

    edge = 2.85
    base = edge / np.sqrt(3)
    cluster = Atoms(
        "Ag4",
        positions=[(-edge / 2, -base / 2, 0), (edge / 2, -base / 2, 0), (0, base, 0), (0, 0, edge * np.sqrt(2 / 3))],
    )
    # Anchor one base atom above a central top-surface O atom: minimum gap 2.5 Å.
    candidates = [a for a in slab if a.symbol == "O" and abs(a.z - slab.positions[:, 2].max()) < 1e-5]
    center = slab.cell.lengths()[:2] / 2
    anchor = min(candidates, key=lambda a: np.linalg.norm(a.position[:2] - center))
    cluster.translate(anchor.position + [0, 0, 2.5] - cluster.positions[0])
    slab += cluster
    slab.center(vacuum=15, axis=2)
    slab.set_constraint(FixAtoms(mask=slab.get_tags() == 1))
    write(here / "initial.xyz", slab)
    return slab


if __name__ == "__main__":
    atoms = build()
    print(atoms.get_chemical_formula(), atoms.cell.lengths())
