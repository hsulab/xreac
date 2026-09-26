"""Build one tetrahedral Pd4 cluster on stoichiometric, unreconstructed rutile (110)."""

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import surface
from ase.constraints import FixAtoms
from ase.io import write
from ase.spacegroup import crystal


def build():
    a, c = 4.594, 2.959
    bulk = crystal(["Ti", "O"], basis=[(0, 0, 0), (0.305, 0.305, 0)], spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
    slab = surface(bulk, (1, 1, 0), layers=3)
    spacing = a / np.sqrt(2)
    # Move the terminal oxygen plane to the opposite face to make three
    # complete, stoichiometric O–TiO–O trilayers with equivalent terminations.
    terminal = slab.positions[:, 2] > 2.5 * spacing
    slab.positions[terminal, 2] -= 3 * spacing
    # Successive (110) repeats also translate laterally by half the surface cell.
    slab.positions[terminal, 0] -= 3 * spacing
    slab.positions[:, 0] %= slab.cell[0, 0]
    slab.set_tags(np.rint(slab.positions[:, 2] / spacing).astype(int) + 1)
    slab = slab.repeat((2, 4, 1))

    edge = 2.75
    base = edge / np.sqrt(3)
    cluster = Atoms(
        "Pd4",
        positions=[
            (-edge / 2, -base / 2, 0),
            (edge / 2, -base / 2, 0),
            (0, base, 0),
            (0, 0, edge * np.sqrt(2 / 3)),
        ],
    )
    cluster.translate([slab.cell[0, 0] / 2, slab.cell[1, 1] / 2, slab.positions[:, 2].max() + 2.3])
    slab += cluster
    slab.center(vacuum=15, axis=2)
    slab.pbc = (True, True, False)
    slab.set_constraint(FixAtoms(mask=slab.get_tags() == 1))
    return slab


if __name__ == "__main__":
    atoms = build()
    write(Path(__file__).with_name("initial.xyz"), atoms, format="extxyz")
    print(f"{atoms.get_chemical_formula()}: {len(atoms)} atoms; {sum(atoms.get_tags() == 1)} fixed")
    print(f"Cell lengths (Å): {atoms.cell.lengths()}")
