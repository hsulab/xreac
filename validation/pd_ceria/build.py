"""Build one small stoichiometric CeO2(111) slab and a tetrahedral Pd4 cluster."""

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell, surface
from ase.constraints import FixAtoms
from ase.io import write


def build():
    a = 5.411  # Angstrom; initial fluorite lattice, not a fitted equilibrium value
    spacing = a / np.sqrt(3)
    slab = surface(bulk("CeO2", "fluorite", a=a), (1, 1, 1), 3)
    # ASE's primitive-cell cut ends with two O planes. Move the outer O plane
    # to the lower periodic image to obtain three complete O-Ce-O trilayers.
    slab.positions[slab.positions[:, 2] > 2.5 * spacing, 2] -= 3 * spacing
    # Twelve primitive surface cells in a compact rectangular box. A literal
    # 3x3 rhombus has a 9.94 A box height, below LAMMPS QEq's 10 A cutoff.
    slab.cell[2, 2] = 3 * spacing
    slab = make_supercell(slab, [[3, 0, 0], [-2, 4, 0], [0, 0, 1]], wrap=False)
    slab.set_tags(np.rint(slab.positions[:, 2] / spacing).astype(int) + 1)
    slab.wrap()
    edge = 2.75
    radius = edge / np.sqrt(3)
    angles = np.arange(3) * 2 * np.pi / 3
    base = np.column_stack((radius * np.cos(angles), radius * np.sin(angles), np.zeros(3)))
    cluster = np.vstack((base, [0, 0, edge * np.sqrt(2 / 3)]))
    cluster += (slab.cell[0] + slab.cell[1]) / 2
    cluster[:, 2] += slab.positions[:, 2].max() + 2.2
    slab += Atoms("Pd4", positions=cluster, tags=[4] * 4)
    slab.center(vacuum=15, axis=2)
    slab.pbc = (True, True, False)
    slab.set_constraint(FixAtoms(mask=slab.get_tags() == 1))
    assert len(slab) == 112 and np.count_nonzero(slab.get_tags() == 1) == 36
    return slab


if __name__ == "__main__":
    atoms = build()
    destination = Path(__file__).with_name("initial.xyz")
    write(destination, atoms, format="extxyz")
    print(f"Wrote {destination}: Ce36O72Pd4; 36 fixed atoms; 15 Angstrom vacuum each side.")
