"""Primitive-cell energies from tied copies in a finite periodic supercell."""
import autograd.numpy as np
import numpy as onp

from .energy import EnergyModel, BOND_CUT
from .geometry import Boundary


def make_model(ff, symbols, cell=None, pbc=None, max_expanded_atoms=512):
    ff.validate_model(symbols)
    boundary = Boundary(cell, pbc)
    repetitions, shifts, expanded_cell = boundary.supercell(
        ff.general[12], min(BOND_CUT, ff.general[12]), len(symbols), max_expanded_atoms)
    if len(shifts) == 1:
        return EnergyModel(ff, symbols, cell, pbc), tuple(map(int, repetitions))
    return ReplicatedModel(ff, symbols, boundary, shifts, expanded_cell), tuple(map(int, repetitions))


class ReplicatedModel:
    """Copies have distinct interaction identities but tied positions/charges.

    Differentiating the energy per copy with respect to primitive coordinates
    folds all image forces automatically. No enlarged QEq solve is needed.
    """

    def __init__(self, ff, symbols, boundary, shifts, expanded_cell):
        self.n = len(symbols)
        self.copies = len(shifts)
        self.shifts = shifts
        self.mass = onp.array([ff.atoms[s]["mass"] for s in symbols])
        self.model = EnergyModel(ff, tuple(symbols)*self.copies, expanded_cell,
                                 boundary.pbc, charge_atoms=self.n)

    def positions(self, x):
        return np.reshape(x[None, :, :]+self.shifts[:, None, :], (-1, 3))

    def components(self, x, fixed_charges=None):
        charges = None if fixed_charges is None else np.tile(fixed_charges, self.copies)
        values, charges = self.model.components(self.positions(x), charges)
        return values/self.copies, charges[:self.n]

    def properties(self, x, charges):
        result = self.model.properties(self.positions(x), np.tile(charges, self.copies))
        # bond_orders[i,j] sums over all images of primitive atom j. A nonzero
        # diagonal denotes bonding to another image of the same atom.
        result["bond_orders"] = np.mean(np.sum(np.reshape(result["bond_orders"],
            (self.copies, self.n, self.copies, self.n)), axis=2), axis=0)
        result["total_bond_orders"] = np.sum(result["bond_orders"], axis=1)
        result["lone_pairs"] = np.mean(np.reshape(result["lone_pairs"], (self.copies, self.n)), axis=0)
        # Apply the bond threshold per image before folding, never to the sum.
        result["bond_counts"] = result["bond_counts"][:self.n]
        center = np.sum(x*self.mass[:, None], axis=0)/np.sum(self.mass)
        result["dipole"] = np.sum((x-center)*charges[:, None], axis=0)
        return result
