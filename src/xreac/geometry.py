"""Fixed-cell periodic geometry with explicit finite-cutoff image sums."""
from itertools import product

import autograd.numpy as np
from autograd.tracer import getval
import numpy as onp


def validate_expansion_limit(value):
    if isinstance(value, bool) or not isinstance(value, (int, onp.integer)) or value < 1:
        raise ValueError("max_expanded_atoms must be a positive integer")


class Boundary:
    """Cell vectors are rows, as in ASE; pbc selects periodic lattice vectors."""

    def __init__(self, cell=None, pbc=None):
        if pbc is None:
            pbc = cell is not None
        flags = onp.asarray(pbc)
        if flags.shape not in ((), (3,)) or flags.dtype.kind not in "biu" or not onp.isin(flags, [0, 1]).all():
            raise ValueError("pbc must be a boolean or three boolean flags")
        self.pbc = onp.broadcast_to(flags.astype(bool), (3,)).copy()
        self.periodic = bool(self.pbc.any())
        self.cell = None
        if not self.periodic:
            return
        if cell is None:
            raise ValueError("A cell is required when pbc is enabled")
        matrix = onp.array(cell, dtype=float, copy=True)
        if matrix.shape == (3,):
            if onp.any(matrix <= 0):
                raise ValueError("Cell lengths must be positive")
            matrix = onp.diag(matrix)
        if matrix.shape != (3, 3) or not onp.isfinite(matrix).all():
            raise ValueError("cell must contain three finite lengths or a finite (3, 3) matrix")
        determinant = onp.linalg.det(matrix)
        if determinant <= 1e-12*onp.prod(onp.linalg.norm(matrix, axis=1)):
            raise ValueError("Periodic cells must be nonsingular and right-handed")
        self.cell = matrix
        self.inverse = onp.linalg.inv(matrix)
        self.heights = 1/onp.linalg.norm(self.inverse, axis=0)
        # With the cutoff/height restriction below, these shifts include every
        # contributing image after centering fractional pair displacements.
        self.translations = onp.array(list(product(*[(-1, 0, 1) if p else (0,) for p in self.pbc]))) @ matrix

    def validate_cutoff(self, cutoff, bond_cutoff):
        if self.periodic:
            minimum = max(cutoff, 2*bond_cutoff)
            if onp.any(self.heights[self.pbc] <= minimum*(1+1e-12)):
                raise ValueError(f"Periodic cell heights must exceed {minimum:g} Angstrom "
                                 "(nonbonded cutoff and twice the bond cutoff); use a larger supercell")

    def supercell(self, cutoff, bond_cutoff, atoms, max_expanded_atoms=512):
        """Return repetitions, translation vectors, and a safe internal cell.

        Replication gives every interacting image a distinct atom identity.
        The public coordinates remain the primitive degrees of freedom.
        """
        validate_expansion_limit(max_expanded_atoms)
        repetitions = onp.ones(3, dtype=int)
        if self.periodic:
            minimum = max(cutoff, 2*bond_cutoff)*(1+1e-12)
            required = onp.floor(minimum/self.heights[self.pbc])+1
            copies = float(onp.prod(required))
            if copies > 1 and copies*atoms > max_expanded_atoms:
                raise ValueError(f"Small-cell replication needs {copies*atoms:g} internal atoms, "
                                 f"exceeding max_expanded_atoms={max_expanded_atoms}; "
                                 "increase this limit explicitly if memory permits")
            repetitions[self.pbc] = required.astype(int)
            shifts = onp.array(list(product(*(range(n) for n in repetitions)))) @ self.cell
            return repetitions, shifts, self.cell*repetitions[:, None]
        return repetitions, onp.zeros((1, 3)), None

    def centered_displacements(self, x):
        delta = x[:, None, :] - x[None, :, :]
        if not self.periodic:
            return delta
        # Image choices are discrete; differentiate only within the selected
        # branch, keeping lattice shifts constant during the derivative.
        images = onp.rint(onp.asarray(getval(delta)) @ self.inverse) * self.pbc
        return delta - images @ self.cell

    def minimum_displacements(self, x):
        delta = self.centered_displacements(x)
        if not self.periodic:
            return delta
        candidates = getval(delta)[None, ...] + self.translations[:, None, None, :]
        closest = onp.argmin(onp.sum(candidates**2, axis=-1), axis=0)
        return delta + self.translations[closest]

    def image_displacements(self, x):
        return self.centered_displacements(x)[None, ...] + self.translations[:, None, None, :]
