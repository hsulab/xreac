"""Builder-independent neighbor arrays and differentiable edge reductions.

An edge (i, j, S) points from i to the image of j at x[j] + S @ cell.
The caller supplies both directions, including nonzero-shift self images.
No neighbor search or ASE import takes place in this module.
"""
import autograd.numpy as np
from autograd.extend import primitive, defvjp
import numpy as onp

from .geometry import Boundary


@primitive
def scatter_sum(values, indices, size):
    """Sum scalar edge values into atom/pair bins, including repeated indices."""
    return onp.bincount(indices, weights=values, minlength=size).astype(float)


defvjp(scatter_sum, lambda ans, values, indices, size: lambda g: g[indices], None, None)


class Neighbors:
    """Validate and copy an explicit full directed list ``(i, j, S)``."""

    def __init__(self, neighbors, atoms, cell=None, pbc=None):
        boundary = Boundary(cell, pbc)
        self.cell = boundary.cell if boundary.periodic else onp.eye(3)
        self.n = atoms
        if not isinstance(neighbors, (tuple, list)) or len(neighbors) != 3:
            raise ValueError("neighbors must be a tuple (i, j, S)")
        i, j, shifts = (onp.asarray(values) for values in neighbors)
        if i.ndim != 1 or j.shape != i.shape or shifts.shape != (len(i), 3):
            raise ValueError("neighbors require i and j with shape (M,) and S with shape (M, 3)")
        if any(values.dtype.kind not in "iu" for values in (i, j, shifts)):
            raise ValueError("Neighbor indices and lattice shifts must be integers")
        if onp.any(i < 0) or onp.any(j < 0) or onp.any(i >= atoms) or onp.any(j >= atoms):
            raise ValueError("Neighbor atom indices are out of range")
        if onp.any(shifts[:, ~boundary.pbc] != 0):
            raise ValueError("Neighbor shifts must be zero in nonperiodic directions")
        if onp.any((i == j) & onp.all(shifts == 0, axis=1)):
            raise ValueError("Zero-shift self neighbors are not allowed")
        # Builders need not sort by j or shift. Canonical sorting gives stable
        # reductions and a deterministic reverse-edge map for many-body terms.
        order = onp.lexsort((shifts[:, 2], shifts[:, 1], shifts[:, 0], j, i))
        self.i, self.j, self.shifts = (values[order].astype(onp.int64) for values in (i, j, shifts))
        self.offsets = self.shifts @ self.cell
        keys = [(int(a), int(b), *map(int, s)) for a, b, s in zip(self.i, self.j, self.shifts)]
        lookup = {key: k for k, key in enumerate(keys)}
        if len(lookup) != len(keys):
            raise ValueError("Duplicate neighbor edges are not allowed")
        reverse_keys = [(b, a, -sx, -sy, -sz) for a, b, sx, sy, sz in keys]
        if any(key not in lookup for key in reverse_keys):
            raise ValueError("neighbors must include both directions: (i, j, S) and (j, i, -S)")
        self.reverse = onp.array([lookup[key] for key in reverse_keys], dtype=int)
        self.half = onp.array([key < rev for key, rev in zip(keys, reverse_keys)], dtype=bool)
        self.rows = [onp.flatnonzero(self.i == atom) for atom in range(self.n)]
        for values in (self.i, self.j, self.shifts, self.offsets):
            values.flags.writeable = False

    def geometry(self, x):
        vectors = x[self.j] - x[self.i] + self.offsets
        return vectors, np.sqrt(np.sum(vectors*vectors, axis=1))

    def atom_sum(self, values):
        return scatter_sum(values, self.i, self.n)

    def pair_sum(self, values):
        return np.reshape(scatter_sum(values, self.i*self.n+self.j, self.n*self.n), (self.n, self.n))
