"""Image-resolved ASE neighbors and differentiable edge reductions.

An edge (i, j, S) points from i to the image of j at x[j] + S @ cell.
Both directions are retained, including nonzero-shift self images. Topology
is rebuilt at each evaluation and held fixed during coordinate derivatives.
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


class ASENeighbors:
    def __init__(self, positions, cell, pbc, cutoff):
        try:
            from ase.neighborlist import primitive_neighbor_list
        except ImportError as exc:
            raise ImportError("Install xreac[ase] to use neighbor_backend='ase'") from exc
        boundary = Boundary(cell, pbc)
        self.cell = boundary.cell if boundary.periodic else onp.eye(3)
        self.n = len(positions)
        i, j, shifts = primitive_neighbor_list(
            "ijS", boundary.pbc, self.cell, positions,
            onp.nextafter(float(cutoff), onp.inf), self_interaction=False)
        # ASE guarantees only ordering by i. Canonical sorting gives stable
        # reductions and a deterministic reverse-edge map for many-body terms.
        order = onp.lexsort((shifts[:, 2], shifts[:, 1], shifts[:, 0], j, i))
        self.i, self.j, self.shifts = i[order], j[order], shifts[order]
        self.offsets = self.shifts @ self.cell
        keys = [(int(a), int(b), *map(int, s)) for a, b, s in zip(self.i, self.j, self.shifts)]
        lookup = {key: k for k, key in enumerate(keys)}
        reverse_keys = [(b, a, -sx, -sy, -sz) for a, b, sx, sy, sz in keys]
        self.reverse = onp.array([lookup[key] for key in reverse_keys], dtype=int)
        self.half = onp.array([key < rev for key, rev in zip(keys, reverse_keys)], dtype=bool)
        self.rows = [onp.flatnonzero(self.i == atom) for atom in range(self.n)]

    def geometry(self, x):
        vectors = x[self.j] - x[self.i] + self.offsets
        return vectors, np.sqrt(np.sum(vectors*vectors, axis=1))

    def atom_sum(self, values):
        return scatter_sum(values, self.i, self.n)

    def pair_sum(self, values):
        return np.reshape(scatter_sum(values, self.i*self.n+self.j, self.n*self.n), (self.n, self.n))
