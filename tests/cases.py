"""Deterministic validation geometries, including active O-only many-body terms."""

import numpy as np


def cluster(n, spacing=2.15):
    side = int(np.ceil(n ** (1 / 3)))
    ijk = np.array(list(np.ndindex(side, side, side)))[:n]
    rng = np.random.default_rng(1729 + n)
    symbols = ["Zn" if row.sum() % 2 == 0 else "O" for row in ijk]
    return symbols, ijk * spacing + rng.normal(0, 0.025, (n, 3))


CASES = {
    "zn_atom": (["Zn"], [[0.0, 0.0, 0.0]]),
    "zno": (["Zn", "O"], [[0, 0, 0], [1.9, 0, 0]]),
    "o4": (["O"] * 4, [[-1.1, 0.6, 0.2], [0, 0, 0], [1.3, 0.1, 0], [1.8, 1.2, 0.5]]),
    "cluster20": cluster(20),
}
