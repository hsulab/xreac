"""Validate the isolated cache experiment without changing the core solver."""

import numpy as np
import pytest

from benchmark_water_qeq import CachedQEq, PairedQEqProbe


def system(size, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(size, size))
    h = a.T @ a + 2 * np.eye(size)
    kkt = np.block([[h, np.ones((size, 1))], [np.ones((1, size)), np.zeros((1, 1))]])
    rhs = np.append(rng.normal(size=size), 0)
    return kkt, rhs


def test_cache_refines_current_matrix_and_rhs():
    matrix, rhs = system(6, 42)
    probe = PairedQEqProbe(np.linalg.solve)
    previous = probe(matrix, rhs)
    assert probe.cache.rebuilds == 1
    for step in range(1, 5):
        current = matrix.copy()
        current[:-1, :-1] += step * 0.001 * np.eye(6)
        updated_rhs = rhs.copy()
        updated_rhs[0] += step * 0.003
        result = probe(current, updated_rhs)
        np.testing.assert_allclose(result, np.linalg.solve(current, updated_rhs), atol=1e-12, rtol=0)
        assert not np.array_equal(result, previous)
    assert probe.cache.rebuilds == 1
    assert probe.cache.corrections > 0
    assert probe.max_charge_difference < 1e-12


def test_cache_refactors_after_large_change_and_atom_count_change():
    cache = CachedQEq(max_corrections=1)
    for size, seed in ((6, 42), (6, 123), (3, 57)):
        matrix, rhs = system(size, seed)
        np.testing.assert_allclose(cache.solve(matrix, rhs), np.linalg.solve(matrix, rhs), atol=1e-12, rtol=0)
    assert cache.rebuilds == 3


def test_cache_does_not_accept_nonfinite_matrix():
    matrix, rhs = system(3, 42)
    matrix[0, 0] = np.nan
    with pytest.raises(ValueError):
        CachedQEq().solve(matrix, rhs)
