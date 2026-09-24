"""Ensure benchmark timings include actual ASE neighbor construction and evaluation."""

import numpy as np

from benchmark_lammps import ase_timing
from validate import PILOT_CASES, validation_cases
from xreac import Calculator, ForceField


def test_ase_timing_bypasses_result_cache(monkeypatch):
    from xreac import ase as adapter
    from xreac import calculator

    calls = []
    original = adapter.PrimitiveNeighborList.build

    def build(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    def forbidden(*args, **kwargs):
        raise AssertionError("ASE benchmarks must not use native neighbor construction")

    monkeypatch.setattr(adapter.PrimitiveNeighborList, "build", build)
    monkeypatch.setattr(calculator, "replicated_neighbors", forbidden)
    filename, symbols, x, cell, pbc = validation_cases()["cluster_water_monomer"]
    calc = Calculator(ForceField.bundled(filename))
    result, timing = ase_timing(calc, symbols, x, cell, pbc, repeats=2, target_seconds=1e-9)
    assert len(calls) == 1 + 2 * timing["calls_per_batch"]
    assert result.neighbor_backend == timing["neighbor_backend"] == "ase"
    assert timing["neighbor_build_included"] and not timing["result_cache_used"]
    assert np.isfinite(result.forces).all()


def test_pilot_cases_are_three_shared_structures():
    cases = validation_cases(include_bulk=True, include_cuo=True)
    assert len(PILOT_CASES) == 3
    assert [len(cases[name][1]) for name in PILOT_CASES] == [192, 128, 96]
    assert not set(PILOT_CASES).intersection(validation_cases())
    for name in PILOT_CASES:
        _, symbols, x, cell, pbc = cases[name]
        heights = 1 / np.linalg.norm(np.linalg.inv(cell), axis=0)
        assert np.all(heights[np.asarray(pbc)] > 10)
        assert np.asarray(x).shape == (len(symbols), 3)
    assert cases["periodic_bulk_zno_128"][1].count("Zn") == 64
    assert cases["periodic_bulk_zno_128"][1].count("O") == 64
