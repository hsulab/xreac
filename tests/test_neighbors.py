import gzip
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.neighborlist import neighbor_list
from autograd import grad

from validate import validation_cases, backend_differences
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator, PrimitiveNeighborList, _directed_neighbors
from xreac.neighbors import Neighbors, scatter_sum, replicated_neighbors
from xreac.energy import EnergyModel

CASES = validation_cases()


@pytest.fixture(scope="module")
def archived_results():
    # Preserve an independent regression target from before model unification.
    root = Path(__file__).resolve().parents[1] / "validation"
    results = {}
    for system in ("water", "zno", "cho"):
        with gzip.open(root / system / "baseline.json.gz", "rt") as stream:
            results.update(json.load(stream))
    return results


def build_neighbors(symbols, x, cell, pbc, cutoff):
    return neighbor_list("ijS", Atoms(symbols, positions=x, cell=cell, pbc=pbc), np.nextafter(cutoff, np.inf))


@pytest.mark.parametrize("name", CASES)
def test_neighbor_backends(name, archived_results, case_results):
    filename, symbols, x, cell, pbc = CASES[name]
    ff = ForceField.bundled(filename)
    native_list, _ = replicated_neighbors(x, ff.general[12], cell, pbc)
    ase_list = build_neighbors(symbols, x, cell, pbc, ff.general[12])
    native_edges = Neighbors(native_list, len(x), cell, pbc)
    ase_edges = Neighbors(ase_list, len(x), cell, pbc)
    for key in ("i", "j", "shifts"):
        np.testing.assert_array_equal(getattr(native_edges, key), getattr(ase_edges, key))
    tree = PrimitiveNeighborList(
        np.full(len(x), np.nextafter(ff.general[12], np.inf) / 2),
        skin=0,
        self_interaction=False,
        bothways=False,
    )
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
    tree.update(atoms.pbc, atoms.cell, atoms.positions)
    tree_edges = Neighbors(_directed_neighbors(tree), len(x), cell, pbc)
    for key in ("i", "j", "shifts"):
        np.testing.assert_array_equal(getattr(tree_edges, key), getattr(ase_edges, key))
    expected = case_results(name).native
    actual = case_results(name).ase
    assert actual.neighbor_backend == "ase"
    assert actual.cell_repetitions == (1, 1, 1)
    report = backend_differences(actual, expected, len(x))
    assert report["passed"], report
    archived = SimpleNamespace(
        **{
            key: np.asarray(value) if isinstance(value, list) else value
            for key, value in archived_results[name].items()
        }
    )
    report = backend_differences(expected, archived, len(x))
    assert report["passed"], report
    supplied = Calculator(ff).evaluate(
        symbols, x, cell=cell, pbc=pbc, neighbors=build_neighbors(symbols, x, cell, pbc, ff.general[12])
    )
    assert supplied.neighbor_backend == "provided"
    assert backend_differences(supplied, expected, len(x))["passed"]
    if name == "periodic_carbon_chain":
        assert abs(actual.components["torsion"]) > 0.1
    if name == "small_water_4A":
        assert actual.components["hydrogen_bond"] < -0.3
    if name == "small_zinc_chain":
        assert actual.bond_counts.tolist() == [2]
        assert actual.bond_orders[0, 0] > 1


@pytest.mark.parametrize("name", ["small_water_4A", "small_partial_pbc_water", "periodic_carbon_chain"])
@pytest.mark.parametrize("full_derivative", [False, True])
def test_neighbor_derivatives(name, full_derivative, monkeypatch):
    from autograd.tracer import Box
    from xreac import energy

    filename, symbols, x, cell, pbc = CASES[name]
    ff = ForceField.bundled(filename)
    solve, solves = energy.np.linalg.solve, []

    def checked(matrix, rhs):
        assert matrix.shape == (len(x) + 1, len(x) + 1)
        solves.append(isinstance(matrix, Box))
        return solve(matrix, rhs)

    monkeypatch.setattr(energy.np.linalg, "solve", checked)
    calc = Calculator(ff)
    actual = calc.evaluate(
        symbols,
        x,
        cell=cell,
        pbc=pbc,
        full_derivative=full_derivative,
        neighbors=build_neighbors(symbols, x, cell, pbc, ff.general[12]),
    )
    assert any(solves) is full_derivative
    expected = calc.evaluate(symbols, x, cell=cell, pbc=pbc, full_derivative=full_derivative)
    assert backend_differences(actual, expected, len(x))["passed"]
    direction = np.random.default_rng(192).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    fixed = None if full_derivative else actual.charges

    def energy(y):
        neighbors = build_neighbors(symbols, y, cell, pbc, ff.general[12])
        model = EnergyModel(ff, symbols, neighbors, cell, pbc)
        return model.components(y, fixed)[0].sum()

    derivative = (energy(x + h * direction) - energy(x - h * direction)) / (2 * h)
    assert derivative == pytest.approx(-np.sum(actual.forces * direction), abs=2e-5, rel=1e-6)


def test_neighbor_shift_completeness():
    # Compare ASE edges to an independent explicit image search. Includes
    # repeated neighbors, nonzero-shift self edges, and a tilted small cell.
    cell = np.array([[3.1, 0.0, 0.0], [0.8, 3.6, 0.0], [0.2, 0.3, 4.1]])
    x = np.array([[0.0, 0.0, 0.0], [0.9, 0.2, 0.3]])
    neighbors = Neighbors(build_neighbors("OO", x, cell, [True, True, False], 7.5), 2, cell, [True, True, False])
    actual = {(int(i), int(j), *s) for i, j, s in zip(neighbors.i, neighbors.j, neighbors.shifts)}
    expected = set()
    for sx in range(-4, 5):
        for sy in range(-4, 5):
            s = np.array([sx, sy, 0])
            for i in range(2):
                for j in range(2):
                    if i == j and not s.any():
                        continue
                    if np.linalg.norm(x[j] - x[i] + s @ cell) < 7.5:
                        expected.add((i, j, *s))
    assert actual == expected
    assert max(abs(sx) for _, _, sx, _, _ in actual) >= 2
    np.testing.assert_array_equal(neighbors.reverse[neighbors.reverse], np.arange(len(neighbors.i)))


def test_scatter_sum_derivative():
    indices = np.array([0, 1, 0, 3])
    x = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(scatter_sum(x, indices, 5), [4, 2, 0, 4, 0])
    derivative = grad(lambda y: np.sum(scatter_sum(y, indices, 5) ** 2))(x)
    np.testing.assert_array_equal(derivative, [8, 4, 8, 8])


@pytest.mark.parametrize("shift", [1, 2**62])
@pytest.mark.parametrize("index_dtype", [np.int64, np.uint64])
def test_neighbor_order_with_large_image_shifts(shift, index_dtype):
    # Exercise ordinary packed keys and the overflow-safe fallback using only
    # integer topology, without introducing another simulation structure.
    i = np.array([1, 0, 0, 1, 0, 1], dtype=index_dtype)
    j = 1 - i
    shifts = np.array([[shift, 0, 0], [0, 0, 0], [-shift, 0, 0], [0, 0, 0], [-shift - 1, 0, 0], [shift + 1, 0, 0]])
    edges = Neighbors((i, j, shifts), 2, np.eye(3), True)
    keys = list(zip(edges.i, edges.j, *edges.shifts.T))
    assert keys == sorted(zip(i, j, *shifts.T))
    np.testing.assert_array_equal(edges.i[edges.reverse], edges.j)
    np.testing.assert_array_equal(edges.shifts[edges.reverse], -edges.shifts)
    np.testing.assert_array_equal(edges.reverse[edges.reverse], np.arange(len(i)))
    with pytest.raises(ValueError, match="both directions"):
        Neighbors((i[:-1], j[:-1], shifts[:-1]), 2, np.eye(3), True)
    with pytest.raises(ValueError, match="Duplicate"):
        Neighbors((np.r_[i, i[:1]], np.r_[j, j[:1]], np.concatenate((shifts, shifts[:1]))), 2, np.eye(3), True)


def test_ase_default_no_replication_and_backend_switch(monkeypatch):
    filename, symbols, x, cell, pbc = CASES["small_water_4A"]
    ff = ForceField.bundled(filename)
    from xreac.geometry import Boundary

    original = Boundary.supercell
    with monkeypatch.context() as patch:

        def forbidden(*args, **kwargs):
            raise AssertionError("ASE neighbor lists must not replicate atoms")

        patch.setattr(Boundary, "supercell", forbidden)
        atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc, calculator=ReaxFFCalculator(ff, max_expanded_atoms=1))
        forces = atoms.get_forces()
        assert atoms.calc.evaluation.neighbor_backend == "ase"
    assert Boundary.supercell is original
    atoms.calc.set(neighbor_backend="replicated")
    assert atoms.calc.evaluation is None and not atoms.calc.results
    with pytest.raises(ValueError, match="max_expanded_atoms"):
        atoms.get_forces()
    atoms.calc.set(max_expanded_atoms=81)
    np.testing.assert_allclose(atoms.get_forces(), forces, atol=1e-11, rtol=0)
    assert atoms.calc.evaluation.neighbor_backend == "replicated"
    assert atoms.calc.evaluation.cell_repetitions == (3, 3, 3)


@pytest.mark.parametrize("bad", ["unknown", None, True])
def test_invalid_neighbor_backend(bad):
    ff = ForceField.bundled("ffield.reax.ZnOH.2010")
    with pytest.raises(ValueError, match="neighbor_backend"):
        ReaxFFCalculator(ff, neighbor_backend=bad)


@pytest.mark.parametrize("name", ["small_partial_pbc_water", "periodic_carbon_chain"])
def test_neighbor_symmetries(name):
    filename, symbols, x, cell, pbc = CASES[name]
    calc = Calculator(ForceField.bundled(filename))

    def evaluate(s, y, lattice):
        return calc.evaluate(
            s, y, cell=lattice, pbc=pbc, neighbors=build_neighbors(s, y, lattice, pbc, calc.force_field.general[12])
        )

    original = evaluate(symbols, x, cell)
    rng = np.random.default_rng(257)
    shifts = rng.integers(-3, 4, x.shape) * pbc
    wrapped = evaluate(symbols, x + shifts @ cell, cell)
    np.testing.assert_allclose(wrapped.forces, original.forces, atol=2e-8, rtol=0)
    np.testing.assert_allclose(wrapped.charges, original.charges, atol=1e-11, rtol=0)
    np.testing.assert_allclose(wrapped.bond_orders, original.bond_orders, atol=1e-10, rtol=0)
    assert wrapped.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(
        wrapped.dipole - original.dipole, np.sum(original.charges[:, None] * (shifts @ cell), axis=0), atol=1e-9
    )
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotated = evaluate(symbols, x @ rotation + [8.0, -5.0, 20.0], cell @ rotation)
    assert rotated.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(rotated.forces, original.forces @ rotation, atol=2e-8, rtol=0)
    order = rng.permutation(len(x))
    permuted = evaluate([symbols[i] for i in order], x[order], cell)
    assert permuted.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(permuted.forces, original.forces[order], atol=2e-8, rtol=0)


def test_neighbor_rebuild_crossing_cutoff():
    atoms = Atoms(
        "ZnO",
        positions=[[0, 0, 0], [10.1, 0, 0]],
        calculator=ReaxFFCalculator(ForceField.bundled("ffield.reax.ZnOH.2010")),
    )
    atoms.get_forces()
    distant_energy = atoms.get_potential_energy()
    atoms.positions[1] = [1.9, 0.1, 0.2]
    assert atoms.get_potential_energy() != distant_energy
    reference = Calculator(atoms.calc.core.force_field).evaluate(atoms.get_chemical_symbols(), atoms.positions)
    assert backend_differences(atoms.calc.evaluation, reference, 2)["passed"]
    atoms.positions[1] = [10.1, 0, 0]
    assert atoms.get_potential_energy() == pytest.approx(distant_energy, abs=1e-12)


def test_ase_builds_arrays_before_core_evaluation(monkeypatch):
    from autograd.tracer import Box
    from xreac import ase as adapter

    ff = ForceField.bundled("ffield.reax.HO.2015")
    filename, symbols, x, cell, pbc = CASES["small_water_4A"]
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc, calculator=ReaxFFCalculator(ff))
    original_builder = adapter.PrimitiveNeighborList.build
    original_evaluate = atoms.calc.core.evaluate
    calls = []

    def build(nl, pbc, cell, positions):
        assert not isinstance(positions, Box)
        calls.append("build")
        return original_builder(nl, pbc, cell, positions)

    def evaluate(*args, **kwargs):
        calls.append("evaluate")
        assert calls == ["build", "evaluate"]
        assert len(kwargs["neighbors"]) == 3
        assert all(v.dtype.kind in "iu" for v in kwargs["neighbors"])
        # No builder call is permitted while evaluating energies or forces.
        with monkeypatch.context() as patch:

            def forbidden(*args, **kwargs):
                raise AssertionError("Neighbor builder called inside core evaluate")

            patch.setattr(adapter.PrimitiveNeighborList, "build", forbidden)
            patch.setattr("ase.neighborlist.primitive_neighbor_list", forbidden)
            return original_evaluate(*args, **kwargs)

    monkeypatch.setattr(adapter.PrimitiveNeighborList, "build", build)
    monkeypatch.setattr(atoms.calc.core, "evaluate", evaluate)
    atoms.get_forces()
    atoms.get_potential_energy()
    assert calls == ["build", "evaluate"]


def test_supplied_arrays_work_without_importing_ase():
    code = """
import sys
import numpy as np
from xreac import Calculator, ForceField
assert 'ase' not in sys.modules
calc = Calculator(ForceField.bundled("ffield.reax.ZnOH.2010"))
shifts = np.array([[s, 0, 0] for s in (-4, -3, -2, -1, 1, 2, 3, 4)])
neighbors = (np.zeros(8, dtype=int), np.zeros(8, dtype=int), shifts)
actual = calc.evaluate(['Zn'], [[0, 0, 0]], cell=[2.5, 12, 12], pbc=[True, False, False], neighbors=neighbors)
expected = calc.evaluate(['Zn'], [[0, 0, 0]], cell=[2.5, 12, 12], pbc=[True, False, False])
assert abs(actual.energy-expected.energy) < 1e-10
assert actual.bond_counts.tolist() == [2]
assert 'ase' not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True)


@pytest.mark.parametrize(
    "neighbors,message",
    [
        (([0], [1]), "tuple"),
        (([0], [1, 0], [[0, 0, 0]]), "shape"),
        (([0], [1], [[0, 0]]), "shape"),
        (([0.0, 1.0], [1, 0], [[0, 0, 0], [0, 0, 0]]), "integers"),
        (([0, 1], [1, 0], [[0.0, 0, 0], [0, 0, 0]]), "integers"),
        (([-1, 1], [1, 0], [[0, 0, 0], [0, 0, 0]]), "out of range"),
        (([0, 1], [2, 0], [[0, 0, 0], [0, 0, 0]]), "out of range"),
        (([0, 1], [1, 0], [[1, 0, 0], [-1, 0, 0]]), "nonperiodic"),
        (([0], [0], [[0, 0, 0]]), "self neighbors"),
        (([0, 0, 1], [1, 1, 0], [[0, 0, 0]] * 3), "Duplicate"),
        (([0], [1], [[0, 0, 0]]), "both directions"),
    ],
)
def test_supplied_neighbor_validation(neighbors, message):
    with pytest.raises(ValueError, match=message):
        Calculator(ForceField.bundled("ffield.reax.ZnOH.2010")).evaluate(
            ["Zn", "O"], [[0, 0, 0], [1.9, 0, 0]], neighbors=neighbors
        )


def test_reuse_supplied_skin_list_for_displacements():
    filename, symbols, x, cell, pbc = CASES["small_water_4A"]
    ff = ForceField.bundled(filename)
    supplied = build_neighbors(symbols, x, cell, pbc, ff.general[12] + 1.0)
    originals = tuple(array.copy() for array in supplied)
    # The same integer arrays can be reused over small displacements, provided
    # the caller keeps the list complete. Distances must be recomputed from x.
    calc = Calculator(ff)
    for shift in (0.0, 0.03, -0.02):
        y = x.copy()
        y[1, 0] += shift
        actual = calc.evaluate(symbols, y, cell=cell, pbc=pbc, neighbors=supplied)
        expected = calc.evaluate(symbols, y, cell=cell, pbc=pbc)
        assert backend_differences(actual, expected, len(x))["passed"]
    for original, array in zip(originals, supplied):
        np.testing.assert_array_equal(original, array)
