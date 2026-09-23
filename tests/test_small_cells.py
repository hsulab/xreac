import json

import numpy as np
import pytest
from ase import Atoms
from ase.units import kcal, mol

from audit_small_cells import differences, image_qeq, repeat_structure
from small_cells import small_cell_cases
from water_cluster import comparison, water_cases
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator
from xreac.energy import EnergyModel
from xreac.neighbors import replicated_neighbors
from xreac.reference import evaluate_lammps

CASES = small_cell_cases()


def test_explicit_small_cell_qeq_matches_supported_supercell():
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()["monomer"]
    cell = np.diag([4.0] * 3)
    primitive = image_qeq(ff, symbols, x, cell)
    ss, sx, sc = repeat_structure(symbols, x, cell, [3] * 3)
    expanded = Calculator(ff).evaluate(ss, sx, cell=sc)
    np.testing.assert_allclose(expanded.charges, np.tile(primitive["charges"], 27), atol=1e-11, rtol=0)
    assert primitive["self_images_within_nonbonded_cutoff"] == [80] * 3
    assert primitive["maximum_integer_shift_within_nonbonded_cutoff"] >= 2
    assert np.all(np.asarray(primitive["qeq_self_image_diagonal"]) > 0)
    assert abs(sum(primitive["charges"])) < 1e-12


def test_small_cell_memory_limit(tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()["monomer"]
    with pytest.raises(ValueError, match="max_expanded_atoms"):
        Calculator(ff, max_expanded_atoms=80).evaluate(symbols, x, cell=[4.0] * 3)
    with pytest.raises(ValueError, match="max_expanded_atoms"):
        evaluate_lammps(ff, symbols, x, cell=[4.0] * 3, directory=tmp_path / "reference", max_expanded_atoms=80)
    assert not (tmp_path / "reference").exists()
    with pytest.raises(ValueError, match="allow_small_cell"):
        evaluate_lammps(ff, symbols, x, cell=[4.0] * 3, allow_small_cell="yes")


@pytest.mark.reference
def test_lammps_small_cell_hydrogen_bond_exclusion(tmp_path):
    """Pinned LAMMPS excludes HB acceptors sharing the donor's original atom ID."""
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()["monomer"]
    cell = np.diag([4.0] * 3)
    primitive = evaluate_lammps(ff, symbols, x, cell=cell, allow_small_cell=True, directory=tmp_path / "primitive")
    ss, sx, sc = repeat_structure(symbols, x, cell, [3] * 3)
    expanded = evaluate_lammps(ff, ss, sx, cell=sc, directory=tmp_path / "supercell")
    actual = Calculator(ff).evaluate(ss, sx, cell=sc)
    report = comparison(actual, expanded, len(sx))
    assert report["passed"], report
    diff = differences(primitive, expanded, 27, len(x))
    assert not diff["consistent"]
    assert diff["max_charges_difference"] < 1e-10
    assert diff["max_forces_difference"] > 1.0
    assert diff["energy_difference_per_primitive_cell"] == pytest.approx(0.3700520183, abs=1e-8)
    assert primitive.components["hydrogen_bond"] == 0
    assert expanded.components["hydrogen_bond"] < -9
    for key, value in diff["component_differences_per_primitive_cell"].items():
        if key != "hydrogen_bond":
            assert abs(value) < 1e-9, key
    metadata = json.loads((primitive.directory / "metadata.json").read_text())
    assert metadata["allow_small_cell"] is True
    assert metadata["within_validated_cell_limits"] is False
    metadata = json.loads((expanded.directory / "metadata.json").read_text())
    assert metadata["within_validated_cell_limits"] is True


@pytest.mark.reference
@pytest.mark.parametrize("name", CASES)
def test_supported_small_cells(name, tmp_path):
    filename, symbols, x, cell, pbc = CASES[name]
    ff = ForceField.bundled(filename)
    result = Calculator(ff).evaluate(symbols, x, cell=cell, pbc=pbc)
    ref = evaluate_lammps(ff, symbols, x, cell=cell, pbc=pbc, directory=tmp_path / name)
    report = comparison(result, ref, len(x))
    assert report["passed"], report
    assert ref.cell_repetitions == result.cell_repetitions
    metadata = json.loads((ref.directory / "metadata.json").read_text())
    assert metadata["within_validated_cell_limits"] is True
    assert metadata["energy_divisor"] == np.prod(result.cell_repetitions)
    assert metadata["input_atoms"] == len(x)
    assert np.loadtxt(ref.directory / "energy.txt")[0] / metadata["energy_divisor"] == ref.energy
    if name == "water_4A":
        assert result.components["hydrogen_bond"] == pytest.approx(-0.3700520183, abs=1e-9)
    if name == "zinc_chain":
        assert result.bond_orders[0, 0] > 1.0
        assert result.bond_counts.tolist() == [2]
        assert result.total_bond_orders[0] == result.bond_orders[0, 0]
        np.testing.assert_allclose(result.forces, 0.0, atol=1e-12)


@pytest.mark.parametrize("name", ["water_4A", "partial_pbc_water"])
@pytest.mark.parametrize("full_derivative", [False, True])
def test_small_cell_derivatives_and_reduced_qeq(name, full_derivative, monkeypatch):
    from autograd.tracer import Box
    from xreac import energy

    filename, symbols, x, cell, pbc = CASES[name]
    ff = ForceField.bundled(filename)
    solve = energy.np.linalg.solve
    solves = []

    def counted(matrix, rhs):
        assert matrix.shape == (len(x) + 1, len(x) + 1), "QEq must solve only primitive charges"
        solves.append(isinstance(matrix, Box))
        return solve(matrix, rhs)

    monkeypatch.setattr(energy.np.linalg, "solve", counted)
    result = Calculator(ff).evaluate(symbols, x, cell=cell, pbc=pbc, full_derivative=full_derivative)
    assert any(solves) is full_derivative
    neighbors, _ = replicated_neighbors(x, ff.general[12], cell, pbc)
    model = EnergyModel(ff, symbols, neighbors, cell, pbc)
    direction = np.random.default_rng(130).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    charges = None if full_derivative else result.charges
    plus = model.components(x + h * direction, charges)[0].sum()
    minus = model.components(x - h * direction, charges)[0].sum()
    assert (plus - minus) / (2 * h) == pytest.approx(-np.sum(result.forces * direction), abs=2e-5, rel=1e-6)


@pytest.mark.parametrize("name", ["zno_4A", "triclinic_water", "partial_pbc_water"])
def test_small_cell_symmetries(name):
    filename, symbols, x, cell, pbc = CASES[name]
    calc = Calculator(ForceField.bundled(filename))
    original = calc.evaluate(symbols, x, cell=cell, pbc=pbc)
    rng = np.random.default_rng(7)
    shifts = rng.integers(-2, 3, x.shape) * pbc
    wrapped = calc.evaluate(symbols, x + shifts @ cell, cell=cell, pbc=pbc)
    assert wrapped.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(wrapped.forces, original.forces, atol=2e-8, rtol=0)
    np.testing.assert_allclose(wrapped.charges, original.charges, atol=1e-10, rtol=0)
    np.testing.assert_allclose(wrapped.bond_orders, original.bond_orders, atol=1e-10, rtol=0)
    np.testing.assert_allclose(
        wrapped.dipole - original.dipole, np.sum(original.charges[:, None] * (shifts @ cell), axis=0), atol=1e-9
    )
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    transformed = calc.evaluate(symbols, x @ rotation + [12.0, -8.0, 3.0], cell=cell @ rotation, pbc=pbc)
    assert transformed.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(transformed.forces, original.forces @ rotation, atol=2e-8, rtol=0)
    order = rng.permutation(len(x))
    permuted = calc.evaluate([symbols[i] for i in order], x[order], cell=cell, pbc=pbc)
    assert permuted.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(permuted.forces, original.forces[order], atol=2e-8, rtol=0)
    np.testing.assert_allclose(permuted.bond_orders, original.bond_orders[np.ix_(order, order)], atol=1e-10)


def test_small_cell_replication_extensivity():
    filename, symbols, x, cell, pbc = CASES["water_4A"]
    calc = Calculator(ForceField.bundled(filename))
    original = calc.evaluate(symbols, x, cell=cell)
    ss, sx, sc = repeat_structure(symbols, x, cell, [2, 1, 1])
    doubled = calc.evaluate(ss, sx, cell=sc)
    # The two calculations use different internal supercell sizes (81/108).
    assert np.prod(original.cell_repetitions) * len(x) != np.prod(doubled.cell_repetitions) * len(sx)
    assert doubled.energy == pytest.approx(2 * original.energy, abs=1e-9)
    np.testing.assert_allclose(doubled.forces, np.tile(original.forces, (2, 1)), atol=1e-9, rtol=0)
    np.testing.assert_allclose(doubled.charges, np.tile(original.charges, 2), atol=1e-11, rtol=0)


@pytest.mark.parametrize("backend", ["ase", "native"])
@pytest.mark.reference
def test_small_cell_relaxation(backend, tmp_path, monkeypatch):
    from autograd.tracer import Box
    from xreac import energy

    filename, symbols, x, cell, pbc = CASES["water_6A"]
    ff = ForceField.bundled(filename)
    solve = energy.np.linalg.solve

    def fixed_only(matrix, rhs):
        assert not isinstance(matrix, Box)
        assert matrix.shape == (4, 4)
        return solve(matrix, rhs)

    monkeypatch.setattr(energy.np.linalg, "solve", fixed_only)
    # Cubic symmetry fixes the molecular orientation, avoiding slow rotation
    # under weak image torques while testing bond/angle relaxation to tolerance.
    half = np.deg2rad(52.25)
    x = np.array(
        [[0.0, 0, 0], [0.97 * np.cos(half), 0.97 * np.sin(half), 0], [0.97 * np.cos(half), -0.97 * np.sin(half), 0]]
    )
    result = Calculator(ff).relax(symbols, (x + 5.6) % 6, cell=cell, backend=backend, force_tolerance=1e-5)
    assert result.converged
    assert result.evaluation.cell_repetitions == ((1, 1, 1) if backend == "ase" else (2, 2, 2))
    assert result.evaluation.neighbor_backend == ("ase" if backend == "ase" else "replicated")
    assert not result.evaluation.full_derivative
    ref = evaluate_lammps(ff, symbols, result.positions, cell=cell, directory=tmp_path / backend)
    report = comparison(result.evaluation, ref, len(x))
    assert report["passed"], report
    assert np.max(abs(ref.forces)) < 1.01e-5


@pytest.mark.parametrize("backend", ["ase", "native"])
def test_small_cell_relaxation_iteration_limit(backend):
    filename, symbols, x, cell, pbc = CASES["water_6A"]
    calc = Calculator(ForceField.bundled(filename))
    result = calc.relax(symbols, x, cell=cell, backend=backend, max_iterations=3)
    assert not result.converged and result.iterations == 3
    assert not result.evaluation.full_derivative
    final = calc.evaluate(symbols, result.positions, cell=cell)
    np.testing.assert_allclose(result.evaluation.forces, final.forces, atol=1e-12)


def test_ase_small_cell_limit_and_caching():
    filename, symbols, x, cell, pbc = CASES["water_4A"]
    ff = ForceField.bundled(filename)
    atoms = Atoms(
        symbols,
        positions=x,
        cell=cell,
        pbc=pbc,
        calculator=ReaxFFCalculator(ff, max_expanded_atoms=80, neighbor_backend="replicated"),
    )
    with pytest.raises(ValueError, match="max_expanded_atoms"):
        atoms.get_forces()
    atoms.calc.set(max_expanded_atoms=81)
    expected = Calculator(ff).evaluate(symbols, x, cell=cell)
    np.testing.assert_allclose(atoms.get_forces(), expected.forces * kcal / mol, atol=1e-12)
    assert atoms.calc.evaluation.cell_repetitions == (3, 3, 3)
    atoms.cell = [12.0] * 3
    atoms.get_forces()
    assert atoms.calc.evaluation.cell_repetitions == (1, 1, 1)


@pytest.mark.parametrize("limit", [0, -1, True, 3.5])
def test_invalid_expansion_limit(limit):
    ff = ForceField.zno()
    with pytest.raises(ValueError, match="max_expanded_atoms"):
        Calculator(ff, max_expanded_atoms=limit)
    with pytest.raises(ValueError, match="max_expanded_atoms"):
        ReaxFFCalculator(ff, max_expanded_atoms=limit)
