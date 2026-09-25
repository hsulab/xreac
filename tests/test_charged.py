"""Net charge, independent constrained QEq, derivatives, and ASE integration."""

import numpy as np
import pytest
from ase import Atoms
from ase.neighborlist import neighbor_list
from ase.units import kcal, mol

from validate import charged_validation_cases, independent_charge_check, backend_differences
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator
from xreac.energy import EnergyModel

CASES = charged_validation_cases()


@pytest.fixture(scope="module")
def calc():
    return Calculator(ForceField.bundled("ffield.reax.HO.2015"))


@pytest.mark.parametrize("name", CASES)
@pytest.mark.parametrize("full_derivative", [False, True])
def test_charged_cases(calc, name, full_derivative):
    _, symbols, x, cell, pbc, charge = CASES[name]
    result = calc.evaluate(symbols, x, total_charge=charge, cell=cell, pbc=pbc, full_derivative=full_derivative)
    check = independent_charge_check(calc.force_field, symbols, x, cell, pbc, charge, result.charges)
    assert check["passed"], check
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc)
    edges = neighbor_list("ijS", atoms, calc.force_field.general[12])
    supplied = calc.evaluate(
        symbols, x, total_charge=charge, cell=cell, pbc=pbc, neighbors=edges, full_derivative=full_derivative
    )
    assert backend_differences(result, supplied, len(x))["passed"]
    atoms.calc = ReaxFFCalculator(calc.force_field, total_charge=charge, full_derivative=full_derivative)
    np.testing.assert_allclose(atoms.get_forces(), result.forces * kcal / mol, atol=1e-9)
    assert backend_differences(result, atoms.calc.evaluation, len(x))["passed"]

    model = EnergyModel(calc.force_field, symbols, edges, cell, pbc, total_charge=charge)
    fixed = None if full_derivative else result.charges
    direction = np.random.default_rng(937).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    derivative = (
        model.components(x + h * direction, fixed)[0].sum() - model.components(x - h * direction, fixed)[0].sum()
    ) / (2 * h)
    assert derivative == pytest.approx(-np.sum(result.forces * direction), abs=2e-5, rel=1e-6)


@pytest.mark.parametrize("charge", [-2, -0.25, 0, np.float64(0.5), 1])
def test_charge_values(calc, charge):
    _, symbols, x, cell, pbc, _ = CASES["hydroxide"]
    result = calc.evaluate(symbols, x, total_charge=charge)
    assert independent_charge_check(calc.force_field, symbols, x, cell, pbc, charge, result.charges)["passed"]
    atom = calc.evaluate(["O"], [[2, 3, 4]], total_charge=charge)
    np.testing.assert_allclose(atom.charges, [charge], atol=1e-12)
    np.testing.assert_allclose(atom.forces, 0, atol=1e-12)
    np.testing.assert_array_equal(atom.dipole, [0, 0, 0])
    if charge == 0:
        default = calc.evaluate(symbols, x)
        assert result.energy == default.energy
        np.testing.assert_array_equal(result.charges, default.charges)


@pytest.mark.parametrize(
    "charge",
    [True, np.bool_(False), "1", None, [1], np.array(1), np.array([1]), 1j, float("nan"), float("inf"), -float("inf")],
)
def test_invalid_charge(calc, charge):
    for method in (calc.evaluate, calc.relax):
        with pytest.raises(ValueError, match="total_charge"):
            method(["O"], [[0, 0, 0]], total_charge=charge)
    with pytest.raises(ValueError, match="total_charge"):
        ReaxFFCalculator(calc.force_field, total_charge=charge)


def test_ase_charge_initialization_and_cache(calc):
    _, symbols, x, _, _, charge = CASES["hydroxide"]
    atoms = Atoms(symbols, positions=x, calculator=ReaxFFCalculator(calc.force_field, total_charge=charge))
    original = atoms.get_potential_energy()
    atoms.set_initial_charges([-0.8, -0.2])
    assert atoms.get_potential_energy() == original
    atoms.set_initial_charges([-0.6, -0.4])
    assert atoms.get_potential_energy() == original
    for initial in ([0, 0], [float("nan"), -1]):
        atoms.set_initial_charges(initial)
        with pytest.raises(ValueError, match="total_charge"):
            atoms.get_charges()
        assert atoms.calc.results == {} and atoms.calc.evaluation is None
    atoms.set_initial_charges(None)
    atoms.calc.set(total_charge=0.5)
    assert atoms.calc.results == {} and atoms.calc.evaluation is None
    assert atoms.get_charges().sum() == pytest.approx(0.5, abs=1e-12)
    assert atoms.get_potential_energy() != original


@pytest.mark.parametrize("backend", ["ase", "native"])
@pytest.mark.parametrize("name", ["hydroxide", "hydronium"])
def test_relaxation_keeps_charge(calc, backend, name, monkeypatch):
    _, symbols, x, _, _, charge = CASES[name]
    evaluate = calc.evaluate
    seen = []

    def checked(*args, **kwargs):
        assert kwargs["total_charge"] == charge
        result = evaluate(*args, **kwargs)
        assert result.charges.sum() == pytest.approx(charge, abs=1e-12)
        seen.append(result)
        return result

    monkeypatch.setattr(calc, "evaluate", checked)
    relaxed = calc.relax(symbols, x, total_charge=charge, backend=backend, force_tolerance=1e-4)
    assert relaxed.converged
    assert len(seen) > 1
    expected = evaluate(symbols, relaxed.positions, total_charge=charge)
    assert backend_differences(relaxed.evaluation, expected, len(x))["passed"]


@pytest.mark.parametrize("name", ["charged_water_4A", "charged_partial_pbc_water"])
def test_periodic_charge_and_dipole(calc, name):
    _, symbols, x, cell, pbc, charge = CASES[name]
    original = calc.evaluate(symbols, x, cell=cell, pbc=pbc, total_charge=charge)
    masses = np.array([calc.force_field.atoms[s]["mass"] for s in symbols])
    center = np.average(x, axis=0, weights=masses)
    np.testing.assert_allclose(original.dipole, np.sum((x - center) * original.charges[:, None], axis=0), atol=1e-12)
    shift = np.array([[1, 0, 0], [-1, 1, 0], [0, -1, 0]]) @ cell
    wrapped = calc.evaluate(symbols, x + shift, cell=cell, pbc=pbc, total_charge=charge)
    assert wrapped.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(wrapped.forces, original.forces, atol=1e-8)
    expected = np.sum(shift * original.charges[:, None], axis=0) - charge * np.average(shift, axis=0, weights=masses)
    np.testing.assert_allclose(wrapped.dipole - original.dipole, expected, atol=1e-9)
    translated = calc.evaluate(symbols, x + [5, 3, -2], cell=cell, pbc=pbc, total_charge=charge)
    np.testing.assert_allclose(translated.dipole, original.dipole, atol=1e-9)
    doubled_cell = cell.copy()
    doubled_cell[0] *= 2
    doubled = calc.evaluate(
        symbols * 2, np.vstack((x, x + cell[0])), cell=doubled_cell, pbc=pbc, total_charge=2 * charge
    )
    assert doubled.energy == pytest.approx(2 * original.energy, abs=1e-9)
    for key in ("charges", "forces", "total_bond_orders", "lone_pairs", "bond_counts"):
        np.testing.assert_allclose(getattr(doubled, key), np.concatenate([getattr(original, key)] * 2), atol=1e-9)


@pytest.mark.parametrize("backend", ["ase", "native"])
def test_periodic_relaxation_keeps_charge(calc, backend):
    _, symbols, x, cell, pbc, charge = CASES["charged_water_4A"]
    relaxed = calc.relax(
        symbols,
        x,
        cell=cell,
        pbc=pbc,
        total_charge=charge,
        backend=backend,
        max_iterations=1,
        force_tolerance=1e-12,
    )
    assert relaxed.iterations == 1
    assert not relaxed.converged
    assert not np.array_equal(relaxed.positions, x)
    expected = calc.evaluate(symbols, relaxed.positions, cell=cell, pbc=pbc, total_charge=charge)
    assert expected.charges.sum() == pytest.approx(charge, abs=1e-12)
    assert backend_differences(relaxed.evaluation, expected, len(x))["passed"]
