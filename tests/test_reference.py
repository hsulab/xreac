import numpy as np
import pytest

from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps
from cases import CASES

pytestmark = pytest.mark.reference


@pytest.mark.parametrize("name", list(CASES))
def test_reference(name, tmp_path):
    ff = ForceField.zno()
    symbols, x = CASES[name]
    actual = Calculator(ff).evaluate(symbols, x)
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path / name)
    n = len(symbols)
    assert actual.energy == pytest.approx(reference.energy, abs=1e-5 * n)
    for term in actual.components:
        assert actual.components[term] == pytest.approx(reference.components[term], abs=1e-5 * n), term
    np.testing.assert_allclose(actual.charges, reference.charges, atol=1e-6, rtol=0)
    np.testing.assert_allclose(actual.forces, reference.forces, atol=1e-4, rtol=0)


@pytest.mark.parametrize("distance", [1.5, 2.5, 3.0, 3.5, 4.0, 4.9999, 5.0001, 9.999, 10.0, 10.001])
def test_stretch_and_cutoffs(distance, tmp_path):
    ff = ForceField.zno()
    symbols, x = ["Zn", "O"], [[0, 0, 0], [distance, 0, 0]]
    actual = Calculator(ff).evaluate(symbols, x)
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path / "scan")
    assert actual.energy == pytest.approx(reference.energy, abs=2e-5)
    np.testing.assert_allclose(actual.forces, reference.forces, atol=1e-4, rtol=0)


@pytest.mark.parametrize("name", ["zno", "cluster20"])
def test_reference_relaxed_geometry(name, tmp_path):
    ff = ForceField.zno()
    calc = Calculator(ff)
    symbols, x = CASES[name]
    relaxed = calc.relax(symbols, x)
    assert relaxed.converged
    reference = evaluate_lammps(ff, symbols, relaxed.positions, directory=tmp_path / "relaxed")
    assert relaxed.evaluation.energy == pytest.approx(reference.energy, abs=1e-5 * len(symbols))
    assert relaxed.evaluation.full_derivative is False
    np.testing.assert_allclose(relaxed.evaluation.forces, reference.forces, atol=1e-4, rtol=0)
    assert np.max(abs(reference.forces)) <= 1e-4 + 1e-7


def test_missing_executable(tmp_path):
    with pytest.raises(FileNotFoundError):
        evaluate_lammps(ForceField.zno(), *CASES["zno"], executable="nonexistent-xreac-lammps", directory=tmp_path)


@pytest.mark.parametrize("name", ["monomer", "dimer"])
def test_lammps_reequilibrated_energy_derivative(name, tmp_path):
    """Differentiate LAMMPS energies, with fresh QEq at each displacement."""
    from water_cluster import water_cases

    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()[name]
    actual = Calculator(ff).evaluate(symbols, x)
    full = Calculator(ff).evaluate(symbols, x, full_derivative=True)
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path / "center")
    step = 1e-5
    displacement = np.zeros_like(x)
    displacement[0, 1] = step
    plus = evaluate_lammps(ff, symbols, x + displacement, directory=tmp_path / "plus")
    minus = evaluate_lammps(ff, symbols, x - displacement, directory=tmp_path / "minus")
    numerical_force = -(plus.energy - minus.energy) / (2 * step)
    assert np.max(abs(plus.charges - minus.charges)) > 1e-6
    assert numerical_force == pytest.approx(full.forces[0, 1], abs=1e-6, rel=0)
    assert reference.forces[0, 1] == pytest.approx(actual.forces[0, 1], abs=1e-8, rel=0)
    assert abs(numerical_force - reference.forces[0, 1]) > 0.1
