import subprocess
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from ase.constraints import FixAtoms
from ase.optimize import BFGS, FIRE
from ase.units import kcal, mol

from water_cluster import water_cases
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps


@pytest.fixture
def atoms():
    symbols, positions = water_cases()["monomer"]
    result = Atoms(symbols, positions=positions)
    result.calc = ReaxFFCalculator(ForceField.bundled("qeq_ff.water"))
    return result


@pytest.mark.parametrize("full_derivative", [False, True])
def test_units_and_properties(atoms, full_derivative):
    atoms.calc.set(full_derivative=full_derivative)
    core = Calculator(ForceField.bundled("qeq_ff.water"))
    expected = core.evaluate(atoms.get_chemical_symbols(), atoms.positions, full_derivative=full_derivative)
    assert atoms.get_potential_energy() == pytest.approx(expected.energy*kcal/mol, abs=1e-12)
    np.testing.assert_allclose(atoms.get_forces(), expected.forces*kcal/mol, atol=1e-12)
    np.testing.assert_array_equal(atoms.get_charges(), expected.charges)
    np.testing.assert_array_equal(atoms.get_dipole_moment(), expected.dipole)
    assert atoms.calc.evaluation.full_derivative is full_derivative
    # Results use eV; retained Evaluation keeps core kcal/mol units.
    assert atoms.calc.evaluation.energy == expected.energy


def test_cache_and_parameter_changes(atoms, monkeypatch):
    evaluations = []
    original = atoms.calc.core.evaluate

    def counted(*args, **kwargs):
        result = original(*args, **kwargs)
        evaluations.append(result)
        return result

    monkeypatch.setattr(atoms.calc.core, "evaluate", counted)
    energy = atoms.get_potential_energy()
    atoms.get_forces()
    charges = atoms.get_charges()
    assert len(evaluations) == 1
    atoms.positions[1, 0] += .03
    assert atoms.get_potential_energy() != energy
    assert not np.allclose(atoms.get_charges(), charges)
    assert len(evaluations) == 2
    fixed = atoms.get_forces()
    atoms.calc.set(full_derivative=True)
    assert atoms.calc.evaluation is None
    assert np.max(abs(atoms.get_forces()-fixed)) > 1e-3
    assert len(evaluations) == 3
    atoms.calc.set(full_derivative=False)
    np.testing.assert_allclose(atoms.get_forces(), fixed, atol=1e-12)


def test_nonperiodic_box_and_unsupported_properties(atoms):
    energy = atoms.get_potential_energy()
    atoms.cell = [20., 20., 20.]
    assert atoms.get_potential_energy() == pytest.approx(energy, abs=1e-12)
    with pytest.raises(PropertyNotImplementedError):
        atoms.get_stress()
    atoms.pbc = True
    assert atoms.get_potential_energy() == pytest.approx(energy, abs=1e-12)
    atoms.cell = [0., 20., 20.]
    with pytest.raises(ValueError):
        atoms.get_potential_energy()
    assert atoms.calc.results == {}
    assert atoms.calc.evaluation is None


def test_charge_validation(atoms):
    energy = atoms.get_potential_energy()
    atoms.set_initial_charges([-.8, .4, .4])
    assert atoms.get_potential_energy() == pytest.approx(energy, abs=1e-12)
    atoms.set_initial_charges([0., .5, .5])
    with pytest.raises(ValueError, match="neutral"):
        atoms.get_charges()
    assert not atoms.calc.results


@pytest.mark.parametrize("kwargs", [
    {"full_derivative": "false"}, {"total_charge": 1}, {"misspelled_option": True},
])
def test_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        ReaxFFCalculator(ForceField.bundled("qeq_ff.water"), **kwargs)


@pytest.mark.parametrize("optimizer_class", [FIRE, BFGS])
def test_ase_optimizers_and_constraints(atoms, optimizer_class, monkeypatch):
    from autograd.tracer import Box
    from xreac import energy

    original = energy.np.linalg.solve

    def solve(matrix, rhs):
        assert not isinstance(matrix, Box), "ASE relaxation must not differentiate QEq"
        return original(matrix, rhs)

    monkeypatch.setattr(energy.np.linalg, "solve", solve)
    atoms.set_constraint(FixAtoms(indices=[0]))
    oxygen = atoms.positions[0].copy()
    with optimizer_class(atoms, logfile=None) as optimizer:
        assert optimizer.run(fmax=1e-5, steps=500)
    np.testing.assert_array_equal(atoms.positions[0], oxygen)
    assert np.max(np.linalg.norm(atoms.get_forces(), axis=1)) < 1e-5
    assert atoms.calc.evaluation.full_derivative is False


@pytest.mark.reference
def test_ase_lammps_reference(tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    symbols, positions = water_cases()["dimer"]
    atoms = Atoms(symbols, positions=positions, calculator=ReaxFFCalculator(ff))
    reference = evaluate_lammps(ff, symbols, positions, directory=tmp_path / "lammps")
    assert atoms.get_potential_energy() == pytest.approx(reference.energy*kcal/mol, abs=1e-8)
    np.testing.assert_allclose(atoms.get_forces(), reference.forces*kcal/mol, atol=1e-8, rtol=0)
    np.testing.assert_allclose(atoms.get_charges(), reference.charges, atol=1e-9, rtol=0)
    np.testing.assert_allclose(atoms.get_dipole_moment(), reference.dipole, atol=1e-9, rtol=0)


def test_native_backend_does_not_import_ase():
    code = """
import sys
from xreac import Calculator, ForceField
assert 'ase' not in sys.modules
result = Calculator(ForceField.zno()).relax(['Zn'], [[0, 0, 0]], backend='native')
assert result.converged
assert 'ase' not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True)
