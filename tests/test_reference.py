import json

from ase import Atoms
import numpy as np
import pytest

from cases import CASES
from validate import validation_cases, charged_validation_cases, backend_differences
from water_cluster import comparison
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps

pytestmark = pytest.mark.reference


@pytest.mark.parametrize("name", charged_validation_cases())
def test_supplied_charged_reference(name, tmp_path):
    filename, symbols, x, cell, pbc, charge = charged_validation_cases()[name]
    ff = ForceField.bundled(filename)
    # Exercise unwrapped branches and rotated cells as well as isolated ions.
    if cell is not None:
        x = x + np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]]) @ cell
    actual = Calculator(ff).evaluate(symbols, x, cell=cell, pbc=pbc, total_charge=charge)
    ref = evaluate_lammps(
        ff, symbols, x, cell=cell, pbc=pbc, supplied_charges=actual.charges, directory=tmp_path / name
    )
    report = comparison(actual, ref, len(x))
    assert report["passed"], report
    script = (ref.directory / "in.lammps").read_text()
    assert "checkqeq no" in script
    assert "fix charges" not in script
    metadata = json.loads((ref.directory / "metadata.json").read_text())
    assert metadata["charge_mode"] == "supplied"
    assert metadata["qeq_tolerance"] is None
    assert metadata["total_charge"] == pytest.approx(charge * np.prod(ref.cell_repetitions), abs=1e-12)
    np.testing.assert_allclose(ref.charges, actual.charges, atol=1e-12, rtol=0)


@pytest.mark.parametrize("charges", [[0], [[0, 0]], [float("nan"), 0], [0, float("inf")]])
def test_invalid_supplied_charges(charges, tmp_path):
    with pytest.raises(ValueError, match="supplied_charges"):
        evaluate_lammps(
            ForceField.bundled("ffield.reax.ZnOH.2010"),
            *CASES["zno"],
            supplied_charges=charges,
            directory=tmp_path,
        )
    assert not list(tmp_path.iterdir())


def test_missing_executable(tmp_path):
    with pytest.raises(FileNotFoundError):
        evaluate_lammps(
            ForceField.bundled("ffield.reax.ZnOH.2010"),
            *CASES["zno"],
            executable="nonexistent-xreac-lammps",
            directory=tmp_path,
        )


@pytest.mark.parametrize("name", ["monomer"])
def test_lammps_reequilibrated_energy_derivative(name, tmp_path):
    """Differentiate LAMMPS energies, with fresh QEq at each displacement."""
    from water_cluster import water_cases

    ff = ForceField.bundled("ffield.reax.HO.2015")
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


@pytest.mark.parametrize("name", validation_cases())
def test_system_reference(name, case_results, tmp_path):
    filename, symbols, x, cell, pbc = validation_cases()[name]
    results = case_results(name)
    reference = evaluate_lammps(results.force_field, symbols, x, cell=cell, pbc=pbc, directory=tmp_path / name)
    for actual in (results.native, results.ase):
        report = comparison(actual, reference, len(x))
        assert report["passed"], report
    if name.startswith("small_"):
        metadata = json.loads((reference.directory / "metadata.json").read_text())
        assert metadata["within_validated_cell_limits"] is True
        assert metadata["energy_divisor"] == np.prod(reference.cell_repetitions)
        assert metadata["input_atoms"] == len(x)
        assert np.loadtxt(reference.directory / "energy.txt")[0] / metadata["energy_divisor"] == reference.energy
    if name == "small_water_4A":
        assert results.native.components["hydrogen_bond"] == pytest.approx(-0.3700520183, abs=1e-9)
    if name == "carbon_dimer":
        assert results.native.components["lone_pair"] > 50  # C2 correction
    if name == "periodic_carbon_chain":
        assert abs(results.native.components["torsion"]) > 0.1
    if name == "small_zinc_chain":
        assert results.native.bond_counts.tolist() == [2]
        assert results.native.bond_orders[0, 0] > 1.0
        assert results.native.total_bond_orders[0] == results.native.bond_orders[0, 0]
        np.testing.assert_allclose(results.native.forces, 0, atol=1e-12)


@pytest.mark.parametrize("distance", [4.9999, 5.0, 5.0001, 9.9999, 10.0, 10.0001])
def test_cutoff_perturbations(distance, tmp_path):
    ff = ForceField.bundled("ffield.reax.ZnOH.2010")
    symbols, x = CASES["zno"]
    x = np.array(x, dtype=float)
    x[1, 0] = distance
    native = Calculator(ff).evaluate(symbols, x)
    atoms = Atoms(symbols, positions=x, calculator=ReaxFFCalculator(ff))
    atoms.get_forces()
    assert backend_differences(native, atoms.calc.evaluation, len(x))["passed"]
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path / "cutoff")
    assert comparison(native, reference, len(x))["passed"]
