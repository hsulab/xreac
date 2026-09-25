"""Published Hur parameter coverage and shared LAMMPS comparison geometries."""

import numpy as np
import pytest

from validate import CHOCL_FORCE_FIELD, chocl_validation_cases, independent_charge_check
from water_cluster import comparison
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def test_hur_parameters():
    ff = ForceField.bundled(CHOCL_FORCE_FIELD)
    assert ff.checksum == "d468f31f71c1ca429aef9e576bf2c5fda3d8eb8ad0c7c4ae433197692ee5825b"
    assert ff.elements == ("C", "H", "O", "Cl", "X")
    assert "10.1039/D1RA04397H" in ff.citation
    assert ff.atoms["Cl"]["chi"] == 6.7204
    assert ff.atoms["Cl"]["eta"] == 2 * 6.1703
    assert ff.pairs["C", "Cl"]["De_s"] == 100.2833
    assert ff.pairs["H", "Cl"]["De_s"] == 103.25
    assert len(ff.angles["H", "C", "Cl"]) == 2
    assert len(ff.hydrogen_bonds) == 3
    Calculator(ff)


@pytest.mark.reference
@pytest.mark.parametrize("name", chocl_validation_cases())
def test_hur_reference(name, tmp_path):
    filename, symbols, x, cell, pbc, charge = chocl_validation_cases()[name]
    ff = ForceField.bundled(filename)
    result = Calculator(ff).evaluate(symbols, x, total_charge=charge)
    assert result.charges.sum() == pytest.approx(charge, abs=1e-12)
    check = independent_charge_check(ff, symbols, x, cell, pbc, charge, result.charges)
    assert check["passed"], check
    reference = evaluate_lammps(
        ff, symbols, x, supplied_charges=result.charges if charge else None, directory=tmp_path / name
    )
    report = comparison(result, reference, len(x))
    assert report["passed"], report
    if name == "chocl_sn2_symmetric":
        assert result.charges[-1] == pytest.approx(result.charges[-2], abs=1e-12)
    full = Calculator(ff).evaluate(symbols, x, total_charge=charge, full_derivative=True)
    direction = np.random.default_rng(79).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    calc = Calculator(ff)
    fd = (
        calc.evaluate(symbols, x + h * direction, total_charge=charge).energy
        - calc.evaluate(symbols, x - h * direction, total_charge=charge).energy
    ) / (2 * h)
    assert fd == pytest.approx(-np.sum(full.forces * direction), abs=2e-5, rel=1e-6)
