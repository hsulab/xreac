"""Optional CuO regression: reuse the one surface fixture and an external force field."""

import os

from ase import Atoms
import pytest

from validate import backend_differences, validation_cases
from water_cluster import comparison
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps


@pytest.mark.reference
def test_cuo_surface_matches_lammps(tmp_path):
    path = os.environ.get("XREAC_CUO_FORCE_FIELD")
    if not path:
        pytest.skip("Set XREAC_CUO_FORCE_FIELD to the external Cu/O/H/Cl supplement")
    ff = ForceField.from_file(path)
    _, symbols, x, cell, pbc = validation_cases(include_cuo=True)["surface_cuo_010"]
    assert len(symbols) == 96 and symbols.count("Cu") == symbols.count("O") == 48
    # This distinction is absent from the earlier bundled validation cases.
    assert ff.atoms["Cu"]["valency_val"] != ff.atoms["Cu"]["valency_boc"]
    native = Calculator(ff).evaluate(symbols, x, cell=cell, pbc=pbc)
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc, calculator=ReaxFFCalculator(ff))
    atoms.get_forces()
    reference = evaluate_lammps(ff, symbols, x, cell=cell, pbc=pbc, directory=tmp_path / "reference")
    assert reference.cell_repetitions == (1, 1, 1)
    assert backend_differences(native, atoms.calc.evaluation, len(x))["passed"]
    for result in (native, atoms.calc.evaluation):
        report = comparison(result, reference, len(x))
        assert report["passed"], report
