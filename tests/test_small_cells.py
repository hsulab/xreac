import json

import numpy as np
import pytest

from audit_small_cells import differences, image_qeq, repeat_structure
from water_cluster import comparison, water_cases
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def test_explicit_small_cell_qeq_matches_supported_supercell():
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()["monomer"]
    cell = np.diag([4.]*3)
    primitive = image_qeq(ff, symbols, x, cell)
    ss, sx, sc = repeat_structure(symbols, x, cell, [3]*3)
    expanded = Calculator(ff).evaluate(ss, sx, cell=sc)
    np.testing.assert_allclose(expanded.charges, np.tile(primitive["charges"], 27), atol=1e-11, rtol=0)
    assert primitive["self_images_within_nonbonded_cutoff"] == [80]*3
    assert primitive["maximum_integer_shift_within_nonbonded_cutoff"] >= 2
    assert np.all(np.asarray(primitive["qeq_self_image_diagonal"]) > 0)
    assert abs(sum(primitive["charges"])) < 1e-12


def test_small_cell_guard_remains_enabled(tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()["monomer"]
    with pytest.raises(ValueError, match="Periodic cell heights"):
        Calculator(ff).evaluate(symbols, x, cell=[4.]*3)
    with pytest.raises(ValueError, match="Periodic cell heights"):
        evaluate_lammps(ff, symbols, x, cell=[4.]*3, directory=tmp_path/"reference")
    assert not (tmp_path/"reference").exists()
    with pytest.raises(ValueError, match="allow_small_cell"):
        evaluate_lammps(ff, symbols, x, cell=[4.]*3, allow_small_cell="yes")


@pytest.mark.reference
def test_lammps_small_cell_hydrogen_bond_exclusion(tmp_path):
    """Pinned LAMMPS excludes HB acceptors sharing the donor's original atom ID."""
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = water_cases()["monomer"]
    cell = np.diag([4.]*3)
    primitive = evaluate_lammps(ff, symbols, x, cell=cell, allow_small_cell=True,
                               directory=tmp_path/"primitive")
    ss, sx, sc = repeat_structure(symbols, x, cell, [3]*3)
    expanded = evaluate_lammps(ff, ss, sx, cell=sc, directory=tmp_path/"supercell")
    actual = Calculator(ff).evaluate(ss, sx, cell=sc)
    report = comparison(actual, expanded, len(sx))
    assert report["passed"], report
    diff = differences(primitive, expanded, 27, len(x))
    assert not diff["consistent"]
    assert diff["max_charges_difference"] < 1e-10
    assert diff["max_forces_difference"] > 1.
    assert diff["energy_difference_per_primitive_cell"] == pytest.approx(.3700520183, abs=1e-8)
    assert primitive.components["hydrogen_bond"] == 0
    assert expanded.components["hydrogen_bond"] < -9
    for key, value in diff["component_differences_per_primitive_cell"].items():
        if key != "hydrogen_bond":
            assert abs(value) < 1e-9, key
    metadata = json.loads((primitive.directory/"metadata.json").read_text())
    assert metadata["allow_small_cell"] is True
    assert metadata["within_validated_cell_limits"] is False
    metadata = json.loads((expanded.directory/"metadata.json").read_text())
    assert metadata["within_validated_cell_limits"] is True
