"""Lossless numerical export must retain fields discarded by model preprocessing."""

import numpy as np
import pytest

from xreac import ForceField


@pytest.mark.parametrize("name", ["ffield.reax.HO.2015", "ffield.reax.CHO.2008", "ffield.AgZnO"])
def test_source_roundtrip_without_source_file(tmp_path, name):
    bundled = ForceField.bundled(name)
    source = tmp_path / "source"
    source.write_bytes(bundled.path.read_bytes())
    parsed = ForceField.from_file(source)
    source.unlink()
    target = tmp_path / "written"
    parsed.to_file(target)
    again = ForceField.from_file(target)
    assert again.elements == parsed.elements
    assert again.section_counts == parsed.section_counts
    # Includes unused atom/bond/torsion fields and unexpanded wildcard records.
    assert again._sections == parsed._sections
    np.testing.assert_array_equal(again.general, parsed.general)
    assert again.atoms == parsed.atoms
    assert again.pairs == parsed.pairs
    for name in ("angles", "torsions", "hydrogen_bonds"):
        a, b = getattr(parsed, name), getattr(again, name)
        assert a.keys() == b.keys()
        for key in a:
            np.testing.assert_array_equal(a[key], b[key])


def test_export_rejects_edits_to_derived_tables(tmp_path):
    parsed = ForceField.bundled("ffield.reax.HO.2015")
    parsed.atoms["O"]["eta"] += 1
    with pytest.raises(ValueError, match="unmodified"):
        parsed.to_file(tmp_path / "modified")
    assert not (tmp_path / "modified").exists()
