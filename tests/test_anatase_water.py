"""Check slab termination, stoichiometry, constraints and adsorption sites."""

import numpy as np
import pytest
from pathlib import Path
from ase.io import read
from ase.neighborlist import neighbor_list

from anatase_water import FORCE_FIELD, adsorption_status, geometry, verify
from xreac import ForceField

INPUTS = Path(__file__).resolve().parents[1] / "validation/anatase"


@pytest.mark.parametrize("waters", [1, 2, 4, 8])
def test_anatase_geometry(waters):
    atoms = geometry(waters)
    fixture = read(INPUTS / f"{waters}water_initial.xyz")
    np.testing.assert_allclose(atoms.positions, fixture.positions, rtol=0, atol=1e-8)
    np.testing.assert_array_equal(atoms.numbers, fixture.numbers)
    slab = atoms[:192]
    assert len(atoms) == 192 + 3 * waters
    assert sum(slab.numbers == 22) == 64
    assert sum(slab.numbers == 8) == 128
    assert np.array_equal(atoms.pbc, [True, True, False])
    assert np.array_equal(np.bincount(slab.get_tags()), [0, 48, 48, 48, 48])
    assert np.array_equal(atoms.constraints[0].get_indices(), np.flatnonzero(atoms.get_tags() == 1))
    assert np.all(atoms.get_tags()[192:] == 0)
    assert np.allclose(atoms.cell.lengths()[:2], [np.hypot(3.784, 9.515), 4 * 3.784])
    i, j = neighbor_list("ij", slab, 2.4)
    opposite = slab.numbers[i] != slab.numbers[j]
    coordination = np.bincount(i[opposite], minlength=len(slab))
    assert set(coordination[slab.numbers == 22]) == {5, 6}
    assert set(coordination[slab.numbers == 8]) == {2, 3}
    occupied = []
    for n in range(waters):
        oxygen = 192 + 3 * n
        distances = atoms.get_distances(oxygen, np.arange(192), mic=True)
        ti = np.flatnonzero(slab.numbers == 22)
        nearest_ti = ti[np.argmin(distances[ti])]
        occupied.append(nearest_ti)
        assert coordination[nearest_ti] == 5
        assert distances[nearest_ti] == pytest.approx(2.3)
        assert atoms.get_distances(oxygen, [oxygen + 1, oxygen + 2], mic=True) == pytest.approx([0.9572] * 2)
    assert len(set(occupied)) == waters
    oxygen_sites = np.arange(192, len(atoms), 3)
    for i, oxygen in enumerate(oxygen_sites[:-1]):
        assert np.min(atoms.get_distances(oxygen, oxygen_sites[i + 1 :], mic=True)) >= 3.784 - 1e-8
    if waters == 2:
        assert atoms.get_distance(192, 195, mic=True) == pytest.approx(3.784)


def test_adsorption_requires_bound_molecular_water():
    atoms = geometry(2)
    assert adsorption_status(atoms, 2)["passed"]
    # A fully intact but detached molecule must fail even if its forces vanish.
    atoms.positions[195:, 2] += 2
    status = adsorption_status(atoms, 2)
    assert not status["passed"]
    assert status["water_geometry"][1]["molecular"]
    assert not status["water_geometry"][1]["ti_bound"]
    atoms = geometry(1)
    atoms.positions[193, 2] += 2
    status = adsorption_status(atoms, 1)
    assert not status["passed"]
    assert status["water_geometry"][0]["ti_bound"]
    assert not status["water_geometry"][0]["molecular"]


@pytest.mark.reference
@pytest.mark.parametrize("waters", [1, 2, 4, 8])
def test_monti_reference(waters, tmp_path):
    atoms = read(INPUTS / f"{waters}water_initial.xyz")
    assert verify(atoms, ForceField.bundled(FORCE_FIELD), tmp_path / "lammps")["passed"]
