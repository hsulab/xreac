import numpy as np
import pytest

from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps
from water_cluster import comparison, water_cases

WATER = water_cases()
CARBON = {
    "methane": (["C"]+["H"]*4, np.vstack(([0, 0, 0], .63*np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])))),
    "carbon_monoxide": (["C", "O"], [[0, 0, 0], [1.15, 0, 0]]),
    "carbon_dimer": (["C", "C"], [[0, 0, 0], [1.2, 0, 0]]),
    "carbon_torsion": (["C"]*4, [[0, 0, 0], [1.5, .1, 0], [2.3, 1.2, .2], [3.7, .9, .8]]),
}


def test_parameter_driven_elements():
    ff = ForceField.bundled("qeq_ff.water")
    assert ff.elements == ("H", "O", "X")
    assert "Zn" not in ff.atoms
    Calculator(ff).evaluate(*WATER["monomer"])
    with pytest.raises(ValueError, match="absent from the force field"):
        Calculator(ff).evaluate(["Zn"], [[0, 0, 0]])
    # Nonbonded parameters are mixed even for types without a bond record.
    assert ff.pairs["X", "O"]["gamma"] > 0
    assert ff.pairs["X", "O"]["De_s"] == 0


def test_wildcard_torsions_and_precedence(tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    np.testing.assert_array_equal(ff.torsions["H", "O", "O", "H"], [2.5, -4., .9, -2.5, -1.])
    np.testing.assert_array_equal(ff.torsions["X", "O", "O", "X"], [.5511, 25.415, 1.133, -5.1903, -1.])
    np.testing.assert_array_equal(ff.torsions["O", "H", "O", "H"], [0, .1, .02, -2.5415, 0])
    lines = ff.path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if "Nr of torsions" in line)+1
    # Putting explicit entries after wildcard defaults must give the same table.
    lines[start:start+6] = lines[start+3:start+6]+lines[start:start+3]
    path = tmp_path / "reordered.ff"
    path.write_text("\n".join(lines)+"\n")
    reordered = ForceField.from_file(path)
    assert ff.torsions.keys() == reordered.torsions.keys()
    for key in ff.torsions:
        np.testing.assert_array_equal(ff.torsions[key], reordered.torsions[key])


def test_standard_format_variants(tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    text = ff.path.read_text().replace("50.0000", "5.00000D+1", 1)
    text = text[:text.index("  1    ! Nr of hydrogen bonds")]
    file = tmp_path / "without_hb.ff"
    file.write_text(text)
    parsed = ForceField.from_file(file)
    np.testing.assert_array_equal(parsed.general, ff.general)
    assert parsed.hydrogen_bonds == {}
    assert parsed.atoms.keys() == ff.atoms.keys()
    file.write_text("")
    with pytest.raises(ValueError, match="Empty"):
        ForceField.from_file(file)
    with pytest.raises(ValueError, match="plain filename"):
        ForceField.bundled("../pyproject.toml")


def test_atom_labels_preserved(tmp_path):
    original = ForceField.bundled("ffield.reax.cho")
    lines = original.path.read_text().splitlines()
    row = next(i for i, line in enumerate(lines) if line.split() and line.split()[0] == "C")
    lines[row] = lines[row].replace("C", "c", 1)
    path = tmp_path / "lowercase.ff"
    path.write_text("\n".join(lines)+"\n")
    parsed = ForceField.from_file(path)
    assert "c" in parsed.elements and "C" not in parsed.elements
    x = [[0, 0, 0], [1.2, 0, 0]]
    expected = Calculator(original).evaluate(["C", "C"], x)
    actual = Calculator(parsed).evaluate(["c", "c"], x)
    assert actual.energy == pytest.approx(expected.energy, abs=1e-12)
    np.testing.assert_allclose(actual.forces, expected.forces, atol=1e-12)


@pytest.mark.parametrize("name", WATER)
def test_water_gradients_and_properties(name):
    calc = Calculator(ForceField.bundled("qeq_ff.water"))
    symbols, x = WATER[name]
    result = calc.evaluate(symbols, x, full_derivative=True)
    h = 1e-5
    direction = np.random.default_rng(410).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    fd = (calc.evaluate(symbols, x+h*direction).energy-calc.evaluate(symbols, x-h*direction).energy)/(2*h)
    assert fd == pytest.approx(-np.sum(result.forces*direction), abs=2e-5, rel=1e-6)
    assert abs(result.charges.sum()) < 1e-12
    np.testing.assert_allclose(result.total_bond_orders, result.bond_orders.sum(axis=1), atol=1e-14)
    np.testing.assert_allclose(result.bond_orders, result.bond_orders.T, atol=1e-14)
    np.testing.assert_allclose(np.diag(result.bond_orders), 0, atol=1e-14)
    np.testing.assert_allclose(result.dipole, (x*result.charges[:, None]).sum(axis=0), atol=1e-12)
    np.testing.assert_allclose(result.forces.sum(axis=0), 0, atol=1e-9)
    if len(x) > 3:
        assert result.components["hydrogen_bond"] < -1


def test_water_symmetries():
    calc = Calculator(ForceField.bundled("qeq_ff.water"))
    symbols, x = WATER["distorted_dimer"]
    result = calc.evaluate(symbols, x)
    rotation, _ = np.linalg.qr(np.random.default_rng(42).normal(size=(3, 3)))
    transformed = calc.evaluate(symbols, x@rotation+[8, -12, 3])
    assert transformed.energy == pytest.approx(result.energy, abs=1e-10)
    np.testing.assert_allclose(transformed.forces, result.forces@rotation, atol=1e-9)
    np.testing.assert_allclose(transformed.dipole, result.dipole@rotation, atol=1e-12)
    order = [4, 3, 0, 2, 1, 5]
    permuted = calc.evaluate([symbols[i] for i in order], x[order])
    assert permuted.energy == pytest.approx(result.energy, abs=1e-10)
    np.testing.assert_allclose(permuted.forces, result.forces[order], atol=1e-9)
    np.testing.assert_allclose(permuted.bond_orders, result.bond_orders[np.ix_(order, order)], atol=1e-13)


@pytest.mark.reference
@pytest.mark.parametrize("name", WATER)
def test_water_reference(name, tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    symbols, x = WATER[name]
    actual = Calculator(ff).evaluate(symbols, x)
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path/name)
    report = comparison(actual, reference, len(x))
    assert report["passed"], report


@pytest.mark.reference
@pytest.mark.parametrize("name", CARBON)
def test_carbon_reference(name, tmp_path):
    ff = ForceField.bundled("ffield.reax.cho")
    symbols, x = CARBON[name]
    actual = Calculator(ff).evaluate(symbols, x)
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path/name)
    report = comparison(actual, reference, len(x))
    assert report["passed"], report
    if name == "carbon_dimer":
        assert actual.components["lone_pair"] > 50  # exercises the C2 correction
    if name == "carbon_torsion":
        assert abs(actual.components["torsion"]) > .1


@pytest.mark.parametrize("name", CARBON)
def test_carbon_gradients(name):
    symbols, x = CARBON[name]
    x = np.array(x, dtype=float)
    calc = Calculator(ForceField.bundled("ffield.reax.cho"))
    result = calc.evaluate(symbols, x, full_derivative=True)
    direction = np.random.default_rng(508).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-6
    fd = (calc.evaluate(symbols, x+h*direction).energy-calc.evaluate(symbols, x-h*direction).energy)/(2*h)
    assert fd == pytest.approx(-np.sum(result.forces*direction), abs=2e-5, rel=1e-6)


def test_collinear_active_torsion_rejected():
    calc = Calculator(ForceField.bundled("ffield.reax.cho"))
    with pytest.raises(ValueError, match="Collinear atoms in an active torsion"):
        calc.evaluate(["C"]*4, [[0, 0, 0], [1.5, 0, 0], [3, 0, 0], [4.5, .3, 0]])


def altered_water_file(path, parameters):
    lines = ForceField.bundled("qeq_ff.water").path.read_text().splitlines()
    for symbol in ("H", "O", "X"):
        start = next(i for i, line in enumerate(lines) if line.split() and line.split()[0] == symbol)
        for index, value in parameters.items():
            line = start+index//8
            fields = lines[line].split()
            fields[index % 8 + (1 if index < 8 else 0)] = str(value)
            lines[line] = " ".join(fields)
    path.write_text("\n".join(lines)+"\n")
    return ForceField.from_file(path)


@pytest.mark.reference
@pytest.mark.parametrize("shield", [False, True])
def test_inner_wall_variants(shield, tmp_path):
    # These synthetic files test the equations against LAMMPS, not physical models.
    ff = altered_water_file(tmp_path/"core.ff", {9: 15. if shield else 0., 29: .8, 30: .2, 31: 3.})
    assert ff.vdw_type == (3 if shield else 2)
    symbols, x = WATER["dimer"]
    actual = Calculator(ff).evaluate(symbols, x)
    reference = evaluate_lammps(ff, symbols, x, directory=tmp_path/"reference")
    report = comparison(actual, reference, len(x))
    assert report["passed"], report


@pytest.mark.reference
def test_water_relaxation(tmp_path):
    ff = ForceField.bundled("qeq_ff.water")
    calc = Calculator(ff)
    symbols, x = WATER["monomer"]
    relaxed = calc.relax(symbols, x, force_tolerance=1e-5)
    assert relaxed.converged, relaxed.message
    assert relaxed.evaluation.energy < calc.evaluate(symbols, x).energy
    reference = evaluate_lammps(ff, symbols, relaxed.positions, directory=tmp_path/"relaxed")
    # Energy minimization uses the full derivative; reference checks use fixed-charge forces.
    fixed = calc.evaluate(symbols, relaxed.positions)
    report = comparison(fixed, reference, len(x))
    assert report["passed"], report
