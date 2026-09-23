import numpy as np
import pytest

from xreac import Calculator, ForceField
from cases import CASES


@pytest.fixture(scope="module")
def calc():
    return Calculator(ForceField.zno())


@pytest.mark.parametrize("name", ["zno", "o2", "ozno", "znozn", "o3", "o4", "cube8"])
def test_gradient(calc, name):
    symbols, positions = CASES[name]
    x = np.array(positions, dtype=float)
    result = calc.evaluate(symbols, x)
    # Directional derivatives exercise all atoms with QEq re-solved each time.
    direction = np.random.default_rng(64).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    fd = (calc.evaluate(symbols, x+h*direction).energy-calc.evaluate(symbols, x-h*direction).energy)/(2*h)
    assert fd == pytest.approx(-np.sum(result.forces*direction), abs=1e-5, rel=1e-6)
    assert abs(result.charges.sum()) < 1e-12
    np.testing.assert_allclose(result.forces.sum(axis=0), 0, atol=1e-10)
    assert result.energy == pytest.approx(sum(result.components.values()), abs=1e-12)


def test_symmetries(calc):
    symbols, x = CASES["cube8"]
    x = np.array(x)
    rotation, _ = np.linalg.qr(np.random.default_rng(33).normal(size=(3, 3)))
    r = calc.evaluate(symbols, x)
    moved = calc.evaluate(symbols, x@rotation + [4.2, -3.1, 7])
    assert moved.energy == pytest.approx(r.energy, abs=1e-10)
    np.testing.assert_allclose(moved.forces, r.forces@rotation, atol=1e-9)
    order = [7, 3, 1, 6, 2, 5, 4, 0]
    perm = calc.evaluate([symbols[i] for i in order], x[order])
    assert perm.energy == pytest.approx(r.energy, abs=1e-10)
    np.testing.assert_allclose(perm.charges, r.charges[order], atol=1e-12)
    np.testing.assert_allclose(perm.forces, r.forces[order], atol=1e-9)


@pytest.mark.parametrize("symbols,x,kwargs", [
    ([], [], {}), (["He"], [[0, 0, 0]], {}),
    (["Zn", "O"], [[0, 0, 0], [0, 0, 0]], {}),
    (["Zn"], [[np.nan, 0, 0]], {}), (["Zn"], [[0, 0]], {}),
    (["Zn"], [[0, 0, 0]], {"total_charge": 1}),
    (["Zn"], [[0, 0, 0]], {"cell": np.eye(3)}),
])
def test_invalid_input(calc, symbols, x, kwargs):
    with pytest.raises(ValueError):
        calc.evaluate(symbols, x, **kwargs)


def test_relaxation(calc):
    symbols, x = ["Zn", "O"], [[0, 0, 0], [2.3, .1, .2]]
    initial = calc.evaluate(symbols, x)
    result = calc.relax(symbols, x, force_tolerance=1e-5)
    assert result.converged, result.message
    assert result.evaluation.energy < initial.energy
    assert np.max(abs(result.evaluation.forces)) <= 1e-5
    short = calc.relax(symbols, x, force_tolerance=1e-12, max_iterations=1)
    assert not short.converged


def test_charge_response_is_explicit(calc):
    result = calc.evaluate(*CASES["zno"])
    assert np.max(abs(result.forces-result.lammps_forces)) > .03


def test_parameter_reader(tmp_path):
    ff = ForceField.zno()
    assert ff.checksum == "b5af5a65573ce60f32405a973e33f7cc89a90bc69e624a8ca89811eb67b1a07b"
    assert ff.atoms["O"]["eta"] == 2*8.3122
    assert ff.pairs["Zn", "O"]["r_vdW"] == 2*2.1414
    assert ff.pairs["Zn", "O"]["r_p"] == pytest.approx((1.0548-1.6836)/2)
    broken = tmp_path / "broken.ff"
    broken.write_text(ff.path.read_text()[:500])
    with pytest.raises(ValueError):
        ForceField.from_file(broken)
