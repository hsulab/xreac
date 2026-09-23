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
    result = calc.evaluate(symbols, x, full_derivative=True)
    # Directional derivatives exercise all atoms with QEq re-solved each time.
    direction = np.random.default_rng(64).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    fd = (calc.evaluate(symbols, x + h * direction).energy - calc.evaluate(symbols, x - h * direction).energy) / (
        2 * h
    )
    assert fd == pytest.approx(-np.sum(result.forces * direction), abs=1e-5, rel=1e-6)
    assert abs(result.charges.sum()) < 1e-12
    np.testing.assert_allclose(result.forces.sum(axis=0), 0, atol=1e-10)
    assert result.energy == pytest.approx(sum(result.components.values()), abs=1e-12)


@pytest.mark.parametrize("full_derivative", [False, True])
def test_symmetries(calc, full_derivative):
    symbols, x = CASES["cube8"]
    x = np.array(x)
    rotation, _ = np.linalg.qr(np.random.default_rng(33).normal(size=(3, 3)))
    r = calc.evaluate(symbols, x, full_derivative=full_derivative)
    moved = calc.evaluate(symbols, x @ rotation + [4.2, -3.1, 7], full_derivative=full_derivative)
    assert moved.energy == pytest.approx(r.energy, abs=1e-10)
    np.testing.assert_allclose(moved.forces, r.forces @ rotation, atol=1e-9)
    order = [7, 3, 1, 6, 2, 5, 4, 0]
    perm = calc.evaluate([symbols[i] for i in order], x[order], full_derivative=full_derivative)
    assert perm.energy == pytest.approx(r.energy, abs=1e-10)
    np.testing.assert_allclose(perm.charges, r.charges[order], atol=1e-12)
    np.testing.assert_allclose(perm.forces, r.forces[order], atol=1e-9)


@pytest.mark.parametrize(
    "symbols,x,kwargs",
    [
        ([], [], {}),
        (["He"], [[0, 0, 0]], {}),
        (["He"], [[0, 0, 0]], {"cell": [4.0] * 3}),
        (["Zn", "O"], [[0, 0, 0], [0, 0, 0]], {}),
        (["Zn"], [[np.nan, 0, 0]], {}),
        (["Zn"], [[0, 0]], {}),
        (["Zn"], [[0, 0, 0]], {"total_charge": 1}),
        (["Zn"], [[0, 0, 0]], {"cell": np.eye(3)}),
        (["Zn"], [[0, 0, 0]], {"full_derivative": "false"}),
    ],
)
def test_invalid_input(calc, symbols, x, kwargs):
    with pytest.raises(ValueError):
        calc.evaluate(symbols, x, **kwargs)


@pytest.mark.parametrize("backend", ["ase", "native"])
def test_relaxation(calc, backend):
    symbols, x = ["Zn", "O"], [[0, 0, 0], [2.3, 0.1, 0.2]]
    initial = calc.evaluate(symbols, x)
    result = calc.relax(symbols, x, force_tolerance=1e-5, backend=backend)
    assert result.converged, result.message
    assert result.evaluation.full_derivative is False
    assert result.evaluation.force_convention == "fixed_charge"
    assert result.evaluation.energy < initial.energy
    assert np.max(abs(result.evaluation.forces)) <= 1e-5
    np.testing.assert_array_equal(result.evaluation.forces, calc.evaluate(symbols, result.positions).forces)
    short = calc.relax(symbols, x, force_tolerance=1e-12, max_iterations=1, backend=backend)
    assert not short.converged
    assert short.iterations == 1
    assert short.evaluation.full_derivative is False
    assert np.max(abs(short.evaluation.forces)) > 1e-12
    np.testing.assert_array_equal(x, [[0, 0, 0], [2.3, 0.1, 0.2]])


@pytest.mark.parametrize("backend", ["ase", "native"])
def test_relaxation_already_converged(calc, backend):
    result = calc.relax(["Zn"], [[1.0, 2.0, 3.0]], max_iterations=1, backend=backend)
    assert result.converged and result.iterations == 0
    assert result.evaluation.full_derivative is False
    np.testing.assert_array_equal(result.positions, [[1.0, 2.0, 3.0]])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"force_tolerance": 0},
        {"force_tolerance": float("nan")},
        {"max_iterations": 0},
        {"max_iterations": 1.5},
        {"max_iterations": True},
        {"backend": "unknown"},
    ],
)
def test_invalid_relaxation_settings(calc, kwargs):
    with pytest.raises(ValueError):
        calc.relax(*CASES["zno"], **kwargs)


def test_charge_response_is_explicit(calc):
    fixed = calc.evaluate(*CASES["zno"])
    full = calc.evaluate(*CASES["zno"], full_derivative=True)
    assert fixed.full_derivative is False
    assert fixed.force_convention == "fixed_charge"
    assert full.full_derivative is True
    assert full.force_convention == "charge_response"
    assert fixed.energy == full.energy
    np.testing.assert_array_equal(fixed.charges, full.charges)
    assert np.max(abs(full.forces - fixed.forces)) > 0.03


@pytest.mark.parametrize("operation", ["evaluate", "relax", "relax_native"])
def test_default_does_not_differentiate_qeq(calc, monkeypatch, operation):
    from autograd.tracer import Box
    from xreac import energy

    original = energy.np.linalg.solve

    def solve_without_derivative(matrix, rhs):
        assert not isinstance(matrix, Box), "Default mode must not differentiate QEq"
        return original(matrix, rhs)

    monkeypatch.setattr(energy.np.linalg, "solve", solve_without_derivative)
    result = (
        calc.relax(*CASES["zno"], backend="native")
        if operation == "relax_native"
        else getattr(calc, operation)(*CASES["zno"])
    )
    if operation.startswith("relax"):
        assert result.converged
        result = result.evaluation
    assert np.isfinite(result.forces).all()


def test_fixed_charge_derivative(calc):
    from xreac.energy import EnergyModel
    from xreac.neighbors import replicated_neighbors

    symbols, positions = CASES["zno"]
    x = np.array(positions, dtype=float)
    result = calc.evaluate(symbols, x)
    neighbors, _ = replicated_neighbors(x, calc.force_field.general[12])
    model = EnergyModel(calc.force_field, symbols, neighbors)
    direction = np.random.default_rng(21).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    plus = model.components(x + h * direction, fixed_charges=result.charges)[0].sum()
    minus = model.components(x - h * direction, fixed_charges=result.charges)[0].sum()
    assert -(plus - minus) / (2 * h) == pytest.approx(np.sum(result.forces * direction), abs=1e-5, rel=0)


def test_parameter_reader(tmp_path):
    ff = ForceField.zno()
    assert ff.checksum == "b5af5a65573ce60f32405a973e33f7cc89a90bc69e624a8ca89811eb67b1a07b"
    assert ff.atoms["O"]["eta"] == 2 * 8.3122
    assert ff.pairs["Zn", "O"]["r_vdW"] == 2 * 2.1414
    assert ff.pairs["Zn", "O"]["r_p"] == pytest.approx((1.0548 - 1.6836) / 2)
    broken = tmp_path / "broken.ff"
    broken.write_text(ff.path.read_text()[:500])
    with pytest.raises(ValueError):
        ForceField.from_file(broken)
