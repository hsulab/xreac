import numpy as np
import pytest
from ase import Atoms
from ase.units import kcal, mol

from periodic_water import periodic_cases
from water_cluster import comparison, water_cases
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator
from xreac.energy import EnergyModel
from xreac.neighbors import replicated_neighbors
from xreac.reference import evaluate_lammps

CASES = periodic_cases()


@pytest.fixture
def calc():
    return Calculator(ForceField.bundled("ffield.reax.HO.2015"))


@pytest.mark.parametrize("name", ["partial_pbc_water"])
def test_wrapping_translation_rotation_and_permutation(calc, name):
    symbols, x, cell, pbc = CASES[name]
    original = calc.evaluate(symbols, x, cell=cell, pbc=pbc)
    rng = np.random.default_rng(900)
    shifts = rng.integers(-3, 4, x.shape) * pbc
    wrapped = calc.evaluate(symbols, x + shifts @ cell, cell=cell, pbc=pbc)
    assert wrapped.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(wrapped.forces, original.forces, atol=2e-8, rtol=0)
    np.testing.assert_allclose(wrapped.charges, original.charges, atol=1e-11, rtol=0)
    np.testing.assert_allclose(wrapped.bond_orders, original.bond_orders, atol=1e-11, rtol=0)
    # Periodic dipoles depend on the chosen coordinate branch.
    np.testing.assert_allclose(
        wrapped.dipole - original.dipole, np.sum(original.charges[:, None] * (shifts @ cell), axis=0), atol=1e-9
    )
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    transformed = calc.evaluate(symbols, x @ rotation + [21.0, -5.0, 40.0], cell=cell @ rotation, pbc=pbc)
    assert transformed.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(transformed.forces, original.forces @ rotation, atol=2e-8, rtol=0)
    order = rng.permutation(len(x))
    permuted = calc.evaluate([symbols[i] for i in order], x[order], cell=cell, pbc=pbc)
    assert permuted.energy == pytest.approx(original.energy, abs=1e-9)
    np.testing.assert_allclose(permuted.forces, original.forces[order], atol=2e-8, rtol=0)
    np.testing.assert_allclose(original.forces.sum(axis=0), 0.0, atol=1e-9)


def test_supercell_extensivity(calc):
    symbols, x, cell, pbc = CASES["boundary_dimer"]
    original = calc.evaluate(symbols, x, cell=cell)
    supercell = cell.copy()
    supercell[0] *= 2
    doubled = calc.evaluate(symbols * 2, np.vstack((x, x + cell[0])), cell=supercell)
    assert doubled.energy == pytest.approx(2 * original.energy, abs=1e-9)
    for name in ("charges", "forces", "total_bond_orders", "lone_pairs", "bond_counts"):
        np.testing.assert_allclose(
            getattr(doubled, name), np.concatenate([getattr(original, name)] * 2), atol=1e-9, rtol=0
        )


@pytest.mark.parametrize("full_derivative", [False, True])
@pytest.mark.parametrize("name", ["boundary_dimer"])
def test_periodic_derivative(calc, name, full_derivative):
    symbols, x, cell, pbc = CASES[name]
    result = calc.evaluate(symbols, x, cell=cell, pbc=pbc, full_derivative=full_derivative)
    neighbors, _ = replicated_neighbors(x, calc.force_field.general[12], cell, pbc)
    model = EnergyModel(calc.force_field, symbols, neighbors, cell, pbc)
    charges = None if full_derivative else result.charges
    direction = np.random.default_rng(581).normal(size=x.shape)
    direction /= np.linalg.norm(direction)
    h = 1e-5
    plus = model.components(x + h * direction, charges)[0].sum()
    minus = model.components(x - h * direction, charges)[0].sum()
    assert (plus - minus) / (2 * h) == pytest.approx(-np.sum(result.forces * direction), abs=2e-5, rel=1e-6)


@pytest.mark.parametrize(
    "cell,pbc",
    [
        (None, True),
        ([12.0, 12.0, 0.0], True),
        (np.zeros((3, 3)), True),
        ([12.0, -12.0, 12.0], True),
        ([12.0, 12.0, float("nan")], True),
        ([12.0, 12.0], True),
        ([12.0] * 3, "true"),
        ([12.0] * 3, [True, False]),
        ([12.0] * 3, [1, 2, 0]),
    ],
)
def test_invalid_cells(calc, cell, pbc):
    with pytest.raises(ValueError):
        calc.evaluate(*water_cases()["monomer"], cell=cell, pbc=pbc)


def test_periodic_input_semantics(calc):
    symbols, x = water_cases()["monomer"]
    original = calc.evaluate(symbols, x)
    display = calc.evaluate(symbols, x, cell=[2.0, 2.0, 2.0], pbc=False)
    assert display.energy == original.energy
    implicit = calc.evaluate(symbols, x, cell=[12.0, 12.0, 12.0])
    explicit = calc.evaluate(symbols, x, cell=np.diag([12.0] * 3), pbc=True)
    assert explicit.energy == implicit.energy
    # Only periodic directions require cell heights larger than the cutoff.
    calc.evaluate(symbols, x, cell=[12.0, 12.0, 3.0], pbc=[True, True, False])
    with pytest.raises(ValueError, match="Coincident"):
        calc.evaluate(["O", "H"], [[0.0, 0, 0], [12.0, 0, 0]], cell=[12.0] * 3)


def test_ase_cell_and_pbc_cache(calc, monkeypatch):
    symbols, x, cell, pbc = CASES["boundary_dimer"]
    atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc, calculator=ReaxFFCalculator(calc.force_field))
    calls = []
    evaluate = atoms.calc.core.evaluate

    def counted(*args, **kwargs):
        calls.append(kwargs)
        return evaluate(*args, **kwargs)

    monkeypatch.setattr(atoms.calc.core, "evaluate", counted)
    original = atoms.get_potential_energy()
    atoms.get_charges()
    assert len(calls) == 1
    atoms.cell[0, 0] += 1
    assert atoms.get_potential_energy() != original
    assert len(calls) == 2
    atoms.pbc = [True, False, False]
    expected = calc.evaluate(symbols, x, cell=atoms.cell.array, pbc=atoms.pbc)
    np.testing.assert_allclose(atoms.get_forces(), expected.forces * kcal / mol, atol=1e-12)
    assert len(calls) == 3


@pytest.mark.reference
@pytest.mark.parametrize("backend", ["ase", "native"])
def test_periodic_relaxation(calc, backend, tmp_path, monkeypatch):
    from autograd.tracer import Box
    from xreac import energy

    solve = energy.np.linalg.solve

    def fixed_only(matrix, rhs):
        assert not isinstance(matrix, Box), "Relaxation must not differentiate through QEq"
        return solve(matrix, rhs)

    monkeypatch.setattr(energy.np.linalg, "solve", fixed_only)
    symbols, x = water_cases()["monomer"]
    cell = np.diag([12.0] * 3)
    x = (x + [11.6, 11.6, 11.6]) % 12
    relaxed = calc.relax(symbols, x, cell=cell, backend=backend, force_tolerance=1e-5)
    assert relaxed.converged
    assert not relaxed.evaluation.full_derivative
    ref = evaluate_lammps(calc.force_field, symbols, relaxed.positions, cell=cell, directory=tmp_path / backend)
    report = comparison(relaxed.evaluation, ref, len(x))
    assert report["passed"], report
    assert np.max(abs(ref.forces)) < 1.01e-5


@pytest.mark.reference
def test_reference_unwrapped_branch(calc, tmp_path):
    symbols, x, cell, pbc = CASES["partial_pbc_water"]
    x = x + [21.0, -15.0, 30.0]
    x = x + (np.random.default_rng(24).integers(-2, 3, x.shape) * pbc) @ cell
    result = calc.evaluate(symbols, x, cell=cell, pbc=pbc)
    ref = evaluate_lammps(calc.force_field, symbols, x, cell=cell, pbc=pbc, directory=tmp_path / "reference")
    report = comparison(result, ref, len(x))
    assert report["passed"], report
