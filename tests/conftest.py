"""Reuse each system's single-point results across property and reference checks."""

from functools import lru_cache
from types import SimpleNamespace

from ase import Atoms
import pytest

from validate import validation_cases
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator


@pytest.fixture(scope="session")
def case_results():
    cases = validation_cases()

    @lru_cache(maxsize=None)
    def evaluate(name):
        filename, symbols, x, cell, pbc = cases[name]
        ff = ForceField.bundled(filename)
        native = Calculator(ff).evaluate(symbols, x, cell=cell, pbc=pbc)
        atoms = Atoms(symbols, positions=x, cell=cell, pbc=pbc, calculator=ReaxFFCalculator(ff))
        atoms.get_forces()
        return SimpleNamespace(force_field=ff, native=native, ase=atoms.calc.evaluation)

    return evaluate
