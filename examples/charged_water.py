"""Demonstrate numerical net-charge support, not validated ionic chemistry.

Run: python examples/charged_water.py
Requires xreac[ase]. Shared structures live in scripts/validate.py.
"""

from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.units import kcal, mol

ROOT = Path(__file__).resolve().parents[1]
for path in ("src", "scripts", "tests"):
    sys.path.insert(0, str(ROOT / path))

from validate import charged_validation_cases
from xreac import Calculator, ForceField
from xreac.ase import ReaxFFCalculator


def main():
    for name, (filename, symbols, x, cell, pbc, charge) in charged_validation_cases().items():
        if cell is not None:
            continue
        ff = ForceField.bundled(filename)
        calc = Calculator(ff)
        result = calc.evaluate(symbols, x, total_charge=charge)
        atoms = Atoms(symbols, positions=x, calculator=ReaxFFCalculator(ff, total_charge=charge))
        np.testing.assert_allclose(atoms.get_forces(), result.forces * kcal / mol, atol=1e-10)
        relaxed = calc.relax(symbols, x, total_charge=charge)
        print(f"{name}: Q={result.charges.sum():+.6f} e, E={result.energy:.6f} kcal/mol")
        print(f"  relaxed Q={relaxed.evaluation.charges.sum():+.6f} e, converged={relaxed.converged}")
        if not relaxed.converged:
            raise RuntimeError(f"{name} relaxation did not converge")


if __name__ == "__main__":
    main()
