# Getting started

## Installation

Python 3.10 or newer is required. From the repository root, install into the
existing development environment:

```sh
mamba run -n catorch3 python -m pip install -e '.[ase,test]'
```

No local virtual environment is needed. For a minimal installation from a
checkout, use `python -m pip install .`. Core dependencies are NumPy and HIPS
Autograd. ASE is optional: `.[ase]` or `.[relax]` enables the adapter and the
default relaxation backend. `.[report]` adds Matplotlib for report generation;
`.[docs]` adds Sphinx, MyST, the theme, and ASE for API documentation.

These installation commands use a local checkout; they do not assume a public
PyPI release or a particular Git hosting URL.

## A single-point calculation

```{testcode} quickstart
from xreac import Calculator, ForceField

ff = ForceField.bundled("qeq_ff.water")
calc = Calculator(ff)
symbols = ["O", "H", "H"]
positions = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0], [-0.243, 0.94, 0.0]]
result = calc.evaluate(symbols, positions)

assert result.forces.shape == (3, 3)
assert abs(result.charges.sum()) < 1e-12
assert result.full_derivative is False
assert np.isfinite(result.energy)
```

Use `result.energy`, `result.forces`, and `result.charges` to access the numerical
results. `result.components` contains the 14 energy-component entries corresponding
to LAMMPS `compute pair reaxff`, including zeros for inactive terms.

| Quantity | Core units | Shape |
| --- | --- | --- |
| Positions | Å | `(N, 3)` |
| Total energy / components | kcal/mol, per input cell for PBC | Scalar / dictionary |
| Forces | kcal/mol/Å | `(N, 3)` |
| Charges | Elementary charge | `(N,)` |
| Dipole | e Å | `(3,)` |
| Bond-order matrix | Dimensionless | `(N, N)` |

Calculations use float64 and require a neutral system. Atom labels must match
the parameter file. Use `ForceField.from_file("path/to/ffield")` for another
supported file; the bundled water file does not contain Zn parameters.

## Run the examples

From the repository root:

```sh
mamba run -n catorch3 python examples/water_cluster.py
mamba run -n catorch3 python examples/periodic_water.py
mamba run -n catorch3 python examples/small_cells.py
```

Add `--verify` to compare against the pinned `lmp_mpi` executable. These scripts
save structures, numerical results, and a summary in a fresh output directory;
verification also retains full LAMMPS inputs and outputs. Use `--output PATH`
to choose a new directory. The example scripts belong to the source checkout,
not the installed Python API.

See [verification](validation.md) for tolerances and the distinction between
implementation agreement and physical accuracy.
