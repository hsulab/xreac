# xreac

A lightweight Python implementation of ReaxFF using NumPy and Autograd for
small-to-medium neutral or charged molecules, clusters, and periodic cells. It computes
energies, QEq charges, forces, and bond properties from standard parameter files,
with optional ASE integration for geometry optimization and fixed-volume MD.

## Install

Requires Python 3.10 or newer. From a source checkout:

```sh
python -m pip install .
```

Use `'.[ase]'` for ASE support, or `-e '.[dev]'` for development tools,
tests, and documentation.

## Example

```python
from xreac import Calculator, ForceField

calc = Calculator(ForceField.bundled("ffield.reax.HO.2015"))
result = calc.evaluate(
    ["O", "H", "H"],
    [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0], [-0.243, 0.94, 0.0]],
)
print(result.energy)  # kcal/mol
print(result.forces)  # kcal/mol/Å
print(result.charges)  # elementary charge
```

Use `ForceField.from_file("path/to/ffield")` for your own parameters. Charges
are equilibrated at each geometry; forces use the LAMMPS fixed-charge convention.

See the [documentation](docs/index.md) for the API and supported features,
[examples](examples/README.md) for ASE and periodic systems, and
[validation results](validation/README.md) for comparisons with `lmp_mpi`.
LAMMPS is needed only for verification.

The code is GPL-2.0-or-later. Bundled `ffield.reax.CuOHCl.2010` and
`ffield.reax.CHOCl.2021` data are separately licensed under **CC BY-NC 4.0**
and **CC BY-NC 3.0**, respectively, including noncommercial restrictions.
See [parameter licenses](data/README.md), [LICENSE](LICENSE),
[code attribution](NOTICE), and [parameter attribution](data/NOTICE)
for the license scope, attribution, and parameter-file provenance.
