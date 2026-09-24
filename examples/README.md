# Molecular and periodic examples

`water_cluster.py` runs entirely from the repository with NumPy and Autograd:

```sh
python examples/water_cluster.py
```

To also verify against the pinned `lmp_mpi` executable:

```sh
python examples/water_cluster.py --verify --output validation/runs/my-water-check
```

Choose a new output directory. The script saves structures in XYZ format,
both sets of numerical results, raw LAMMPS inputs/outputs, and a summary with
explicit tolerances. Nonmatching reference results cause a nonzero exit status.
The default data source is the bundled LAMMPS QEq water parameter file; it
does not depend on the ZnOH example or on any Zn-specific calculator logic.

The example uses a monomer and a distorted dimer. The same geometries support
single-point energy, force, charge, property, symmetry, and derivative checks.
These are unoptimized validation fixtures.

Verified properties are charges, dipole vectors (e Å), per-atom total bond
orders, lone-pair counts, and bond counts at the LAMMPS default threshold of
0.3. Force comparisons use the default fixed-charge `forces`
(`full_derivative=False`), matching LAMMPS. The example computes no
charge-response forces. Those are available through
`calc.evaluate(symbols, positions, full_derivative=True)` and are checked
separately by finite differences in the tests.

See [saved results](../validation/water/summary.json) and the
[package README](../README.md) for supported formats and limitations.

`ase_water.py` uses the optional `ReaxFFCalculator` adapter with ASE optimizers:

```sh
python examples/ase_water.py --verify
python examples/ase_water.py --optimizer BFGS --verify
```

It writes the initial/final XYZ structures, an ASE trajectory, optimizer log,
and results in a fresh directory. ASE energies/forces and optimizer tolerance
use eV and eV/Angstrom; the saved core `python.json` and reference comparison
use kcal/mol and kcal/mol/Angstrom. Both optimizers use fixed-charge forces
and equilibrate charges at every evaluated geometry.

`Calculator.relax()` uses ASE FIRE by default. The original FIRE implementation
remains available through `backend="native"` and in the native relaxation benchmark.

`periodic_water.py` verifies fixed-cell periodic calculations:

```sh
python examples/periodic_water.py --verify
python examples/periodic_water.py --verify --include-bulk
```

It covers a boundary-crossing dimer and a rotated triclinic water slab.
The optional `--include-bulk` adds the 192-atom box. Structures are deterministic
fixtures, not equilibrated liquid snapshots. See the
[water results](../validation/water/README.md).

`small_cells.py` demonstrates repeated images and verification against normalized
larger LAMMPS supercells:

```sh
python examples/small_cells.py --verify
```

Four cases cover the 4 Å water monomer, a small triclinic water slab, periodic
ZnO, and a zinc atom bonded to its own images. Both neighbor builders feed one
energy model; QEq solves for input-cell charges. The native builder uses
replicated search cells; ASE builds image-resolved neighbors directly.

The shared system validation runner combines these structures with the C/H/O
and Zn/O checks, retaining each structure only once:

```sh
python scripts/validate.py --verify
python scripts/validate.py --verify --system cho
```

See [validation](../validation/README.md) for the coverage and retained records.

`cuo_surface.py` adds one optional 96-atom CuO(010) slab. It always checks
energies, forces, charges, and properties against LAMMPS, and measures
single-CPU performance. The Cu/O parameter file is an external published
supplement; see [download instructions and saved results](../validation/cuo/README.md).

```sh
python examples/cuo_surface.py --ffield path/to/jp102272z_si_001.txt
```
