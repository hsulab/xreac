# Molecular and periodic examples

`water_cluster.py` runs entirely from the repository with NumPy and Autograd:

```sh
mamba run -n catorch3 python examples/water_cluster.py
```

To also verify against the pinned `lmp_mpi` executable:

```sh
mamba run -n catorch3 python examples/water_cluster.py --verify --output validation/my-water-check
```

Choose a new output directory. The script saves structures in XYZ format,
both sets of numerical results, raw LAMMPS inputs/outputs, and a summary with
explicit tolerances. Nonmatching reference results cause a nonzero exit status.
The default data source is the bundled LAMMPS QEq water parameter file; it
does not depend on the ZnOH example or on any Zn-specific calculator logic.

The example includes nonoptimized monomer, dimer, distorted dimer, trimer,
and hexamer geometries. Ring donor bonds are tilted by 4 degrees away from
exact O-H...O collinearity to avoid undefined dihedral derivatives in weak
intermolecular torsions. These are validation fixtures, not proposed
equilibrium structures or benchmark predictions of water-cluster stability.

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
mamba run -n catorch3 python examples/ase_water.py --verify
mamba run -n catorch3 python examples/ase_water.py --optimizer BFGS --verify
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
mamba run -n catorch3 python examples/periodic_water.py --verify
```

It covers wrapped molecules, multiple images within the interaction cutoffs,
rotated triclinic cells, partial periodicity, and a 192-atom bulk water box.
Structures are deterministic test fixtures, not equilibrated liquid snapshots.
The output includes cell/PBC metadata, extended XYZ structures, fixed-charge
results, timings, and complete LAMMPS reference runs. See
[saved results](../validation/periodic-water/summary.json). Periodic cell heights
must exceed 10 Å for the bundled parameter files. The dipole follows the
supplied coordinate branch; wrapping atoms can change it.
