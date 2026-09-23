# Molecular examples

`water_cluster.py` runs entirely from the repository with NumPy and Autograd:

```sh
python examples/water_cluster.py
```

To also verify against the pinned `lmp_mpi` executable:

```sh
python examples/water_cluster.py --verify --output validation/my-water-check
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
0.3. Force comparisons use `lammps_forces`; the primary `forces` include the
QEq response and are checked separately by finite differences in the tests.

See [saved results](../validation/water/summary.json) and the
[package README](../README.md) for supported formats and limitations.
