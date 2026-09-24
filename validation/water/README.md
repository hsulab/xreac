# Water validation

All 7 retained cases pass fresh LAMMPS comparisons. The largest force
difference is 7.941e-10 kcal/mol/Å. Both Python neighbor builders agree.

| Case | Atoms | Checks |
| --- | ---: | --- |
| `cluster_water_monomer` | 3 | QEq, bond/angle forces, units, relaxation, charge-response audit |
| `cluster_water_distorted_dimer` | 6 | Hydrogen bonds, molecular properties, symmetry and derivatives |
| `periodic_boundary_dimer` | 6 | Boundary-crossing bonds, periodic forces, extensivity and neighbor rebuilding |
| `periodic_partial_pbc_water` | 24 | Rotated triclinic slab, partial PBC and dipole coordinate branches |
| `periodic_bulk_water_192` | 192 | Optional 64-molecule bulk-water size check |
| `small_water_4A` | 3 | Repeated images, self-image QEq and hydrogen-bond exclusions |
| `small_partial_pbc_water` | 3 | Small triclinic cell; full/partial PBC symmetry and derivative checks |

See [summary.json](summary.json), [structures.json](structures.json),
[full numerical results](results.json.gz), and [raw LAMMPS runs](reference.tar.gz).
The [independent baseline](baseline.json.gz) predates energy-model consolidation.

```sh
python scripts/validate.py --verify --system water
```

Use `--include-bulk` to include the optional 192-atom box.

The same box is the water performance pilot. Current benchmarks build fresh
ASE neighbors on every timed evaluation, bypass result caching, and verify
the result against LAMMPS. See [pilot timings](cpu_ase_pilot.json),
[numerical results](cpu_ase_pilot_results.json.gz), and
[raw reference runs](cpu_ase_pilot_reference.tar.gz).

```sh
python scripts/benchmark_lammps.py
```

- [Monomer/dimer PDF](report.pdf) and [CSV](report.csv).
- [QEq audit](qeq.json): the shared monomer, with fresh QEq at each finite-difference displacement.
- [Hydrogen-bond image diagnostic](hbond-images.json): the shared 4 Å cell versus its supercell.
- [Raw diagnostic runs](diagnostics.tar.gz).

```sh
python scripts/water_report.py
python scripts/audit_qeq.py
python scripts/audit_small_cells.py
```

[Relaxation](relaxation.json) compares ASE and native FIRE using the same
starting geometry as the single-point suite. [Raw optimizer results](relaxation.tar.gz)
include initial/final structures and LAMMPS verification.

```sh
python scripts/validate_relaxation.py --backend ase
python scripts/validate_relaxation.py --backend native
```
