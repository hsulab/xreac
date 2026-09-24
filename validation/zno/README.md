# Zn/O validation

All 6 routine cases pass fresh LAMMPS comparisons. The largest force
difference is 1.505e-10 kcal/mol/Å. Both Python neighbor builders agree.

| Case | Atoms | Checks |
| --- | ---: | --- |
| `zno_zn_atom` | 1 | Isolated atom and empty neighbor list |
| `zno_zno` | 2 | Dimer, force conventions, cutoff perturbations and relaxation |
| `zno_o4` | 4 | O-only angle, conjugation and torsion terms |
| `zno_cluster20` | 20 | Mixed-element coordination and cluster symmetries |
| `small_zno_4A` | 2 | Periodic mixed-element bonds |
| `small_zinc_chain` | 1 | Bonding to nonzero images of the same atom |

See [summary.json](summary.json), [structures.json](structures.json),
[full numerical results](results.json.gz), and [raw LAMMPS runs](reference.tar.gz).
The [independent baseline](baseline.json.gz) predates energy-model consolidation.

```sh
python scripts/validate.py --verify --system zno
```

[Relaxation](relaxation.json) compares ASE and native FIRE using the same
starting geometry as the single-point suite. [Raw optimizer results](relaxation.tar.gz)
include initial/final structures and LAMMPS verification.

```sh
python scripts/validate_relaxation.py --backend ase
python scripts/validate_relaxation.py --backend native
```

[Archived benchmarks](benchmarks.json) retain historical 20/100/200-atom timings.
Their iteration-capped native FIRE timings are not convergence claims.

The current performance pilot is **one 128-atom wurtzite ZnO bulk cell**
(64 Zn and 64 O), defined in `scripts/validate.py`. It uses an orthorhombic
wurtzite cell with representative a=3.25 Angstrom, c=5.21 Angstrom, u=0.382,
repeated (4, 2, 2) and perturbed by seeded 0.01 Angstrom displacements.
Its dimensions are approximately 13.00 x 11.26 x 10.42 Angstrom, all larger
than the 10 Angstrom cutoff. It is an unrelaxed fixture, not an equilibrium
prediction. Existing isolated clusters cannot exercise bulk periodicity.

See [ASE pilot timings and LAMMPS checks](cpu_ase_pilot.json),
[numerical results](cpu_ase_pilot_results.json.gz), and
[raw reference runs](cpu_ase_pilot_reference.tar.gz). Coordinates are stored
under `periodic_bulk_zno_128` in [structures.json](structures.json).
The newer [tree/reuse comparison](cpu_ase_neighbors.json),
[numerical results](cpu_ase_neighbors_results.json.gz), and
[reference archive](cpu_ase_neighbors_reference.tar.gz) record both optimization
commits and LAMMPS checks after movement, including a forced rebuild.

```sh
python scripts/benchmark.py
python scripts/benchmark.py --compare-neighbors
python scripts/validate.py --verify --system zno --include-bulk
```

`scripts/benchmark.py` now runs the same three ASE-neighbor pilots as
`scripts/benchmark_lammps.py`; the old cluster/FIRE driver is retained in Git
history through `b0ceb61`.
