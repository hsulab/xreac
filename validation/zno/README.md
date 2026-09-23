# Zn/O validation

All 6 retained cases pass fresh LAMMPS comparisons. The largest force
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

[Archived benchmarks](benchmarks.json) retain 20/100/200-atom timings.
These larger sizes are performance checks outside the default regression suite.
Iteration-capped native FIRE timings are not convergence claims.

```sh
python scripts/benchmark.py
python scripts/benchmark.py --relax-iterations 20
```
