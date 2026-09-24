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
The original bin-based baseline is retained there. The newer
[tree/reuse comparison](cpu_ase_neighbors.json),
[numerical results](cpu_ase_neighbors_results.json.gz), and
[reference archive](cpu_ase_neighbors_reference.tar.gz) record both optimization
commits and LAMMPS checks after movement, including a forced rebuild.

```sh
python scripts/benchmark_lammps.py
python scripts/benchmark_lammps.py --compare-neighbors
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

## Bulk-water Berendsen MD

The same 192-atom box is used by `examples/water_md.py`. A 4,000-step run at
0.25 fs/step measures **1 ps of moving MD**, after 100 untimed warmup steps.
The target is 300 K with a 100 fs thermostat time constant. Both codes use
the force-field masses, the same seeded velocities, 3N temperature degrees
of freedom, and a 0.6 Å pair buffer (ASE's per-atom skin is 0.3 Å).
ASE uses `NVTBerendsen`; LAMMPS uses `fix nve` plus `fix temp/berendsen`.

| Engine | Total timed seconds | Mean ms/step | ns/day |
| --- | ---: | ---: | ---: |
| xreac with ASE neighbors | 951.07 | 237.77 | 0.091 |
| LAMMPS, one MPI rank | 127.72 | 31.93 | 0.676 |

xreac is **7.45x slower** in this run. Timings include integration, charge
equilibration every step, and displacement-triggered neighbor rebuilds;
startup, warmup, and reference validation are excluded. All numerical
libraries use one thread, without OS CPU affinity. xreac rebuilds its ASE
neighbor list 192 times during the 4,000 production steps.
xreac computes all reported bond properties on every force call; the LAMMPS
MD loop uses its usual energy/force/charge evaluation, with property checks
outside the timer. This compares the normal MD paths of the two programs.

These are whole-run means, not best-case fixed-geometry timings. The first
and second ASE blocks averaged 137.13 and 338.41 ms/step, respectively;
process CPU/wall ratios were approximately 0.99 and 0.88. The within-run
variation limits the precision of throughput comparisons and is not evidence
of a QEq-cache benefit. The separate cache experiment pairs solvers on each
identical matrix to isolate their costs.

Three snapshots from each engine pass fresh LAMMPS checks of energies,
components, forces, charges, dipoles, and bond properties. Live MD forces and
charges also agree with fresh evaluations. The largest sampled force error
is below 7e-10 kcal/mol/Å. ASE's sampled temperatures after warmup, halfway,
and at the end are 593.7, 362.0, and 327.3 K. This unrelaxed starting fixture
settles during the run; 1 ps is not a converged equilibrium-water study.
ASE and LAMMPS apply their Berendsen velocity scaling at different points in
velocity Verlet, so geometry-matched results, not identical trajectories,
are the correctness criterion.

```sh
python examples/water_md.py --steps 4000 --warmup 100
python scripts/benchmark_water_qeq.py --steps 4000 --warmup 100
```

### QEq caching result

The separate 1 ps experiment uses cached LU factors as a preconditioner and
the previous two solutions as a linear predictor. It rebuilds the current
QEq matrix at every geometry and requires an absolute KKT residual of at most
1e-12. If six refinement corrections are insufficient, it refactors.
Every cached solution is also compared with a direct NumPy solve of the
identical matrix; solver order alternates to reduce ordering bias.

| Solver | Total over 4,000 timed steps | Median per solve |
| --- | ---: | ---: |
| Fresh NumPy direct solve | 1.736 s | 0.329 ms |
| Cached SciPy LU with refinement/history | 2.504 s | 0.427 ms |

**This cache does not improve speed for the 192-atom water example.** It uses
1.44x as much total solver time, or 1.30x by the per-solve median. Across the
initial evaluation, warmup, and production, it makes 1,034 LU builds and
23,423 refinement corrections. All 4,101 charge comparisons pass, with maximum
charge difference 2.23e-13 e and maximum KKT residual below 1e-12.
Three cached-trajectory samples also pass the fresh-LAMMPS checks.

The direct linear solve accounts for only about **0.18%** of estimated MD
time in the paired experiment. Eliminating it completely would therefore
save at most about 0.18%; this cache instead adds an estimated **0.08%**.
Matrix assembly plus solving accounts for about **3.1%**. These fractions
subtract the extra cached-solver measurement from the experimental total;
the experiment's wall timer includes both solvers and charge-verification
overhead. They are cost estimates, not an independently measured production
speedup. Do not compare its whole-run wall time with the direct run to infer
a cache effect. The production calculator keeps its direct solver.

Retained files:

- [MD timings, QEq measurements, settings, checks, and source hashes](md_berendsen.json).
- [Sampled numerical results](md_berendsen_results.json.gz).
- [Raw LAMMPS MD and single-point inputs/outputs](md_berendsen_reference.tar.gz).

The initial geometry remains in `structures.json`; MD velocities and sampled
coordinates are included in the raw inputs. Archives omit duplicate parameter
files: restore `data/ffield.reax.HO.2015` as `ffield` in each extracted reference
directory. `committed_revision` identifies the implementation verified against
all recorded source hashes. The example is commit `47add41`; the cache experiment
is `020d9ab` with its final-correction check in `2852301`.
