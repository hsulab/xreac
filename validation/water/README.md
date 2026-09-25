# Water validation

| Location | Retained records |
| --- | --- |
| Top level | Shared structures, numerical results, summary, independent baseline, and raw LAMMPS reference archive |
| [charged/](charged/) | Charged-water checks with supplied-charge LAMMPS comparisons |
| [diagnostics/](diagnostics/) | QEq, hydrogen-bond images, relaxation, and monomer/dimer reports |
| [performance/](performance/) | CPU benchmarks, ASE-neighbor pilots, and optimization stages |
| [md/](md/) | Bulk-water dynamics and QEq-cache measurements |
| [surfaces/](surfaces/) | Two-layer Pt(111)/Ni(111) NEB examples, energy curves, raw outputs, and force-field provenance |

The shared cases remain at the top level for regression tests and report
generation. Files were moved without changing numerical data or archive
contents. Paths inside historical JSON records and archives describe their
original runs; current links and reproduction commands are provided here.
Fresh calculations still go under ignored `validation/runs/`.

All 7 core water cases pass fresh LAMMPS comparisons. The largest force
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

## Charged systems

[Charged validation](charged/charged.json) retains four numerical checks: hydroxide
(−1 e), hydronium (+1 e), and −1 e variants of the existing small cubic and
partially periodic water cells. It includes structures, force-field checksums,
numerical results, neighbor-backend comparisons, and independent QEq checks.
All four pass charge conservation and chemical-potential checks, with positive
curvature for charge-conserving redistributions. The same run passed six neutral
LAMMPS cases. All four charged cases also pass LAMMPS energy, component,
fixed-charge force, and property comparisons using identical supplied charges
with LAMMPS QEq disabled. This does not independently validate the charges.
[Raw charged LAMMPS inputs and outputs](charged/charged_reference.tar.gz) retain the
replicated cells, supplied charges, logs, and metadata. These checks do not establish ionic chemical accuracy.

```sh
python scripts/validate.py --system water --include-charged --verify
python examples/charged_water.py
```

Unit tests additionally cover both force derivatives, isolated and periodic
relaxation, ASE charge initialization, center-of-mass dipoles, and supercell scaling.

## Performance and other records

The same box is the water performance pilot. Current benchmarks build fresh
ASE neighbors on every timed evaluation, bypass result caching, and verify
the result against LAMMPS. See [pilot timings](performance/cpu_ase_pilot.json),
[numerical results](performance/cpu_ase_pilot_results.json.gz), and
[raw reference runs](performance/cpu_ase_pilot_reference.tar.gz).
The original bin-based baseline is retained there. The newer
[tree/reuse comparison](performance/cpu_ase_neighbors.json),
[numerical results](performance/cpu_ase_neighbors_results.json.gz), and
[reference archive](performance/cpu_ase_neighbors_reference.tar.gz) record both optimization
commits and LAMMPS checks after movement, including a forced rebuild.

The [CPU benchmark](performance/cpu_benchmark.json) and
[optimization-stage record](performance/cpu_commit_stages.json) retain later
measurements; their matching raw archives are in the same directory.

```sh
python scripts/benchmark_lammps.py
python scripts/benchmark_lammps.py --compare-neighbors
```

- [Monomer/dimer PDF](diagnostics/report.pdf) and [CSV](diagnostics/report.csv).
- [QEq audit](diagnostics/qeq.json): the shared monomer, with fresh QEq at each finite-difference displacement.
- [Hydrogen-bond image diagnostic](diagnostics/hbond-images.json): the shared 4 Å cell versus its supercell.
- [Raw diagnostic runs](diagnostics/diagnostics.tar.gz).

```sh
python scripts/water_report.py
python scripts/audit_qeq.py
python scripts/audit_small_cells.py
```

[Relaxation](diagnostics/relaxation.json) compares ASE and native FIRE using the same
starting geometry as the single-point suite. [Raw optimizer results](diagnostics/relaxation.tar.gz)
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

- [MD timings, QEq measurements, settings, checks, and source hashes](md/md_berendsen.json).
- [Sampled numerical results](md/md_berendsen_results.json.gz).
- [Raw LAMMPS MD and single-point inputs/outputs](md/md_berendsen_reference.tar.gz).

The initial geometry remains in `structures.json`; MD velocities and sampled
coordinates are included in the raw inputs. Archives omit duplicate parameter
files: restore `data/ffield.reax.HO.2015` as `ffield` in each extracted reference
directory. `committed_revision` identifies the implementation verified against
all recorded source hashes. The example is commit `47add41`; the cache experiment
is `020d9ab` with its final-correction check in `2852301`.

## Two-layer Pt(111) and Ni(111) water dissociation

The retained examples use `data/ffield.reax.PtNiCHO.2026`, a Gai2016 derivative
with only O–H–Pt `p_val4` changed from 1.0000 to 1.0250. Ni/O/H parameters are
unchanged from the original. These are not Assowe2012 results. See
[data/README.md](../../data/README.md#2026-experimental-revision) and `NOTICE`
for provenance and CC BY-NC 4.0 terms. No DFT fitting is claimed.

Both models contain p(2x2), two-layer slabs (eight metal atoms plus O/H/H),
with the bottom four metal atoms fixed, seven NEB images, and total vacuum
padding of 15 Angstrom. The bottom layer is at z=2.5 Angstrom: the slab and
adsorbate have been translated together down 5 Angstrom, leaving more space
above the surface. Direct checks show unchanged starting energies and forces.
Only x/y are periodic. Fixed lattice parameters are Pt 3.95 and Ni 3.52 Angstrom.

Endpoints are independently relaxed to 0.02 eV/Angstrom before IDPP and
ordinary/climbing-image NEB at 0.05 eV/Angstrom. Optimization includes QEq
charge-response derivatives; LAMMPS comparisons use fixed-charge derivatives
at equilibrated charges. FIRE uses dtmax=0.1 for NEB.

| Quantity | Ni(111) | Pt(111) |
| --- | ---: | ---: |
| Reactant/product relaxation steps | 177 / 257 | 107 / 688 |
| Ordinary NEB / climbing-image steps | 77 / 57 | 1,462 / 348 |
| NEB time, excluding endpoints and checks | 7.01 s | 83.71 s |
| Forward barrier | 0.52115 eV | 0.98303 eV |
| Final CI-NEB maximum force | 0.04919 eV/Angstrom | 0.04889 eV/Angstrom |
| Final LAMMPS comparisons | All pass | Reactant and peak pass; product fails |

Both bands meet their force targets, also checked by independent force
reconstruction. Pt's product retains the LAMMPS equivalent-copy force
consistency failure; its raw output and failed status are preserved. Numerical
convergence does not establish DFT accuracy, a global adsorption minimum,
or a global minimum-energy path. Molecular intermediates may lie below the
chosen locally relaxed reactant. Slab convergence and saddle frequencies
have not been checked.

Earlier two-layer, 13-image barriers were 0.51873 eV for Ni and 0.98022 eV for
Pt. Their compact comparison is retained in the current JSON reports; the
superseded bands, plots, and raw runs have been removed. The seven-image
agreement is specific to these examples, not a general resolution guarantee.

- [Ni results and audit](surfaces/ni_2026_7images_lowered.json), [structures](surfaces/ni_2026_7images_lowered_band.extxyz), [raw archive](surfaces/ni_2026_7images_lowered_raw.tar.gz)
- [Pt results and audit](surfaces/pt_2026_7images_lowered.json), [structures](surfaces/pt_2026_7images_lowered_band.extxyz), [raw archive](surfaces/pt_2026_7images_lowered_raw.tar.gz)
- [Energy curves](surfaces/metal_2026_7images_lowered.png)
- [Compact parameter-tuning evidence and original cutoff diagnostic](surfaces/pt_2026_parameter_tuning.json)
- [Original Gai supplement and publisher metadata](surfaces/gai2016_source.tar.gz)

The raw archives contain full-precision final bands, inputs, outputs, logs,
parameter files, and script snapshots. Historical paths inside numerical reports
identify the original runs; those scratch directories are not retained.

```sh
python scripts/neb_pt_water.py --metal Pt --ffield data/ffield.reax.PtNiCHO.2026 \
  --layers 2 --vacuum 15 --bottom-height 2.5 --images 7 --neb-dtmax 0.1 --steps 2500 \
  --output validation/runs/pt-water-example
python scripts/neb_pt_water.py --metal Ni --ffield data/ffield.reax.PtNiCHO.2026 \
  --layers 2 --vacuum 15 --bottom-height 2.5 --images 7 --neb-dtmax 0.1 --steps 2500 \
  --output validation/runs/ni-water-example
```

Use fresh output directories. The script retains its historical defaults;
the explicit options above define the kept examples. Further usage is described
in [examples/README.md](../../examples/README.md#water-dissociation-on-pt111-and-ni111).
