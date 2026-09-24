# Validation

Results are organized by chemical system. The default suite uses **16 structures**
for energy, forces, QEq charges, dipoles, bond properties, and neighbor-list
agreement. Each structure is stored once. The 192-atom water box and
128-atom wurtzite ZnO bulk cell are optional.
An additional optional [96-atom CuO(010) slab](cuo/README.md) checks transfer
of the performance improvements to a copper-oxide surface, always with fresh
LAMMPS verification and the bundled Cu/O/H/Cl parameters, separately licensed
under [CC BY-NC 4.0](../data/README.md).
Its [follow-up benchmark](cuo/README.md#further-single-cpu-optimizations)
records a further 1.40x speedup in three separate commits with LAMMPS verification.

| System | Default cases | Purpose |
| --- | ---: | --- |
| [Water](water/README.md) | 6 | Molecules, hydrogen bonds, boundary crossings, small cells, triclinic/partial periodicity |
| [Zn/O](zno/README.md) | 6 | Atom, dimer, oxygen many-body terms, cluster, periodic/self-image bonds |
| [C/H/O](cho/README.md) | 4 | Methane, CO, C2 correction, periodic torsions |

```sh
python scripts/validate.py --verify
python scripts/validate.py --verify --system water --include-bulk
python -m pytest -q
```

Use `--output PATH` for a new output directory. Fresh runs go under ignored
`validation/runs/`; only selected records belong in the system folders.
Omit `--verify` for a comparison of the two Python neighbor backends alone.
Reference checks require `lmp_mpi` (LAMMPS 22 Jul 2025, Update 4).

The three default systems have the same compact files:

- `structures.json`: coordinates, atom labels, cell, and PBC for each case.
- `summary.json`: tolerances, discrepancies, timings, versions, and checksums.
- `results.json.gz`: full `ase`, `replicated`, and `reference` numerical results.
- `reference.tar.gz`: complete LAMMPS inputs and outputs for those cases.
- `baseline.json.gz`: selected independent results from before energy-model
  consolidation, used by regression tests. These values were not regenerated.

Read compressed results with `json.load(gzip.open(path, "rt"))`. The baseline
contains one native result per default case. The optional bulk cases have fresh
LAMMPS verification but is outside the routine baseline suite.

Water and Zn/O also retain ASE/native FIRE relaxation on their shared monomer
and dimer. Water contains two targeted diagnostics and a [PDF report](water/report.pdf).
Historical Zn/O cluster sizes are retained in [performance benchmarks](zno/benchmarks.json).

## ASE-neighbor performance pilots

```sh
python scripts/benchmark_lammps.py --repeats 5 --batch-seconds 0.5 --lammps-calls 100
# Equivalent entry point:
python scripts/benchmark.py --repeats 5 --batch-seconds 0.5 --lammps-calls 100
```

The default runs exactly **three structures**. All are shared cases from
`scripts/validate.py`: the existing water box and CuO slab, plus one new ZnO
wurtzite bulk cell. Both bulk cells have full PBC; CuO has PBC in the surface
plane. Cell heights exceed the cutoff, so both codes use identical atom counts.
These are deterministic, unrelaxed fixtures, not equilibrium predictions.

The default xreac timer includes **fresh ASE neighbor construction on every call**
(`neighbor_skin=0`),
QEq, energy, fixed-charge forces, and all reported properties. It explicitly
calls `ReaxFFCalculator.calculate()` to bypass ASE's result cache. No neighbor
list is reused in this mode. Force-field loading, imports, and process startup are excluded.
Both programs use one numerical thread; LAMMPS uses one MPI rank. No OS core
affinity is imposed. The original bin-based ASE baseline, retained for comparison,
gave these medians over five batches on the Apple M1 Pro:

| Pilot | Atoms | xreac/ASE (ms) | Fresh LAMMPS (ms) | xreac / LAMMPS |
| --- | ---: | ---: | ---: | ---: |
| Bulk water | 192 | 535.94 | 45.45 | 11.79x |
| Bulk ZnO | 128 | 251.28 | 39.43 | 6.37x |
| CuO(010) surface | 96 | 117.10 | 8.23 | 14.24x |

All three pass fresh LAMMPS checks of energies, components, forces, charges,
dipoles, and bond properties. Timed final energies, charges, and forces also
pass. The LAMMPS fresh and steady timing modes have the meanings described
below; steady timings are retained separately. These ASE measurements are a
new baseline and do not replace or extend the historical native speedup tables.

Each system has `cpu_ase_pilot.json`, `cpu_ase_pilot_results.json.gz`, and
`cpu_ase_pilot_reference.tar.gz`:
[water](water/cpu_ase_pilot.json), [ZnO](zno/cpu_ase_pilot.json),
[CuO](cuo/cpu_ase_pilot.json). Coordinates remain in the existing system
`structures.json` files. Raw archives omit duplicate `ffield` files; restore
the named bundled file from `data/` before replaying their LAMMPS inputs.

`--skip-lammps-timing` retains mandatory LAMMPS numerical verification while
skipping the timing loops. The CuO example also defaults to ASE timing;
`--neighbor-backend replicated` explicitly reproduces historical CuO timings.
Optional bulk validation uses `--include-bulk`. The routine 16-case validation
suite remains unchanged.

### Tree construction and neighbor reuse

Two separate improvements use the same three structures:

- `3e68ee6`: use ASE's tree-based `PrimitiveNeighborList`, then expand its half
  list into directed arrays with NumPy. Small periodic boxes at a 10 Å cutoff
  gave the previous bin-based builder too few bins to prune pairs efficiently.
- `b2e35a8`: reuse ASE topology with a 0.3 Å per-atom skin. Rebuild after any atom
  moves more than the skin from its build position, or when cell, PBC, atom
  identity/count, cutoff, or parameters change. QEq and all physical results
  are recomputed on every call. Zero skin explicitly requests fresh builds.

Single-thread medians in milliseconds, five batches per mode:

| Pilot | Original ASE | Tree commit, fresh | Skin commit, fresh | Skin commit, reuse | Native, fresh | Total speedup, original / reuse |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Water, 192 atoms | 535.94 | 120.91 | 122.53 | 109.82 | 165.90 | 4.88x |
| ZnO, 128 atoms | 251.28 | 68.45 | 73.16 | 66.14 | 89.21 | 3.80x |
| CuO surface, 96 atoms | 117.10 | 29.28 | 30.85 | 27.62 | 30.68 | 4.24x |

The reuse and native columns were measured in the same run as the skin commit's
fresh column. Reuse gives another **1.11–1.12x** over fresh tree construction in
that run, and **1.11–1.51x** over fresh native construction. Timing variation
between runs explains part of the difference between the two fresh columns.
Reuse timings exclude the initial build and hold coordinates fixed; these are
steady throughput measurements, not moving-atom MD or relaxation speedups.
Actual savings during motion depend on how often rebuilding is necessary.
ASE result caching is bypassed; QEq has no history reuse in xreac.

LAMMPS remains faster. In the skin-commit run, fresh LAMMPS times were
44.73 / 38.79 / 8.08 ms (water / ZnO / CuO); its fixed-geometry reuse times were
17.25 / 15.55 / 4.38 ms. Against the latter, xreac/ASE reuse is still
6.37x / 4.25x / 6.31x slower. The LAMMPS modes also differ in setup and QEq
history, so these ratios do not isolate neighbor construction alone.

All initial evaluations pass LAMMPS checks of energies, components, forces,
charges, dipoles, and bond properties. Each pilot is also checked after a seeded
displacement of at most 0.03 Å (list reuse), then a 0.35 Å translation that forces
a rebuild. Both states match fresh native and LAMMPS evaluations. These are
perturbations of the three existing fixtures, not additional structure families.
The second commit passes all 182 tests, including pair-cutoff crossings,
cumulative displacement, periodic images, and cache invalidation.

Reproduce the comparison:

```sh
python scripts/benchmark_lammps.py --compare-neighbors --repeats 5 --batch-seconds 0.3 --lammps-calls 30
```

Per-commit timings, source hashes, rebuild counts, and numerical checks:
[water](water/cpu_ase_neighbors.json), [ZnO](zno/cpu_ase_neighbors.json),
[CuO](cuo/cpu_ase_neighbors.json). Each system also retains
`cpu_ase_neighbors_results.json.gz` and `cpu_ase_neighbors_reference.tar.gz`.
The archives contain initial, reuse, and rebuild reference inputs/outputs;
restore the bundled parameter file as `ffield` in each extracted directory.
`git_revision` records HEAD when a run started; `committed_revision` identifies
the exact implementation, verified against all recorded core source hashes.

## Bulk-water moving MD

The [Berendsen MD example and results](water/README.md#bulk-water-berendsen-md)
reuse the 192-atom water pilot for a 1 ps trajectory. Unlike the fixed-geometry
tables above, this comparison includes atom movement and neighbor rebuilding.
It also tests experimental QEq factorization/history reuse with a fresh matrix,
checked convergence, and direct-solve agreement at every step.

```sh
python examples/water_md.py --steps 4000 --warmup 100
python scripts/benchmark_water_qeq.py --steps 4000 --warmup 100
```

## Historical native-neighbor benchmarks

Single-CPU comparisons with `lmp_mpi` are recorded for
[water](water/cpu_benchmark.json), [Zn/O](zno/cpu_benchmark.json), and
[C/H/O](cho/cpu_benchmark.json). Each record contains the baseline and optimized
runs, their source checksums, and the measured speedup. These records used
native neighbors. Reproduce them with the driver at the recorded historical
commit (the current driver's default is ASE):

```sh
# At the recorded historical commit:
python scripts/benchmark_lammps.py --include-large
```

Both programs use one numerical thread, and LAMMPS uses one MPI rank. No OS
core affinity is imposed. The benchmark excludes imports, force-field loading,
process startup, and final dumps. It reports median wall time over five batches
for energy, charge equilibration, and fixed-charge forces. xreac uses the native
`Calculator.evaluate` API, including neighbor construction and all properties.
The structures and force fields are shared with validation; the larger Zn/O
clusters reuse the optional benchmark generator. All timed final energies,
forces, and charges are checked against xreac.

The **fresh** LAMMPS comparison rebuilds neighbors and recreates the QEq fix
for every `run 0`, clearing charge history. It also computes energy components,
bond properties, and dipole each call. Its timing includes command parsing and
run setup, so it is a conservative comparison with the stateless xreac API.
The **steady** comparison keeps the geometry fixed, reuses neighbors and QEq
history, and computes properties only at batch endpoints. This gives an
optimistic LAMMPS throughput baseline; it is not a moving-atom MD benchmark.
Each system's `cpu_benchmark.tar.gz` archives both runs under `baseline/` and
`optimized/`, together with the optimized source snapshot.

The optimization replaces per-edge Python parameter lookups with indexing
into small element-pair tables, constructs reverse-edge maps with array
sorting, and rejects unsupported torsion types before enumerating full atom
chains. The calculator obtains energy components and forces from one traced
energy pass. Fixed-charge evaluation still solves QEq outside differentiation;
charge-response evaluation differentiates through its QEq solve. No
geometry-dependent state is cached between evaluations, and the same force
field equations, cutoffs, and output properties are used.

On the Apple M1 Pro, with one numerical thread and five timed batches:

| Case | Atoms | Before (ms) | After (ms) | Speedup | After / fresh LAMMPS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Zn/O cluster | 20 | 13.77 | 5.36 | 2.6x | 21.3x |
| Zn/O cluster | 100 | 166.43 | 20.78 | 8.0x | 7.2x |
| Zn/O cluster | 200 | 455.51 | 52.73 | 8.6x | 5.6x |
| Periodic water | 192 | 725.82 | 242.60 | 3.0x | 5.2x |

Water monomer/dimer and methane timings are also retained. Their smaller
workloads improve by about 1.0–1.2x. The optimized larger cases remain about
12–14x slower than LAMMPS with neighbor and charge-history reuse.

The four implementation commits were also measured cumulatively in a fresh
run. Each row below adds only the named optimization; the step speedup is
relative to the preceding row. These are separate measurements from the
complete LAMMPS timing comparison above.

| Commit | Added optimization | Zn/O 200 (ms) | Step speedup | Water 192 (ms) | Step speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| `c21c74b` | Baseline with benchmark harness | 466.23 | — | 741.35 | — |
| `c65b199` | Array parameter lookup | 316.85 | 1.47x | 361.39 | 2.05x |
| `9e756fe` | Array neighbor bookkeeping | 285.22 | 1.11x | 281.47 | 1.28x |
| `33a53d4` | Early torsion filtering | 69.73 | 4.09x | 272.42 | 1.03x |
| `0f9e92f` | Single traced energy pass | 54.33 | 1.28x | 243.69 | 1.12x |

The cumulative speedups in this run are **8.58x** for Zn/O 200 and **3.04x**
for water 192. All stages passed the complete test suite and formatting check
before their implementation commit: 163 tests for parameter lookup and
neighbor bookkeeping, and 164 after adding the inactive-torsion regression.
All seven benchmark cases pass fresh LAMMPS verification at every stage.

Detailed per-commit timings, sample ranges, incremental and cumulative
speedups, source checksums, and numerical discrepancies are retained for
[water](water/cpu_commit_stages.json), [Zn/O](zno/cpu_commit_stages.json), and
[C/H/O](cho/cpu_commit_stages.json). Each system's `cpu_commit_stages.tar.gz`
contains the corresponding raw reference inputs and outputs. Reproduce the
measurement at any listed commit with:

```sh
# At the recorded historical commit:
python scripts/benchmark_lammps.py --include-large --skip-lammps-timing
```

This option times xreac and verifies each case with a LAMMPS single point; it
skips the long LAMMPS timing loops. Step speedups depend on the commit order,
and changes of only a few percent should be interpreted with timing variation
in mind.

Cutoff checks perturb the ZnO dimer around 5 and 10 Å. Symmetry, finite-difference,
QEq, cache, and optimizer checks reuse the same geometries. Distinct checks for
invalid inputs, inner-wall variants, and singular torsions remain in the tests.

The previous feature/version-based records and larger overlapping case lists
are available in Git history. Numerical tolerances are unchanged; see the
[verification guide](../docs/validation.md). These unoptimized fixtures verify
implementation consistency, not physical accuracy of a parameterization.
