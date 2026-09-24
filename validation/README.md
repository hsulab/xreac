# Validation

Results are organized by chemical system. The default suite uses **16 structures**
for energy, forces, QEq charges, dipoles, bond properties, and neighbor-list
agreement. Each structure is stored once. The 192-atom water box is optional.
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
contains one native result per default case. The optional bulk case has fresh
LAMMPS verification but is outside the routine baseline suite.

Water and Zn/O also retain ASE/native FIRE relaxation on their shared monomer
and dimer. Water contains two targeted diagnostics and a [PDF report](water/report.pdf).
Large Zn/O sizes are confined to optional [performance benchmarks](zno/benchmarks.json).

Single-CPU comparisons with `lmp_mpi` are recorded for
[water](water/cpu_benchmark.json), [Zn/O](zno/cpu_benchmark.json), and
[C/H/O](cho/cpu_benchmark.json). Each record contains the baseline and optimized
runs, their source checksums, and the measured speedup. Reproduce a fresh run with:

```sh
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
