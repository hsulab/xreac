# CuO surface example

This check uses **one structure**: a stoichiometric 96-atom CuO(010) slab
(48 Cu and 48 O), built from monoclinic tenorite. It has two conventional-cell
layers, a 2x3 surface repeat, 12 Å vacuum on each side, and periodicity only
in the surface plane. A seeded 0.01 Å displacement removes exact geometric
degeneracies. It is an unrelaxed implementation/performance fixture, not a
prediction of a relaxed surface or surface energy.

The bulk cell follows [Yang et al., Phys. Rev. B 39, 4343 (1989)](https://doi.org/10.1103/PhysRevB.39.4343):
space group C2/c, a=4.6837 Å, b=3.4226 Å, c=5.1288 Å, beta=99.54°,
Cu at (1/4,1/4,0), and O at (0,0.416,1/4). Both periodic slab heights
exceed the 10 Å cutoff, so xreac and LAMMPS evaluate exactly the same 96 atoms.

The bundled `data/ffield.reax.CuOHCl.2010` parameter file is the unmodified supporting data for
[van Duin et al., J. Phys. Chem. A 114, 9507–9514 (2010)](https://doi.org/10.1021/jp102272z),
[supplement DOI](https://doi.org/10.1021/jp102272z.s001). It includes H/O/Cu/Cl/X;
this example uses Cu and O only. ACS Figshare lists the supplement under
CC BY-NC 4.0. This file retains that separate license, including its
noncommercial restriction; see [parameter license notes](../../data/README.md),
[full license](../../LICENSES/CC-BY-NC-4.0.txt), and
[provenance and checksum](provenance.json).

## Run the example

From the repository root, with the ASE extra installed:

```sh
python examples/cuo_surface.py
```

The timed backend is **ASE**, with neighbor construction included on every
call and result caching bypassed. The native backend is evaluated separately
for consistency. Current ASE timings are retained in
[cpu_ase_pilot.json](cpu_ase_pilot.json); the earlier timing sections below
used native neighbors and are historical baselines.

LAMMPS verification is mandatory. The example checks energy components,
forces, charges, dipoles, and bond properties for both native and ASE neighbor
builders. A mismatch exits with an error. It saves the slab as extended XYZ,
full numerical results, and raw LAMMPS files under ignored `validation/runs/`.
Use `--ffield path/to/ffield` to override the bundled file. The supplied
force field must support Cu/O; the retained results use SHA256
`8b1a57a6945be8b1d9df3b69df8f329dc32561393ed1f02221a5ef0c21e8c07e`.

The shared validation runner reuses this exact fixture:

```sh
python scripts/validate.py --system cuo --verify
python -m pytest -q tests/test_cuo.py
```

CuO remains optional in the shared validation runner and adds no structures
to its default set. The CuO reference test uses the bundled file by default;
`XREAC_CUO_FORCE_FIELD` can override it. Like the other reference tests, it
requires LAMMPS and is excluded by `pytest -m 'not reference'`.

## Initial optimization (v0.6.1)

The new case exposed a pre-existing distinction in LAMMPS's bond-order
correction: the uncorrected coordination uses `valency_val`, while the
corrected coordination uses `valency_boc`. Cu has values 4 and 1 respectively
in this parameter file. xreac previously used `valency_boc` for both.
The correction follows the pinned
[LAMMPS bond-order implementation](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/REAXFF/reaxff_bond_orders.cpp).

For a meaningful speed comparison, the **same correctness correction** is
applied to both the pre-optimization baseline (`c21c74b`) and the optimized
implementation. [The baseline patch](baseline_correction.patch) contains no
performance changes. Both implementations pass all reference checks and agree
with each other within the repository's strict numerical tolerances.

Single-CPU medians over five batches on the Apple M1 Pro, excluding startup:

| Implementation | Time per evaluation |
| --- | ---: |
| Baseline with valency correction | 235.41 ms |
| Optimized with the same correction | 41.51 ms |
| LAMMPS, fresh neighbors and QEq | 8.09 ms |
| LAMMPS, fixed geometry with neighbor/QEq-history reuse | 4.22 ms |

The optimization gives **5.67x** speedup. Optimized xreac remains **5.13x**
slower than fresh LAMMPS calls and **9.83x** slower than the fixed-geometry
reuse benchmark. The latter is not a moving-atom MD benchmark. Each LAMMPS
timing batch contains 100 calls; both xreac and LAMMPS use one numerical thread,
and LAMMPS uses one MPI rank. No core affinity is imposed.

Retained records:

- [Structure](structures.json), [validation summary](summary.json), and
  [full numerical results](results.json.gz).
- [Baseline/optimized timings, source hashes, and all consistency checks](cpu_benchmark.json).
- [Raw LAMMPS inputs and outputs](reference.tar.gz). The duplicate `ffield` is
  omitted; copy `data/ffield.reax.CuOHCl.2010` to `ffield` in each extracted run
  directory before replaying it.

To reproduce the historical baseline without changing the working tree:

```sh
mkdir -p validation/runs/cuo-baseline
git archive c21c74b src | tar -x -C validation/runs/cuo-baseline
patch -d validation/runs/cuo-baseline -p1 < validation/cuo/baseline_correction.patch
python examples/cuo_surface.py \
  --neighbor-backend replicated \
  --source-root validation/runs/cuo-baseline --skip-lammps-timing
```

`--skip-lammps-timing` still performs fresh LAMMPS numerical verification.

## Further single-CPU optimizations

The next round reuses the **same 96-atom slab**, parameters, and force convention.
Each stage is measured in a fresh process in both forward and reverse commit
order. Each process runs seven batches targeting 0.5 seconds per batch;
the table reports medians over all 14 batch times per evaluation, excluding
startup. One numerical thread and one LAMMPS MPI rank are used, without core
affinity. These measurements are independent of the earlier timing session.

| Commit | Change | Time per evaluation | Speedup over previous stage | Speedup over this round's baseline |
| --- | --- | ---: | ---: | ---: |
| `691483c` | Previous implementation | 42.52 ms | — | 1.000x |
| `082920b` | Reuse bond-order values for properties | 40.14 ms | 1.059x | 1.059x |
| `7008b33` | Limit bond-order work to short-range neighbors | 35.77 ms | 1.122x | 1.188x |
| `5b22be2` | Pack neighbor sorting keys | 30.42 ms | 1.176x | 1.397x |

This reduces evaluation time by **28.4%**. Fresh LAMMPS evaluations in the
same run take **8.19 ms**, leaving xreac **3.71x** slower on this fixture.

Bond-order reuse is confined to a single evaluation: properties consume the
numerical values already computed for that geometry's energy and forces.
Bond candidates are selected again at every geometry using the existing
5 Angstrom cutoff (capped by the force-field cutoff). Nonbonded terms retain
their full neighbor lists. Packed integer keys preserve the canonical edge
order and reverse-edge checks; large shifts fall back to the previous sorting
method if packing could overflow. No cross-geometry cache is introduced.

All stages and both sweeps have **zero differences** from the baseline in
energy, components, forces, charges, bond orders, coordination, lone pairs,
bond counts, and dipoles. Each run independently verifies native and ASE
neighbors against LAMMPS. Maximum differences from LAMMPS are
3.8e-14 kcal/mol/atom in energy, 4.6e-10 kcal/mol/Angstrom in forces, and
5.3e-12 e in charges. All **169 tests** pass, including four integer-topology
checks for ordinary/large shifts and signed/unsigned indices. These checks
add no simulation structures.

- [Per-commit timing batches, source hashes, environment, and all comparisons](cpu_followup.json)
- [Final numerical results](cpu_followup_results.json.gz)
- [Final raw LAMMPS inputs, outputs, and timing logs](cpu_followup_reference.tar.gz)

The raw archive omits the duplicate `ffield`; restore it from
`data/ffield.reax.CuOHCl.2010` as for the original archive above.

To measure an individual stage with the current example driver:

```sh
revision=5b22be2
source_dir="validation/runs/cuo-followup-source-$revision"
mkdir -p "$source_dir"
git archive "$revision" src | tar -x -C "$source_dir"
python examples/cuo_surface.py --source-root "$source_dir" \
  --neighbor-backend replicated --repeats 7 --batch-seconds 0.5 --skip-lammps-timing
```

Repeat for the four revisions in table order, then reverse order, to obtain
14 batches per stage. The JSON records each batch's call count and duration;
pool duration/count across both sweeps before taking the median. Omit
`--skip-lammps-timing` to measure LAMMPS as well.
