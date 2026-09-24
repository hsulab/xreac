# Changelog

## 0.7.0

Improve ASE neighbor performance and add a verified bulk-water molecular
dynamics example. QEq remains fully converged at each geometry.

- Use ASE's tree-based neighbor builder and vectorized reverse edges. Add
  `neighbor_skin=0.3` Angstrom per atom by default, reusing topology until
  displacement, cell, periodicity, atom identity/count, cutoff, or calculator
  settings require rebuilding. Set `neighbor_skin=0` for fresh construction
  on every evaluation. Charges, energies, and forces are always recomputed
  at changed geometries.
- Across the three single-CPU pilots, tree construction improves full
  evaluation throughput by 3.7–4.4x over the original ASE bin-based path;
  topology reuse adds about 1.11x in the same-run comparisons. Reuse timings
  hold geometry fixed and exclude the initial build; moving-MD gains depend
  on rebuild frequency. Retain source hashes and numerical checks per commit.
- Add `examples/water_md.py`: ASE Berendsen MD for the existing 192-atom
  water box, with matched masses, initial velocities, and settings in LAMMPS.
  A 1 ps run gives 237.8 ms/step for xreac versus 31.9 ms/step for LAMMPS
  on one CPU thread, with substantial within-run timing variation. Sampled
  states from both trajectories pass fresh LAMMPS checks. The initial box is
  an unrelaxed fixture, not an equilibrated-water production model.
- Add an isolated QEq factorization/history-cache experiment. All 4,101
  charge comparisons pass, but cached solves take 2.50 s versus 1.74 s for
  direct solves over 4,000 timed steps. The direct solve accounts for only
  about 0.18% of MD time in this case; keep the production direct solver.

- Use fresh ASE neighbor lists inside timed benchmark evaluations, explicitly
  bypassing ASE result caching. Make bulk water (192 atoms), wurtzite ZnO bulk
  (128 atoms), and the CuO(010) surface (96 atoms) the three default pilots,
  with mandatory single-CPU LAMMPS verification. Share both bulk fixtures with
  optional validation and retain older native timings as historical records.
- Make `scripts/benchmark.py` an alias for the pilot benchmark driver and
  switch the CuO example's timing default to ASE neighbors.
- Reduce single-CPU CuO evaluation time from 42.52 ms to 30.42 ms (1.40x)
  with three separate optimizations: reuse bond-order values for properties,
  restrict bond-order work to the existing short-range cutoff, and pack
  neighbor sorting keys with an overflow-safe fallback. The same 96-atom
  surface gives identical results at every stage and passes LAMMPS checks.
  Retain per-commit timings under `validation/cuo/cpu_followup.json`.
- Bundle the unmodified Cu/O/H/Cl supplement as `data/ffield.reax.CuOHCl.2010`
  under its separate CC BY-NC 4.0 license. Include attribution, the complete
  license text, and file-specific license notes in source and wheel distributions.
- Use the bundled parameters by default in the CuO example and reference test;
  `--system cuo --verify` runs the existing slab without a separate download.

Validation: all 186 tests pass, including LAMMPS reference checks. The sampled
water-MD force differences are below 7e-10 kcal/mol/Angstrom. Formatting,
documentation, and package checks are part of the release workflow.

Licensing: code remains GPL-2.0-or-later. The bundled Cu/O/H/Cl parameter
file is separately CC BY-NC 4.0, including its noncommercial restriction;
see `NOTICE`, `LICENSES/CC-BY-NC-4.0.txt`, and `data/README.md`.

## 0.6.1

Improve single-CPU performance and correct a bond-order parameter selection
exposed by the new CuO surface example. The public calculator API is unchanged.

- Add one optional CuO(010) surface example with mandatory LAMMPS verification
  and a before/after performance comparison using an external published parameter file.
- Match LAMMPS's use of `valency_val` in uncorrected bond-order corrections,
  fixing Cu/O energies and forces when it differs from `valency_boc`.
- Speed up evaluations with array-based parameter lookup and neighbor
  bookkeeping, early filtering of unsupported torsions, and a single energy
  pass for reported values and forces. Force conventions and parameter
  equations are unchanged.
- Add reproducible single-CPU comparisons with `lmp_mpi`, including retained
  timings before and after optimization and per-commit speedups under `validation/`.
- On the 96-atom CuO(010) fixture, reduce evaluation time from 235.41 ms to
  41.51 ms (5.67x) on one Apple M1 Pro CPU thread. Both versions include the
  same correctness fix and produce identical results. Fresh LAMMPS evaluations
  take 8.09 ms; timings exclude process startup.
- Include the CuO baseline correction patch in source distributions so the
  retained before/after benchmark can be reproduced.

Validation: all 165 tests pass with the optional external Cu/O parameter file
enabled, including LAMMPS reference checks. The single CuO surface agrees with
LAMMPS to within 4.6e-10 kcal/mol/Angstrom in force components. Cu/O parameters
remain an external download; see `validation/cuo/README.md`.

## 0.6.0

First public release of xreac, a NumPy and Autograd ReaxFF implementation for
small-to-medium neutral molecules, clusters, and periodic cells.

- Parameter-driven energies, QEq charges, forces, dipoles, and bond properties.
- ASE calculator and optimizers, with native FIRE retained for benchmarks.
- Fixed-cell periodic boundaries, including small cells and repeated images.
- One energy model accepting explicit `(i, j, S)` neighbor arrays; ASE and native
  neighbor builders produce matching results.
- Bundled `ffield.reax.HO.2015`, `ffield.reax.CHO.2008`, and
  `ffield.reax.ZnOH.2010` parameter files with original citations and contents.
- System-organized validation: 16 shared routine structures and an optional
  192-atom water box. All 163 regression tests and fresh comparisons with
  LAMMPS 22 Jul 2025, Update 4 pass locally.
- Markdown/Sphinx documentation and retained numerical records and water PDF.

Requires Python 3.10 or newer. Install the wheel attached to this release,
or install from the tag:

```sh
python -m pip install 'git+https://github.com/hsulab/xreac.git@v0.6.0'
```

The `ase` extra enables ASE integration; `dev` adds testing, documentation,
and formatting tools. Core dependencies are NumPy and HIPS Autograd.

Forces default to the LAMMPS fixed-charge convention, with fresh QEq at each
geometry. Full differentiation through QEq is an explicit single-point option.
Relaxation always uses fixed-charge forces. Net charge, stress, variable-cell
relaxation, and MD integration are not supported in this release.
