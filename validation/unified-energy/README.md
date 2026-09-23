# Single energy model verification

There is one `EnergyModel`, in `src/xreac/energy.py`. Native replicated searches
and ASE neighbor lists both produce `(i, j, S)` for that model. The duplicate
energy implementation and expanded-cell energy wrapper have been removed.
Native replication is now only a neighbor-search operation: each selected
image maps back to an input index and lattice shift before energy evaluation.

All **45 cases** in [summary.json](summary.json) pass fresh `lmp_mpi` comparisons
(LAMMPS 22 Jul 2025, Update 4), using fixed-charge forces. Native and ASE builders
give exactly identical total energies, all 14 components, forces, charges,
dipoles, full bond-order matrices, total bond orders, lone pairs, and bond counts.
Tests separately verify equality of the complete `(i, j, S)` arrays, including
self images and repeated neighbors. The largest discrepancy from LAMMPS is
2.48e-8 kcal/mol/Å in a force component; all existing reference tolerances pass.
Small-cell references still use normalized larger LAMMPS supercells.

The full suite passes **240 tests**. All 45 cases also compare against the
[pre-refactor numerical archive](../ase-neighbors/results.json.gz), which
retains the old independent dense/replicated implementation's results. This
checks the shared physics against a fixed baseline in addition to checking
the builders against each other. Other tests cover finite-difference derivatives,
QEq size and force conventions, image identities, wrapping/rotation/permutation,
neighbor-list reuse, input validation, and both relaxation optimizers.

The adapted QEq audit passes both monomer and dimer controls and LAMMPS finite
differences. The benchmark script runs a 20-atom calculation and two native FIRE
steps as a smoke check; that capped relaxation is not a convergence claim.

## Retained results

- [summary.json](summary.json): per-case comparisons, tolerances, timings,
  environment, force-field checksums, and reference versions.
- [structures.json](structures.json): all input structures, cells, and PBC flags.
- [results.json.gz](results.json.gz): full numerical results keyed by case and
  then `ase`, `replicated`, or `reference`.

Read the compressed JSON with `json.load(gzip.open(path, "rt"))`.
Raw LAMMPS runs remain locally in `validation/runs/unified-energy/`, including
the successful audit in `qeq-audit-check/`. Previous archives remain unchanged.
`cell_repetitions` in current Python results describes neighbor-search copies;
the energy model always receives input-cell atoms. In LAMMPS results it still
describes the actual reference simulation expansion.

Reproduce in a new directory:

```sh
python scripts/validate_neighbors.py --verify
python -m pytest -q
```
