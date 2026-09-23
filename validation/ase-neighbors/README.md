# ASE neighbor-list verification

These are the retained results from before energy-model consolidation. They
remain an independent numerical regression target. See the
[single-model results](../unified-energy/README.md) for the current architecture.

The ASE calculator now defaults to `neighbor_backend="ase"`, using ASE's
image-resolved neighbor list on the input cell. `neighbor_backend="replicated"`
retains the native dense image search and small-cell replication. The core
calculator retains the latter as its default. Both use fixed-charge forces.
The ASE adapter now builds `(i, j, S)` before calling the core and passes them
as `evaluate(..., neighbors=(i, j, S))`. The core consumes supplied arrays without
ASE or a neighbor search. Autograd differentiates edge vectors through positions
while treating indices and image shifts as constants.

**45 single-point cases pass** comparison between the two backends and against
fresh `lmp_mpi` runs (LAMMPS 22 Jul 2025, Update 4). Small-cell LAMMPS references
use equivalent expanded cells normalized to the input cell, preserving the
previous hydrogen-bond verification convention.

| Quantity | ASE vs replication | ASE vs LAMMPS |
| --- | ---: | ---: |
| Total energy per atom, kcal/mol | 3.23e-13 | 4.29e-12 |
| Energy component per atom, kcal/mol | 2.03e-12 | 2.06e-11 |
| Force component, kcal/mol/Å | 1.97e-12 | 2.48e-8 |
| Charge, e | 1.71e-14 | 1.14e-11 |
| Dipole component, e Å | 1.43e-14 | 1.70e-10 |
| Per-atom total bond order | 2.67e-15 | 4.22e-15 |
| Lone-pair count | 1.78e-15 | 1.78e-15 |
| Bond count | Exact match | Exact match |

Values are maxima over the retained cases, rounded upward. The full bond-order
matrices also agree between backends within 1.12e-15. LAMMPS does not export a
matching full matrix through this harness; its per-atom properties are checked.
Every result satisfies the strict backend tolerances and the existing LAMMPS
tolerances recorded in [summary.json](summary.json).

Coverage includes isolated atoms; Zn/O clusters up to 200 atoms; water monomer,
dimer, distorted dimer, trimer, and hexamer; methane, CO, C2, and active carbon
torsions; boundary-crossing molecules; multiple images; rotated triclinic and
partially periodic cells; 192-atom bulk water; 3.12–9 Å small water cells; Zn/O
small cells; a one-atom zinc chain with self-image bonds; a carbon chain with
torsions across repeated images; and cutoff crossings at 5 and 10 Å. These
are implementation checks, not physically equilibrated structures.

The full suite passes **240 tests**, including 75 neighbor tests. Additional
checks cover both force derivatives by finite differences, primitive-sized QEq,
independent explicit enumeration of neighbor shifts, wrapping/rotation/atom
permutation, backend switching and cache invalidation, list rebuilding across
cutoffs, no replication on the ASE path, both inner-wall vdW variants, and
ASE/native fixed-charge relaxation verified with LAMMPS. Explicit-array tests
also guard the builder/core boundary, verify evaluation without importing ASE,
reject malformed/duplicate/half lists, and check reuse of supplied skin lists.

## Retained data

- [summary.json](summary.json): environment, force-field checksums, per-case
  errors/tolerances, timings, LAMMPS version, and pass/fail results.
- [structures.json](structures.json): all input structures, cells, and PBC flags.
- [results.json.gz](results.json.gz): full numerical input-cell results keyed
  by case, then `ase`, `replicated`, and `reference`. Python entries include
  full bond-order matrices; all entries include forces and charges.

The compressed JSON can be read with `json.load(gzip.open(path, "rt"))`.
Raw LAMMPS inputs and outputs remain locally under
`validation/runs/supplied-neighbors/` for the latest explicit-array rerun;
the compact numerical records above are
tracked by Git. Per-case timings are single evaluations, not a statistical
performance benchmark.

Reproduce in a new output directory:

```sh
python scripts/validate_neighbors.py --verify
python -m pytest -q
```

The validation script writes complete per-case inputs, both Python results,
and LAMMPS inputs/outputs. Its optional `--output` selects a new directory;
`--executable` selects the reference binary. Omit `--verify` to compare only
the two Python backends.
