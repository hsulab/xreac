# Small-cell diagnostic

Tested September 23, 2026 with the pinned single-rank `lmp_mpi` (22 Jul 2025,
Update 4), using `mamba run -n catorch3`. All structures are fixed-geometry
single points with fresh QEq. **139 tests pass.**

Each primitive cell is repeated exactly until every cell height exceeds 10 Å.
The larger structure represents the same infinite periodic arrangement.
Supercell energies are divided by the number of copies; per-atom forces,
charges, and bond properties are compared across every copy. Dipoles are not
compared across these different wrapping branches.

## Results

The following differences compare **LAMMPS primitive cell minus LAMMPS
supercell**, with energy normalized per primitive cell. They do not compare
an unsupported Python primitive calculation with LAMMPS.

| Primitive structure | Side (Å) | Repetitions | Expanded atoms | Energy difference (kcal/mol/cell) | Maximum force difference (kcal/mol/Å) |
| --- | ---: | --- | ---: | ---: | ---: |
| Water monomer control | 12 | 1×1×1 | 3 | 0 | 0 |
| Water monomer | 9 | 2×2×2 | 24 | -2.84e-13 | 3.55e-12 |
| Water monomer | 6 | 2×2×2 | 24 | 8.33e-7 | 3.78e-6 |
| Water monomer | 4 | 3×3×3 | 81 | 0.370052 | 1.05241 |
| Rotated water monomer | 3.12 | 4×4×4 | 192 | 3.70351 | 2.61656 |
| Distorted water dimer | 6 | 2×2×2 | 48 | 3.32e-5 | 9.89e-5 |
| Zn/O two-atom cell | 4 | 3×3×3 | 54 | 2.97e-12 | 5.68e-13 |

All seven **supported supercells** pass Python/LAMMPS comparisons for energy,
every energy component, fixed-charge forces, charges, dipole, bond-order sums,
lone pairs, and bond counts. Their maximum force discrepancy is 5.2e-12
kcal/mol/Å. Small-cell QEq charges agree with the replicated-cell charges to
better than 8e-14 e in these runs. An independent dense periodic QEq sum,
including nonzero self images, agrees with both sets of charges within 8e-14 e.

The water energy differences above come from **hydrogen bonding**, not QEq.
The pinned [LAMMPS hydrogen-bond source](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/REAXFF/reaxff_hydrogen_bonds.cpp#L98)
requires the donor and acceptor to have different `orig_id` values. In a
one-water primitive cell, an acceptor oxygen in another periodic image has the
same original ID as the donor oxygen, so that hydrogen bond is excluded. After
explicit replication, these oxygens have different IDs and the term is present.
The 4 Å case has zero hydrogen-bond energy in the primitive representation but
-0.370052 kcal/mol per primitive cell in the replicated representation; other
energy components agree within 1e-9 kcal/mol per cell.

This does **not** override the [documented LAMMPS QEq small-cell restriction](https://docs.lammps.org/fix_qeq_reaxff.html#restrictions).
These are specific single-rank tests of the pinned executable. They establish
that QEq was not the source of the observed discrepancies here, rather than
establishing small-cell correctness for arbitrary systems, versions, MPI
decompositions, or accelerator styles.

## Implications for xreac

The main implementation work is an image-aware interaction representation:

1. Identify a neighbor by `(atom_index, lattice_shift)`. Different images of
   the same atom may be separate bonded neighbors; a single `(N, N)` bond-order
   matrix cannot represent them. The dense water and Zn/O cases each have two
   active bond images for at least one atom pair.
2. Choose the number of image layers from cutoff/cell heights, rather than
   relying on the current large-cell ±1 shell.
3. Include `i=j` at **nonzero** lattice shifts in vdW/Coulomb and in QEq's
   diagonal. Exclude only the zero-displacement self term. At 3.12 Å, each atom
   has 146 nonzero self images inside the 10 Å cutoff.
4. Build angles, torsions, coordination sums, and exclusions from those image
   identities, and accumulate image forces back onto their primitive atoms.
   To preserve equivalence under replication, a donor and its translated copy
   must be distinguished when deciding hydrogen-bond exclusions.

A neighbor list alone is insufficient if later code collapses image identities
or discards all `i=j` terms. Small-cell implementation should be checked against
equivalent **large supercells**, because primitive LAMMPS results can themselves
depend on the representation, as demonstrated above.

The public Python calculator still rejects unsupported cell sizes. Only the
development reference harness accepts `allow_small_cell=True`, and its metadata
marks whether a run is outside the validated cell range.

## Reproduce

```sh
mamba run -n catorch3 python scripts/audit_small_cells.py
```

Use `--output` with a new directory to retain another run. The default is a
timestamped directory under `validation/runs/`. The script saves the primitive
structure and replication factors, both raw LAMMPS runs, supported-supercell
Python results, independent image-QEq results, and [summary.json](summary.json).
It fails if a supported-supercell comparison or the independent QEq check fails.
An inconsistent primitive LAMMPS result is recorded as a diagnostic finding.

These synthetic structures are not equilibrated liquid or crystal structures.
The 3.12 Å monomer is rotated to avoid an exactly collinear active torsion;
all copies retain the same orientation and geometry.
