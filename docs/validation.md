# Verification with LAMMPS

The reference executable is `lmp_mpi`, pinned to **LAMMPS 22 Jul 2025, Update 4**
with REAXFF support. The harness invokes one rank directly, without `mpirun` or
LAMMPS Python bindings. It writes zero initial charges and runs a single-point
`run 0` with active QEq, a tolerance of `1e-12`, and up to 2000 iterations.
Coordinates do not move during this calculation.

## Run a reference

```python
from xreac.reference import evaluate_lammps

reference = evaluate_lammps(ff, symbols, positions, cell=[4.0] * 3, directory="reference-water")
print(reference.energy)
print(reference.forces)
```

The directory must be new or empty. The harness retains the force field,
structure, script, dump, log, stdout/stderr, numerical outputs, executable
version, and parameter-file checksum. Missing executables, QEq failures, and
unexpected versions raise errors. Use the `executable` argument or
`XREAC_LAMMPS` environment variable to select an executable. A different version
requires an explicit `expected_version` override and revalidation.

For small cells, the harness runs an equivalent larger supercell and returns
results normalized to the input cell. Raw LAMMPS files refer to the expanded
system. `primitive.json` and metadata preserve the original structure,
replication factors, and energy divisor. Equivalent-copy forces and properties
are checked before averaging. Dipoles use the input coordinate branch.

The reference-only `allow_small_cell=True` option instead runs the raw primitive
cell for diagnostics. It is not the reference for validating small-cell support.

## Acceptance criteria

| Quantity | Maximum absolute difference |
| --- | --- |
| Total/component energy per atom | 1e-5 kcal/mol |
| Force component, fixed-charge convention | 1e-4 kcal/mol/Å |
| Charge | 1e-6 e |
| Dipole component on the same coordinate branch | 1e-6 e Å |
| Total bond order / lone-pair count | 1e-8 |
| Bond count | Exact match |

Tests also check symmetry, wrapping, extensivity, both force derivatives, ASE
integration, and both relaxation backends. Finite differences use fixed charges
when checking fixed-charge forces and fresh QEq when checking charge-response
forces. Numerical derivative checks are tests, not the production force method.

```sh
python -m pytest -q
python -m pytest -q -m 'not reference'
python -m pytest -q -m reference
```

## Results by system

The default suite reuses **16 representative structures**: six water cases,
six Zn/O cases, and four C/H/O cases. Each structure checks both neighbor
builders, all energy components, charges, forces, dipoles, and bond properties.
Independent pre-consolidation baselines are retained for those same structures.
The 192-atom water box is an optional seventeenth case.

| System | Retained results | Coverage |
| --- | --- | --- |
| Water | {download}`summary <../validation/water/summary.json>` | Monomer, distorted dimer, boundary crossing, rotated triclinic/partial PBC, small-cell hydrogen bonds, optional bulk |
| Zn/O | {download}`summary <../validation/zno/summary.json>` | Isolated atom, dimer, O-only many-body terms, 20-atom cluster, small cell, zinc self-image chain |
| C/H/O | {download}`summary <../validation/cho/summary.json>` | Methane, CO, C2 correction, periodic carbon torsions |

Reproduce all cases or select a system. The output path must be new:

```sh
python scripts/validate.py --verify
python scripts/validate.py --verify --system water --include-bulk
```

Each system directory contains `structures.json`, `summary.json`, full numerical
results in `results.json.gz`, and complete LAMMPS inputs/outputs in
`reference.tar.gz`. The results are keyed by case and then `ase`, `replicated`,
or `reference`. Read them with `json.load(gzip.open(path, "rt"))`.
Fresh runs default to ignored `validation/runs/` directories. The tracked
system folders contain the selected reproducible records.

Download the water {download}`structures <../validation/water/structures.json>`,
{download}`results <../validation/water/results.json.gz>`, and
{download}`illustrated monomer/dimer report <../validation/water/report.pdf>`.
The {download}`water relaxation summary <../validation/water/relaxation.json>` and
{download}`Zn/O relaxation summary <../validation/zno/relaxation.json>` compare
ASE and native FIRE on the same monomer and dimer used by the single-point checks.

Cutoff tests perturb the ZnO dimer around 5 and 10 Å. Symmetry, finite-difference,
cache, neighbor-list reuse, and relaxation tests reuse the shared structures.
Large Zn/O clusters belong to the optional benchmark script, not the routine
reference suite. Input validation and targeted edge-case tests remain separate.

These are deterministic verification geometries, not equilibrated liquid
snapshots or predictions of physical stability. Agreement establishes
implementation consistency, not agreement with quantum chemistry.

## Small-cell verification

The shared 4 Å water cell tests interactions between multiple images of the
same input atom. Python agrees with normalized larger LAMMPS supercells.
The pinned LAMMPS primitive-cell implementation excludes hydrogen-bond
acceptors sharing the donor's original atom ID. This omits a contribution of
-0.370052 kcal/mol/cell present in the equivalent supercell. QEq charges agree.
The retained {download}`hydrogen-bond diagnostic <../validation/water/hbond-images.json>`
uses only this cell; the full-cell/partial-cell symmetry tests reuse the same
small triclinic monomer. See the
[LAMMPS QEq restriction](https://docs.lammps.org/fix_qeq_reaxff.html#restrictions).

## QEq audit

The {download}`QEq audit <../validation/water/qeq.json>` reuses the isolated
monomer. It compares LAMMPS forces with finite differences of its energies,
solving QEq again at each displaced geometry. Controls repeat QEq at fixed
coordinates and change the initial neutral charge guess. The charge-response
difference persists and follows the conversion-constant mismatch described in
[force conventions](calculations.md#force-convention).
