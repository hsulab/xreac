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
The 192-atom water box, 128-atom wurtzite ZnO bulk cell, and 96-atom CuO(010)
surface are optional cases outside the routine suite. They form the three
performance pilots:

```sh
python scripts/benchmark_lammps.py
```

By default every timed xreac evaluation builds fresh ASE neighbors
(`neighbor_skin=0`) and computes QEq, energy, forces, and properties.
Use `--compare-neighbors` to additionally time native builds and ASE topology
reuse at fixed geometry, and check displaced copies of the same pilots against
LAMMPS on both a reuse and a rebuild. Reuse timings exclude the initial build;
they measure steady throughput, not moving MD. The driver calls `calculate()` explicitly
to bypass ASE's result cache at unchanged coordinates. All three pilots are
verified against LAMMPS using one MPI rank and one numerical thread. Bulk
ZnO is an unrelaxed wurtzite fixture with representative lattice parameters;
the water box is not an equilibrated liquid snapshot.

For moving bulk-water MD with ASE Berendsen and single-rank LAMMPS, run:

```sh
python examples/water_md.py --steps 1000 --warmup 100
```

This reuses the 192-atom water pilot at fixed volume, a 0.25 fs timestep,
300 K target, and a 100 fs thermostat time constant. Timings exclude warmup
and reference checks. Three samples from each trajectory are verified against
fresh LAMMPS evaluations. Matched masses, initial velocities, and temperature
degrees of freedom make the settings comparable; different thermostat ordering
means trajectories need not coincide. Increase `--steps` for longer runs.

| System | Retained results | Coverage |
| --- | --- | --- |
| Water | {download}`summary <../validation/water/summary.json>` | Monomer, distorted dimer, boundary crossing, rotated triclinic/partial PBC, small-cell hydrogen bonds, optional bulk |
| Zn/O | {download}`summary <../validation/zno/summary.json>` | Isolated atom, dimer, O-only many-body terms, 20-atom cluster, small cell, zinc self-image chain |
| C/H/O | {download}`summary <../validation/cho/summary.json>` | Methane, CO, C2 correction, periodic carbon torsions |

Reproduce all cases or select a system. The output path must be new:

```sh
python scripts/validate.py --verify
python scripts/validate.py --verify --system water --include-bulk
python scripts/validate.py --verify --system zno --include-bulk
```

Each system directory contains `structures.json`, `summary.json`, full numerical
results in `results.json.gz`, and complete LAMMPS inputs/outputs in
`reference.tar.gz`. The results are keyed by case and then `ase`, `replicated`,
or `reference`. Read them with `json.load(gzip.open(path, "rt"))`.
Fresh runs default to ignored `validation/runs/` directories. The tracked
system folders contain the selected reproducible records.

Download the water {download}`structures <../validation/water/structures.json>`,
{download}`results <../validation/water/results.json.gz>`, and
{download}`illustrated monomer/dimer report <../validation/water/diagnostics/report.pdf>`.
The {download}`water relaxation summary <../validation/water/diagnostics/relaxation.json>` and
{download}`Zn/O relaxation summary <../validation/zno/relaxation.json>` compare
ASE and native FIRE on the same monomer and dimer used by the single-point checks.

Cutoff tests perturb the ZnO dimer around 5 and 10 Å. Symmetry, finite-difference,
cache, neighbor-list reuse, and relaxation tests reuse the shared structures.
Historical larger Zn/O clusters are available with the benchmark's explicit
`--suite legacy --include-large` option; its timing backend is now ASE too.
Input validation and targeted edge-case tests remain separate.

These are deterministic verification geometries, not equilibrated liquid
snapshots or predictions of physical stability. Agreement establishes
implementation consistency, not agreement with quantum chemistry.

## Small-cell verification

The shared 4 Å water cell tests interactions between multiple images of the
same input atom. Python agrees with normalized larger LAMMPS supercells.
The pinned LAMMPS primitive-cell implementation excludes hydrogen-bond
acceptors sharing the donor's original atom ID. This omits a contribution of
-0.370052 kcal/mol/cell present in the equivalent supercell. QEq charges agree.
The retained {download}`hydrogen-bond diagnostic <../validation/water/diagnostics/hbond-images.json>`
uses only this cell; the full-cell/partial-cell symmetry tests reuse the same
small triclinic monomer. See the
[LAMMPS QEq restriction](https://docs.lammps.org/fix_qeq_reaxff.html#restrictions).

## QEq audit

The {download}`QEq audit <../validation/water/diagnostics/qeq.json>` reuses the isolated
monomer. It compares LAMMPS forces with finite differences of its energies,
solving QEq again at each displaced geometry. Controls repeat QEq at fixed
coordinates and change the initial neutral charge guess. The charge-response
difference persists and follows the conversion-constant mismatch described in
[force conventions](calculations.md#force-convention).

## Charged-system checks

`python scripts/validate.py --system water --include-charged` adds hydroxide,
hydronium, and charged variants of existing periodic water fixtures. It compares
neighbor backends and independently assembles the QEq matrix using ASE image
distances, eliminating one charge to enforce the total. Reports check charge
conservation, chemical-potential equality, and positive reduced curvature.
These are numerical checks, not validation against ionic experimental data.

`--verify` runs LAMMPS for both neutral and charged cases. Neutral cases use
fresh LAMMPS QEq. Charged cases supply the xreac charges and disable LAMMPS
QEq with `checkqeq no`, comparing energies, energy components, fixed-charge
forces, and properties. Charge agreement here checks preservation of the
supplied values; the independent constrained solve validates QEq itself.
The compact charged results and raw reference archive are retained under
`validation/water/charged/charged.json` and `validation/water/charged/charged_reference.tar.gz`. Unit tests also check both force conventions against
finite differences, charge propagation during relaxation, and supercell scaling.

## Hur C/H/O/Cl parameters

`python scripts/validate.py --system cho --include-chocl --verify` checks five
shared Hur-force-field fixtures. Neutral molecules use fresh LAMMPS QEq; the
charged SN2 geometries use supplied charges. The
{download}`results <../validation/cho/hur2021.json>` and
{download}`raw LAMMPS archive <../validation/cho/hur2021_reference.tar.gz>`
record numerical implementation checks, not a validated reaction barrier.
The {download}`source parameter pages and provenance <../validation/cho/hur2021_source.tar.gz>`
document extraction from the published supplementary information.
