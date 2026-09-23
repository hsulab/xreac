# Validation record

Verified locally on September 22, 2026 with LAMMPS 22 Jul 2025, Update 4,
invoked through `/opt/homebrew/bin/lmp_mpi`. The force-field checksums and
original citations are recorded in [NOTICE](../NOTICE).

Local environment names and Python executable paths have been omitted from
retained summaries. Package versions and numerical results are preserved.

## Single energy model

Both native replication and ASE now supply `(i, j, S)` to one `EnergyModel`.
The duplicated energy equations and the separate replicated-energy wrapper
have been removed. Native replication is used only to find neighbors, then
mapped back to input indices and lattice shifts. QEq, energy, and properties
are evaluated directly on the input atoms for both paths.

**240 tests pass**. All [45 retained comparisons](unified-energy/README.md)
give identical arrays and exactly identical numerical results between neighbor
builders. Fresh `lmp_mpi` checks also pass, with a maximum force discrepancy
of 2.48e-8 kcal/mol/Å. Regression tests compare with the pre-refactor numerical
archive below as well as testing neighbor-set equality. The QEq audit and native
FIRE benchmark script also run with the shared model.

## ASE neighbor lists (before model consolidation)

**240 tests pass**. ASE calculators now default to direct,
image-resolved ASE neighbors; native replication remains available explicitly
and remains the default for core single points and native FIRE.
ASE constructs `(i, j, S)` outside the core and passes `neighbors=(i, j, S)`
to `evaluate`. The core accepts these arrays from any builder without importing
ASE or calling neighbor search during evaluation or differentiation.

The [45-case comparison](ase-neighbors/README.md) passes both backend equivalence
and fresh LAMMPS checks. Maximum ASE/native differences are 3.23e-13 kcal/mol/atom
for total energy and 1.97e-12 kcal/mol/Å for forces. Full numerical results and
structures are retained with the report. Both backends include all periodic
images within their cutoffs, self-image bonds, and hydrogen-bond image identities.

```sh
python scripts/validate_neighbors.py --verify
```

## Small-cell support (0.6)

**163 tests pass**.

Small primitive cells now work through internal replication, with positions
and charges tied across equivalent copies. Energies and forces are returned
per input cell; QEq solves for only the input atom count. Nonzero self-image
interactions and hydrogen bonds between translated copies are retained.
`bond_orders` sums image contributions, including diagonal self-image bonds;
`bond_counts` applies its threshold to each image separately.

The nine [retained small-cell results](small-cell-support/summary.json) all pass
against normalized LAMMPS supercells. Coverage includes 3.12–9 Å water, a water
dimer, rotated triclinic and partially periodic cells, Zn/O, and a one-atom zinc
chain bonded to two images of itself. The maximum force difference is
4.2e-12 kcal/mol/Å; the maximum energy difference per atom is 4.3e-12 kcal/mol.

Verification also covers both force derivatives, reduced QEq matrix dimensions,
wrapping/rotation/permutation invariance, extensivity under different internal
replications, ASE caching and memory-limit changes, and fixed-charge relaxation
with both ASE and native FIRE. Relaxation convergence is checked on a
symmetry-controlled water orientation; a separate test checks iteration-limit
reporting for the general starting geometry. Weak image torques can make
orientation relaxation much slower than bond relaxation.

`cell_repetitions` reports the expansion. Automatic expansion defaults to a
512-atom limit, configurable with `max_expanded_atoms`. Pair work still scales
with the expanded system, while the QEq solve uses only the primitive atoms.
No dependencies were added. `evaluate_lammps()` now automatically expands small
cells and normalizes results, retaining the raw expanded run and its primitive
input/mapping metadata. Its diagnostic `allow_small_cell=True` mode remains
available to reproduce the original LAMMPS discrepancy below.

```sh
python examples/small_cells.py --verify
```

## Small-cell audit

The September 23, 2026 [small-cell audit](small-cells/README.md) compares seven
primitive cells with exactly replicated supercells using single-rank `lmp_mpi`.
All supported Python supercells match LAMMPS. The primitive water cells at
4 Å and 3.12 Å differ from their replicated representations in hydrogen-bond
energy and forces. QEq charges agree in all tested cases; an independent
explicit-image QEq sum confirms that QEq is not the discrepancy here.

The pinned LAMMPS hydrogen-bond code excludes donor/acceptor images sharing
the same original atom ID. The audit documents this representation dependence
and the image-aware neighbor/QEq work required for small-cell support in xreac.
This historical audit predates the internal-replication implementation above;
the calculator rejected small cells at that point. The reference-only
`allow_small_cell=True` switch is for diagnostics. **139 tests passed.**

## Periodic cells (0.5)

Verified September 23, 2026: **136 tests pass**. Fixed orthorhombic
and triclinic cells, rotated cell vectors, and partial periodicity are supported.
Periodic face heights must exceed 10 Å for the bundled files; all contributing
nonbonded and hydrogen-bond images are summed within their finite cutoffs.

The five [retained periodic water cases](periodic-water/summary.json) pass
total/component energy, fixed-charge force, QEq charge, dipole, bond-order,
lone-pair, and bond-count comparisons against `lmp_mpi`. They include
boundary-crossing molecules, multiple interacting images, rotated triclinic
and partially periodic cells, and 64 waters / 192 atoms in a 12.48 Å box.
These deterministic structures are not equilibrated liquid snapshots.

The largest discrepancies across the retained cases are approximately:

| Quantity | Maximum absolute difference |
| --- | ---: |
| Total energy per atom | 2.6e-12 kcal/mol |
| Energy component per atom | 2.1e-11 kcal/mol |
| Force component | 8.0e-10 kcal/mol/Å |
| Charge | 1.2e-11 e |
| Dipole component on the supplied coordinate branch | 1.7e-10 e Å |
| Per-atom total bond order | 3.6e-15 |
| Bond count | Exact match |

The tests additionally verify lattice wrapping, translation/rotation and atom
permutation invariance, supercell energy extensivity, derivatives under both
force conventions, periodic C/H/O and Zn/O cases, and ASE cell/PBC cache
invalidation. Both ASE FIRE and native FIRE converge a boundary-crossing water
molecule with no differentiation through QEq; final forces match LAMMPS.

The 192-atom evaluation took about 0.71 s in the retained run, including energy,
QEq, fixed-charge forces, and properties. This is a single local timing, not
an averaged benchmark. Results retain environment metadata, JSON/extended XYZ
structures, and full reference inputs and outputs. Periodic dipoles depend on
the input coordinate branch; the reference preserves it with image flags.

```sh
python examples/periodic_water.py --verify
python -m pytest -q
```

## ASE integration (0.4)

Validation used Python 3.10.13, NumPy 2.0.2, Autograd 1.9.1,
and ASE 3.27.0. **102 tests pass**, including adapter units, result caching,
parameter changes, unsupported inputs, ASE constraints, FIRE/BFGS integration,
and both relaxation backends. Tests confirm that neither relaxation backend
differentiates through QEq and that native FIRE does not import ASE.

`ReaxFFCalculator` exposes energy (eV), forces (eV/Å), charges (e), and dipole
(e Å) to ASE. `Calculator.relax()` reuses ASE FIRE by default; the original
implementation remains available with `backend="native"`. The convenience
method retains the max Cartesian component stopping criterion in kcal/mol/Å,
while direct ASE optimizers use max atomic vector norm in eV/Å.

Both backends pass final-geometry comparisons with `lmp_mpi`:

| System | ASE FIRE iterations | Native FIRE iterations | Force tolerance (kcal/mol/Å) |
| --- | ---: | ---: | ---: |
| ZnO dimer | 76 | 85 | 1e-5 |
| 20-atom Zn/O cluster | 205 | 209 | 1e-4 |
| Water monomer | 75 | 85 | 1e-5 |

See the retained [ASE results](relaxation-ase/summary.json) and
[native results](relaxation-native/summary.json). Direct ASE examples
also pass: [FIRE](ase-water-fire/summary.json) converges in 61 steps and
[BFGS](ase-water-bfgs/summary.json) in 6 steps, with a 1e-5 eV/Å vector-norm
tolerance. Their directories contain trajectories, optimizer logs, structures,
and raw LAMMPS verification inputs/outputs.

The [September 23 calculation benchmark](benchmark-20260923.json) records fixed-charge
evaluations only. The [native FIRE benchmark](benchmark-native-fire.json)
adds 20-step relaxation timings using our retained implementation. These runs
use a different interpreter/NumPy environment from the archived benchmarks;
do not attribute timing differences solely to code changes. Capped relaxation
timings do not imply convergence.

```sh
python -m pytest -q
python scripts/validate_relaxation.py --backend ase --output validation/runs/ase-check
python scripts/validate_relaxation.py --backend native --output validation/runs/native-check
python examples/ase_water.py --verify
```

## Fixed-charge relaxation (0.3.1)

`relax()` now uses FIRE with fixed-charge forces throughout, including its
returned evaluation. Charges are re-equilibrated at every geometry, and no
derivative through QEq is taken. FIRE needs no SciPy dependency. Convergence
is assessed using the largest absolute Cartesian force component.

**86 tests pass**, including a guard against differentiating the QEq solve
during relaxation, iteration-limit reporting, already-converged geometries,
and reference checks of the final forces. The retained September 23, 2026
[relaxation results](relaxation-fixed-charge/summary.json) all pass:

| System | FIRE iterations | Maximum force (kcal/mol/Å) | Requested tolerance |
| --- | ---: | ---: | ---: |
| ZnO dimer | 85 | 6.454e-6 | 1e-5 |
| 20-atom Zn/O cluster | 209 | 8.903e-5 | 1e-4 |
| Water monomer | 85 | 8.311e-6 | 1e-5 |

LAMMPS single-point checks at the final geometries confirm each force
tolerance. Maximum Python/LAMMPS force discrepancy is 2.0e-10 kcal/mol/Å.
These checks compare forces at the Python-relaxed geometries; they do not
assert identical optimization paths or minima from separate LAMMPS relaxations.
The earlier full-derivative relaxation records below remain historical results.

Reproduce in a new output directory:

```sh
python scripts/validate_relaxation.py --backend native --output validation/runs/relaxation-check
```

## Force selection (0.3)

The default `evaluate()` now returns fixed-charge forces, matching LAMMPS,
in `result.forces`. `full_derivative=True` selects charge-response forces
instead. Only one derivative is computed per evaluation. The separate
`lammps_forces` field is retired; archived 0.1/0.2 JSON files retain their
original field meanings.

**79 tests pass**, including both force modes, default-path checks that QEq
is not differentiated, and energy finite differences from fresh LAMMPS runs.
The new [water results](water-fixed-charge/summary.json) verify all five
cases with the default mode. The [updated report](water-fixed-charge/report.pdf)
uses only fixed-charge forces; optional charge-response entries are marked
as not evaluated. The archived report below includes both historical arrays.

The [new benchmark](benchmark-fixed-charge.json) measures only fixed-charge
evaluation and energy evaluation. It excludes relaxation and charge-response
calculations. The original benchmark below includes both force modes and
relaxation in peak RSS, so its memory figures are not directly comparable.

## General calculator and water validation (0.2)

**73 tests pass**, retaining all Zn/O regressions and adding parameter-driven
element selection, torsion-wildcard precedence, hydrogen bonding, C/H/O
systems, inner-wall vdW variants, water properties, and relaxation checks.
Synthetic inner-wall parameter modifications test mathematical compatibility;
they are not proposed physical force fields.

The dedicated water example uses LAMMPS's `ffield.reax.HO.2015` (Achtyl et al.,
2015), containing H/O and a dummy X type, with no Zn parameters. All five
retained cases pass: monomer, dimer, distorted dimer, trimer, and hexamer.
See [water/summary.json](water/summary.json) for tolerances and full results.

The [water PDF report](water/report.pdf) shows all five structures, both total
energies, signed energy differences (Python minus LAMMPS), and maximum/RMS
Cartesian force differences. It includes both matched-convention forces and
full energy-gradient forces, with per-atom plots and energy component tables.
The numerical overview is also available as [CSV](water/report.csv).
These use the retained results, which are identical to the subsequent
`water-20260923T032422Z` verification run. These geometries are not optimized;
energy differences in this report compare implementations, not binding energies.

Regenerate the report without rerunning calculations:

```sh
python -m pip install '.[ase]'
python scripts/water_report.py
# Archived two-force report:
python scripts/water_report.py --input validation/water --output validation/water/report.pdf
```

Matplotlib comes with the ASE extra and is not required by the core calculator. To report a new
verified run, pass `--input path/to/results --output path/to/report.pdf`.

| Compared quantity | Maximum absolute discrepancy |
| --- | --- |
| Total energy per atom | 7.6e-14 kcal/mol |
| Energy component per atom | 3.1e-13 kcal/mol |
| LAMMPS-convention force component | 8.2e-12 kcal/mol/Å |
| Atomic charge | 4.7e-14 e |
| Dipole vector component | 1.1e-12 e Å |
| Per-atom total bond order | 2.3e-16 |
| Lone-pair count | 6.4e-25 |
| Bond count (bond order > 0.3) | Exact match |

Each case retains its XYZ geometry, complete Python result, reference result,
and raw LAMMPS run files. Water-cluster directional finite differences verify
the full energy gradient independently. Rotation, translation, atom
permutation, zero net charge, and zero total force are also checked. A water
monomer relaxation converges and is re-evaluated against LAMMPS.

The separate Chenoweth C/H/O tests exercise methane, carbon monoxide,
the C2 correction, and a four-carbon torsion. These use an independent
parameter file and confirm the calculator does not require Zn/O parameters.
Exactly collinear active torsions are rejected because their dihedral
derivative is undefined; the water rings use 4-degree donor tilts.

Reproduce the retained water results in a new directory:

```sh
python examples/water_cluster.py --verify --output validation/my-water-check
python -m pytest -q
```

## Original Zn/O validation (0.1)

- **46 tests passed**, including finite-difference gradients, invariance,
  charge conservation, input validation, and reference comparisons.
- **27 retained reference cases passed**, including isolated atoms, dimers,
  all Zn/O angular environments in the file, an O4 torsional geometry,
  8/20/100/200-atom clusters, bond-stretch/cutoff scans, and relaxed structures.
- Both the ZnO dimer and 20-atom cluster converged under the full energy
  gradient. The 20-atom cluster took 100 L-BFGS iterations and reached a
  maximum force component of approximately 6.53e-5 kcal/mol/Å.

The original Zn/O report is [verified/summary.json](verified/summary.json).
Each case directory retains the reference input, exact parameter file, log,
atom dump, Python results, and executable metadata. Earlier exploratory runs
are retained separately and are not the current acceptance report.

## Performance

These are the original 0.1 baseline timings, before the additional property
outputs in 0.2. Measured medians of three repeats on macOS x86_64 Python 3.12
environment, with NumPy 2.5.3 and Autograd 1.9.1. Evaluation includes both the
full energy-gradient forces and the LAMMPS-convention forces.

| Atoms | Energy only | Full evaluation | Peak process RSS |
| --- | --- | --- | --- |
| 20 | 0.006 s | 0.036 s | 75.5 MiB |
| 100 | 0.094 s | 0.361 s | 130.7 MiB |
| 200 | 0.227 s | 0.935 s | 203.2 MiB |

See [benchmark.json](benchmark.json) for all timings and environment details.
Peak RSS includes the interpreter, numerical libraries, and relaxation, not
just array storage. Benchmark relaxation runs were capped at 20 iterations
and did not converge; they are timing measurements, separate from the
converged relaxation checks above. These results are geometry- and machine-
dependent, not a performance guarantee.

## Force comparison convention

Current reference agreement is measured using default fixed-charge `forces`.
Charge-response forces (`full_derivative=True`) pass independent energy
finite-difference tests with charges re-equilibrated. The distinction is
necessary because the reference QEq and
energy routines use inconsistent electrostatic conversion constants; see
the [force-convention guide](../docs/calculations.md#force-convention).
