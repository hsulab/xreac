# Validation record

Verified locally on September 22, 2026 with LAMMPS 22 Jul 2025, Update 4,
invoked through `/opt/homebrew/bin/lmp_mpi`. The force-field checksums and
original citations are recorded in [NOTICE](../NOTICE).

## ASE and catorch3 (0.4)

The local virtual environment was removed. Development and validation now use
`mamba run -n catorch3`, with Python 3.10.13, NumPy 2.0.2, Autograd 1.9.1,
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

See the retained [ASE results](relaxation-ase-catorch3/summary.json) and
[native results](relaxation-native-catorch3/summary.json). Direct ASE examples
also pass: [FIRE](ase-water-fire/summary.json) converges in 61 steps and
[BFGS](ase-water-bfgs/summary.json) in 6 steps, with a 1e-5 eV/Å vector-norm
tolerance. Their directories contain trajectories, optimizer logs, structures,
and raw LAMMPS verification inputs/outputs.

The [catorch3 calculation benchmark](benchmark-catorch3.json) records fixed-charge
evaluations only. The [native FIRE benchmark](benchmark-native-fire-catorch3.json)
adds 20-step relaxation timings using our retained implementation. These runs
use a different interpreter/NumPy environment from the archived benchmarks;
do not attribute timing differences solely to code changes. Capped relaxation
timings do not imply convergence.

```sh
mamba run -n catorch3 python -m pytest -q
mamba run -n catorch3 python scripts/validate_relaxation.py --backend ase --output validation/runs/ase-check
mamba run -n catorch3 python scripts/validate_relaxation.py --backend native --output validation/runs/native-check
mamba run -n catorch3 python examples/ase_water.py --verify
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

The dedicated water example uses LAMMPS's `qeq_ff.water` (Achtyl et al.,
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
python -m pip install '.[report]'
python scripts/water_report.py
# Archived two-force report:
python scripts/water_report.py --input validation/water --output validation/water/report.pdf
```

Matplotlib is optional and is not required by the calculator. To report a new
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
the [main README](../README.md#force-derivative-option).
