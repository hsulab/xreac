# Validation record

Verified locally on September 22, 2026 with LAMMPS 22 Jul 2025, Update 4,
invoked through `/opt/homebrew/bin/lmp_mpi`. The force-field checksums and
original citations are recorded in [NOTICE](../NOTICE).

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

Reference agreement is measured using `lammps_forces`. The primary `forces`
field passes independent energy finite-difference tests with charges
re-equilibrated. The distinction is necessary because the reference QEq and
energy routines use inconsistent electrostatic conversion constants; see
the [main README](../README.md#two-force-conventions).
