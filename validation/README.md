# Validation record

Verified locally on September 22, 2026 with LAMMPS 22 Jul 2025, Update 4,
invoked through `/opt/homebrew/bin/lmp_mpi`. The calculator uses the exact
bundled `ffield.reax.ZnOH` checksum recorded in [NOTICE](../NOTICE).

- **46 tests passed**, including finite-difference gradients, invariance,
  charge conservation, input validation, and reference comparisons.
- **27 retained reference cases passed**, including isolated atoms, dimers,
  all Zn/O angular environments in the file, an O4 torsional geometry,
  8/20/100/200-atom clusters, bond-stretch/cutoff scans, and relaxed structures.
- Both the ZnO dimer and 20-atom cluster converged under the full energy
  gradient. The 20-atom cluster took 100 L-BFGS iterations and reached a
  maximum force component of approximately 6.53e-5 kcal/mol/Å.

The complete current report is [verified/summary.json](verified/summary.json).
Each case directory retains the reference input, exact parameter file, log,
atom dump, Python results, and executable metadata. Earlier exploratory runs
are retained separately and are not the current acceptance report.

## Performance

Measured medians of three repeats on the installed macOS x86_64 Python 3.12
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
