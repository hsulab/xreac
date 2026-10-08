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

[Validation](../validation/README.md) retains only extended XYZ input files
with `.xyz` extensions and brief READMEs for water, Zn/O, C/H/O, CuO,
Ag/ZnO, anatase and the pending Pd/oxide examples. Parameter provenance stays
in `data/`; calculation outputs are not retained in these system folders.
Historical detailed reports remain in Git history.

```sh
python scripts/validate.py --verify --output validation/runs/check
python scripts/validate.py --system water --include-charged --verify
python scripts/validate.py --system cho --include-chocl --verify
python scripts/validate.py --system cuo --verify
python scripts/validate.py --system anatase --verify
```

Use fresh output directories. Generated system directories contain
`structures.json`, `summary.json`, `results.json.gz` and `reference.tar.gz`
for inspecting energies, forces, charges, dipoles and bond properties.
Charged inputs carry their net charge in XYZ metadata; supply that charge
explicitly when constructing a calculator. Shared fixtures are defined in
`scripts/validate.py`. `--include-bulk` enables optional larger structures.

To generate the water report from a fresh verified run:

```sh
python scripts/water_report.py --input validation/runs/check/water \
  --output validation/runs/check/water-report.pdf
```

ASE/native relaxation, finite differences, cutoff perturbations, symmetries,
cache checks and neighbor rebuilding reuse these fixtures. They establish
implementation consistency, not equilibrated liquid structure or agreement
with quantum chemistry. Performance runs remain optional.

## Small-cell verification

The 4 Angstrom water cell tests repeated periodic images. Comparisons use
normalized larger LAMMPS supercells: the primitive-cell implementation excludes
hydrogen-bond acceptors sharing a donor's original atom ID, omitting
-0.370052 kcal/mol/cell present in the equivalent supercell. QEq charges agree.
See the [LAMMPS QEq restriction](https://docs.lammps.org/fix_qeq_reaxff.html#restrictions).

## QEq and charged-system checks

The isolated monomer audit compares fixed-charge forces with finite
differences after resolving QEq at each displacement. The difference follows
the conversion-constant mismatch in [force conventions](calculations.md#force-convention).

Charged fixtures include hydroxide, hydronium and periodic charged water.
An independent constrained QEq solve checks charge conservation, chemical
potential equality and positive reduced curvature. LAMMPS receives supplied
charges and disables QEq (`checkqeq no`) for nonneutral comparisons; matching
those supplied charges does not independently validate QEq. Neutral fixtures
use fresh LAMMPS QEq. Hur2021 adds three neutral molecules and two charged
SN2 geometries. Numerical agreement does not validate reaction barriers.
