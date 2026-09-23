# Verification with LAMMPS

The reference executable is `lmp_mpi`, pinned to **LAMMPS 22 Jul 2025, Update 4**
with REAXFF support. The harness invokes one rank directly, without `mpirun` or
LAMMPS Python bindings. It writes zero initial charges and runs a single-point
`run 0` with active QEq, a tolerance of `1e-12`, and up to 2000 iterations.
Coordinates do not move during this calculation.

## Run a reference

```python
from xreac.reference import evaluate_lammps

reference = evaluate_lammps(ff, symbols, positions, cell=[4.0]*3,
                            directory="reference-water")
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
mamba run -n catorch3 python -m pytest -q
mamba run -n catorch3 python -m pytest -q -m 'not reference'
mamba run -n catorch3 python -m pytest -q -m reference
```

The 0.6 implementation passed **163 tests** in the recorded `catorch3` run.
This is a historical validation count, not a dynamically executed docs build.
LAMMPS tests are not run while building documentation.

## Small-cell verification

Nine retained cases cover 3.12–9 Å water cells, a water dimer, tilted and partially
periodic cells, Zn/O, and a one-atom zinc chain with bonds to its own images.
The maximum fixed-charge force discrepancy is **4.2e-12 kcal/mol/Å** against
normalized LAMMPS supercells; the maximum energy discrepancy per input atom is
**4.3e-12 kcal/mol**. Download the
{download}`small-cell summary <../validation/small-cell-support/summary.json>`.

The earlier primitive-cell audit found that the pinned LAMMPS hydrogen-bond
implementation excludes donor/acceptor images sharing an original atom ID.
At a 4 Å water cell, that omits a -0.370052 kcal/mol/cell contribution present in
the equivalent supercell. The Python small-cell implementation includes it.
QEq charges agreed in the audited single-rank cases; this does not supersede
the [LAMMPS QEq small-cell restriction](https://docs.lammps.org/fix_qeq_reaxff.html#restrictions).
Download the {download}`original diagnostic <../validation/small-cells/summary.json>`.

## Other retained results

- {download}`Periodic water, including 192 atoms <../validation/periodic-water/summary.json>`.
- {download}`Water clusters with fixed-charge forces <../validation/water-fixed-charge/summary.json>`.
- {download}`Illustrated water report (PDF) <../validation/water-fixed-charge/report.pdf>`.
- {download}`ASE FIRE relaxation <../validation/relaxation-ase-catorch3/summary.json>`.
- {download}`Native FIRE relaxation <../validation/relaxation-native-catorch3/summary.json>`.

These are deterministic verification geometries, not equilibrated liquid
snapshots or predictions of physical stability. Agreement establishes
implementation consistency, not agreement with quantum chemistry.

## QEq audit

The retained {download}`QEq audit <../validation/qeq-audit/summary.json>`
compares LAMMPS forces with finite differences of its own energies, solving QEq
again at each displaced geometry. Controls repeat QEq at fixed coordinates
and change the initial neutral charge guess. The charge-response difference
persists and follows the conversion-constant mismatch described in
[force conventions](calculations.md#force-convention).
