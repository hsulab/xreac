# xreac

A small NumPy + Autograd ReaxFF calculator for **neutral, isolated Zn/O
clusters**, targeting roughly 1–200 atoms. It computes energies, charges,
forces, and optionally relaxes geometries using SciPy. No LAMMPS installation
is required for calculations; `lmp_mpi` is used for independent verification.

## Install

```sh
python -m pip install .
python -m pip install '.[relax]'  # optional geometry optimization
```

Python 3.10 or newer is required. The core dependencies are NumPy and HIPS
Autograd; SciPy is only imported when relaxation is requested.

The source parameter file lives in `data/ffield.reax.ZnOH` at the repository
root. Distributions bundle it as package data so `ForceField.zno()` also works
outside a source checkout.

## Calculate

```python
from xreac import Calculator, ForceField

ff = ForceField.zno()  # bundled Raymand 2010 ZnOH parameters
# In a source checkout: ff = ForceField.from_file("data/ffield.reax.ZnOH")
calc = Calculator(ff)
symbols = ["Zn", "O"]
positions = [[0.0, 0.0, 0.0], [2.1, 0.0, 0.0]]
result = calc.evaluate(symbols, positions)
print(result.energy)
print(result.forces)
print(result.charges)
print(result.components)

relaxed = calc.relax(symbols, positions, force_tolerance=1e-4, max_iterations=500)
print(relaxed.converged, relaxed.message)
print(relaxed.positions)
```

Coordinates and forces have shape `(N, 3)`; charges have shape `(N,)`.
Units are Å, kcal/mol, kcal/mol/Å, and elementary charge. Calculations use
float64. Energy components correspond to the 14 LAMMPS `compute pair reaxff`
entries, with zero entries retained for inactive terms. Bond energies,
coordination and lone-pair energies, angles, angle penalties, O-only torsions
and conjugation, shielded vdW, Coulomb energy, and QEq self-energy are included.

## Two force conventions

`result.forces` is **minus the derivative of the reported energy**, including
the response of QEq charges. This is the force used for geometry relaxation
and finite-difference checks.

`result.lammps_forces` differentiates with converged charges held fixed,
matching the convention in the reference LAMMPS implementation. These are
not exactly the same: LAMMPS uses `14.4` in the QEq coupling, `332.06371` in
the Coulomb energy, and `23.02` in the self-energy conversion;
`14.4 * 23.02 != 332.06371`. Thus its QEq solution is not precisely stationary
for its reported energy. For a ZnO dimer at 1.9 Å, the force difference is
approximately 0.038 kcal/mol/Å. This is recorded explicitly rather than
absorbed into test tolerances. A geometry stationary under `forces` can have
nonzero `lammps_forces`.

## Verify with lmp_mpi

```sh
python -m pip install '.[test]'
python -m pytest -q
```

Reference tests require `lmp_mpi` on PATH, or an explicit override:

```sh
XREAC_LAMMPS=/opt/homebrew/bin/lmp_mpi python -m pytest -q -m reference
python -m pytest -q -m 'not reference'  # checks without LAMMPS
```

The reference harness invokes one MPI rank directly, without `mpirun` or
LAMMPS's Python bindings. Missing executables, failed calculations, and QEq
nonconvergence cause failures rather than skipped validation. The initial
reference is LAMMPS **22 Jul 2025, Update 4** with the REAXFF package.
The harness checks this version; using a different reference requires an
explicit `expected_version` argument to `evaluate_lammps` and revalidation.

```python
from xreac.reference import evaluate_lammps
ref = evaluate_lammps(ff, symbols, positions, directory="reference-zno")
```

Use a new or empty directory. It retains the exact parameter file, structure,
input script, log, atom dump, energies, stdout/stderr, version, and SHA256.
Reference tests use temporary directories. The validation script retains a
complete reproducible suite and a JSON summary:

```sh
python scripts/validate.py
python scripts/benchmark.py
```

Acceptance targets are 1e-5 kcal/mol/atom for total/component energy,
1e-6 e for charges, and 1e-4 kcal/mol/Å for `lammps_forces`. Independent
finite-difference tests check `forces` with QEq re-solved on each displacement.

The initial results (46 passing tests, 27 saved reference cases, and timings
through 200 atoms) are recorded in [validation/README.md](validation/README.md).

## Limits

- The bundled **2010 ZnOH** file is the initial reference model, not the 2008
  ZnO set. Reference agreement establishes implementation consistency, not
  accuracy against quantum chemistry or suitability for a particular cluster.
  The reader accepts compatible standard parameter files, but other parameter
  sets have not been validated and do not inherit this model's compatibility claim.
- Only Zn/O calculations are enabled; H in the file is not supported. No
  periodic cells, stress, MD, net charge, external fields, parameter fitting,
  inner-wall vdW, lgvdW, or wildcard torsion variants are supported.
- The reference model has branch thresholds in bond orders and interaction
  selection. Derivatives are local to the active branches, not derivatives of
  discrete list membership. Exact collinear angles/torsions are singular;
  geometric denominators are guarded, and exact singular geometries should
  not be used to assess derivative agreement.
- Dense charge equilibration takes cubic time in atom count and quadratic
  memory. Pair arrays are dense; angle/torsion lists contain local interactions.
  Highly compressed or densely connected systems may be expensive or unstable.
- Global QEq permits charge transfer between separated fragments. This is a
  property of the reference model, not a guarantee of physical dissociation.

## Attribution

GPL-2.0-or-later. Equations and conventions are adapted from the LAMMPS/PuReMD
implementation. The unmodified parameter file retains its original citation.
See [NOTICE](NOTICE) and [LICENSE](LICENSE) for attribution and provenance.
