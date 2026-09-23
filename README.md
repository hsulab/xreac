# xreac

A small NumPy + Autograd ReaxFF calculator for **neutral, isolated molecules
and clusters**, targeting roughly 1–200 atoms. Elements and interactions come
from a standard ReaxFF parameter file. It computes energies, charges, forces,
bond properties, and dipoles, and optionally relaxes geometries using SciPy.
No LAMMPS installation is required for calculations; `lmp_mpi` is used for
independent verification.

## Install

```sh
python -m pip install .
python -m pip install '.[relax]'  # optional geometry optimization
```

Python 3.10 or newer is required. The core dependencies are NumPy and HIPS
Autograd; SciPy is only imported when relaxation is requested.

Source parameter files live in `data/` at the repository root. Distributions
bundle them as package data so `ForceField.bundled(name)` also works outside
a source checkout. Bundled examples are `qeq_ff.water` (Achtyl QEq water),
`ffield.reax.cho` (Chenoweth C/H/O), and `ffield.reax.ZnOH` (Raymand Zn/O/H).
`ForceField.zno()` remains a backward-compatible shortcut for the last file.

## Calculate

```python
from xreac import Calculator, ForceField

ff = ForceField.bundled("qeq_ff.water")
# Or use your own supported parameter file:
# ff = ForceField.from_file("my_force_field.ff")
print(ff.elements)  # ('H', 'O', 'X'); labels are defined by the parameter file
calc = Calculator(ff)
symbols = ["O", "H", "H"]
positions = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0], [-0.243, 0.94, 0.0]]
result = calc.evaluate(symbols, positions)
print(result.energy)
print(result.forces)
print(result.charges)
print(result.components)
print(result.dipole)  # e Angstrom
print(result.bond_orders)  # symmetric (N, N) corrected bond-order matrix

relaxed = calc.relax(symbols, positions, force_tolerance=1e-4, max_iterations=500)
print(relaxed.converged, relaxed.message)
print(relaxed.positions)
```

Coordinates and forces have shape `(N, 3)`; charges have shape `(N,)`.
Units are Å, kcal/mol, kcal/mol/Å, and elementary charge. Calculations use
float64. Energy components correspond to the 14 LAMMPS `compute pair reaxff`
entries, with zero entries retained for inactive terms. Bond energies,
coordination and lone-pair energies, angles, penalties and conjugation,
torsions, hydrogen bonding, vdW, Coulomb energy, and QEq self-energy are included.
The standard C2 correction and terminal triple-bond stabilization are included.

Additional results are `total_bond_orders` (per-atom row sums), `lone_pairs`,
and `bond_counts` (per-atom counts of corrected bond orders **greater than 0.3**,
matching the default LAMMPS bond-graph cutoff). They have shape `(N,)`.
`dipole` is a three-component vector in e Å; for a neutral system it is
independent of the coordinate origin. Dummy labels such as `X` remain in
parameter metadata; they do not denote an additional chemical element.

## Water-cluster example

```sh
python examples/water_cluster.py
python examples/water_cluster.py --verify --output validation/my-water-check
```

The example evaluates a monomer, dimer, distorted dimer, trimer, and hexamer
with the dedicated `qeq_ff.water` file distributed in LAMMPS's water example.
The structures are reproducible test geometries, not optimized clusters.
`--verify` checks all energy components, every force component, charges,
dipole vectors, per-atom bond-order sums, lone pairs, and bond counts against
`lmp_mpi`. It writes XYZ structures, Python results, reference inputs/outputs,
and a comparison summary. The output directory must not already exist.
Use `--ffield path/to/file` to select another compatible parameter set.

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
ref = evaluate_lammps(ff, symbols, positions, directory="reference-water")
```

Use a new or empty directory. It retains the exact parameter file, structure,
input script, log, atom dump, energies, stdout/stderr, version, and SHA256.
Reference tests use temporary directories. The original Zn/O regression and
benchmark scripts remain available:

```sh
python scripts/validate.py
python scripts/benchmark.py
```

Acceptance targets are 1e-5 kcal/mol/atom for total/component energy,
1e-6 e for charges, and 1e-4 kcal/mol/Å for `lammps_forces`. Independent
finite-difference tests check `forces` with QEq re-solved on each displacement.

Property comparisons additionally require dipoles within 1e-6 e Å,
bond-order sums and lone pairs within 1e-8, and identical bond counts.
The results and retained water and Zn/O reference runs are recorded in
[validation/README.md](validation/README.md).

## Supported parameter format

The reader supports the standard **39-global-parameter, four-line-atom**
LAMMPS ReaxFF format. It preserves atom-type labels and reads explicit bonds,
off-diagonal overrides, repeated angle entries, hydrogen-bond parameters,
and explicit or terminal `0-i-j-0` wildcard torsions. Explicit torsions take
precedence over wildcard defaults, including when they appear later in the
file. Nonbonded mixing covers every atom-type pair. Fortran `D` exponents and
an omitted final hydrogen-bond block are accepted.

Calculation support is determined by parameters, not by an element whitelist.
Standard shielded vdW, inner-wall vdW, and their combination are implemented.
Missing atom labels and unsupported formats/settings raise errors. The
dedicated water, C/H/O, and Zn/O regressions establish the tested coverage;
a successfully parsed arbitrary parameter set still needs its own validation.

## Limits

- Reference agreement establishes implementation consistency, not accuracy
  against quantum chemistry or suitability for a particular system. The bundled
  ZnOH file is the **2010** parameterization, not the earlier 2008 ZnO set.
- No periodic cells, stress, MD, net charge, external fields, parameter fitting,
  lgvdW/five-line-atom extensions, nonzero lower taper radius, or alternative
  charge models such as ACKS2 are supported. The standard LAMMPS defaults of
  5 Å for bond candidates and 7.5 Å for hydrogen bonds are used (also limited
  by the force field's nonbonded cutoff).
- The reference model has branch thresholds in bond orders and interaction
  selection. Derivatives are local to the active branches, not derivatives of
  discrete list membership. Exactly collinear **active torsions** raise an
  explicit error because their dihedral derivative is undefined. Angle
  denominators are guarded; singular geometries should not be used to assess
  derivative agreement. Even weak intermolecular bonds can activate torsions.
- Dense charge equilibration takes cubic time in atom count and quadratic
  memory. Pair arrays are dense; angle/torsion lists contain local interactions.
  Highly compressed or densely connected systems may be expensive or unstable.
- Global QEq permits charge transfer between separated fragments. This is a
  property of the reference model, not a guarantee of physical dissociation.

## Attribution

GPL-2.0-or-later. Equations and conventions are adapted from the LAMMPS/PuReMD
implementation. The unmodified parameter files retain their original citations.
See [NOTICE](NOTICE) and [LICENSE](LICENSE) for attribution and provenance.
