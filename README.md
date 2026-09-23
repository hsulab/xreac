# xreac

A small NumPy + Autograd ReaxFF calculator for **neutral, isolated molecules
and clusters**, targeting roughly 1–200 atoms. Elements and interactions come
from a standard ReaxFF parameter file. It computes energies, charges, forces,
bond properties, and dipoles, and relaxes geometries using fixed-charge forces.
No LAMMPS installation is required for calculations; `lmp_mpi` is used for
independent verification.

## Install

```sh
mamba run -n catorch3 python -m pip install -e '.[ase,test]'
```

Development uses the existing `catorch3` mamba environment; there is no local
virtual environment. Python 3.10 or newer is required. Core calculations and
the native FIRE optimizer need only NumPy and HIPS Autograd. ASE is an optional
dependency for the adapter and default relaxation backend; install `.[ase]`
or `.[relax]` to enable it. A minimal installation is `python -m pip install .`.

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

# Optional: include the response of the equilibrated charges in the derivative.
full = calc.evaluate(symbols, positions, full_derivative=True)
print(full.forces)

# Default relaxation reuses ASE FIRE with fixed-charge forces.
relaxed = calc.relax(symbols, positions, force_tolerance=1e-4, max_iterations=500)
print(relaxed.converged, relaxed.message)
print(relaxed.positions)

# Keep the original NumPy FIRE implementation available for benchmarks.
native = calc.relax(symbols, positions, backend="native")
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

## ASE calculator and optimizers

```python
from ase import Atoms
from ase.optimize import FIRE  # BFGS and other ASE optimizers can also use the adapter
from xreac import ForceField
from xreac.ase import ReaxFFCalculator

atoms = Atoms("OH2", positions=[[0, 0, 0], [0.97, 0, 0], [-0.243, 0.94, 0]])
atoms.calc = ReaxFFCalculator(ForceField.bundled("qeq_ff.water"))
print(atoms.get_potential_energy())  # eV
print(atoms.get_forces())           # eV/Angstrom, fixed-charge convention
print(atoms.get_charges())          # elementary charge
print(atoms.get_dipole_moment())    # e Angstrom
with FIRE(atoms, logfile="relax.log", trajectory="relax.traj") as optimizer:
    converged = optimizer.run(fmax=1e-5, steps=500)
```

The adapter caches results and recalculates after geometry or parameter changes.
It supports ASE position constraints, neutral initial charge guesses, and
nonperiodic display boxes. Periodic systems, net charges, and stress are not
supported. ASE chemical symbols must match parameter-file atom labels; use
the core API and `backend="native"` for nonstandard labels.
`atoms.calc.evaluation` retains the last core result, including
bond properties, in **kcal/mol units**, while ASE energy/force methods use
**eV units**. Single-point charge-response forces can be selected explicitly
with `ReaxFFCalculator(ff, full_derivative=True)` or
`atoms.calc.set(full_derivative=True)`; relaxation examples use the default False.

Direct ASE optimizers stop on the largest atomic **force-vector norm** in
eV/Angstrom. The `Calculator.relax()` convenience method retains the existing
maximum **Cartesian component** criterion in kcal/mol/Angstrom for both backends.
Use the adapter directly for optimizer selection, constraints, and trajectories.

```sh
mamba run -n catorch3 python examples/ase_water.py --verify
mamba run -n catorch3 python examples/ase_water.py --optimizer BFGS --verify
```

## Water-cluster example

```sh
mamba run -n catorch3 python examples/water_cluster.py
mamba run -n catorch3 python examples/water_cluster.py --verify --output validation/my-water-check
```

The example evaluates a monomer, dimer, distorted dimer, trimer, and hexamer
with the dedicated `qeq_ff.water` file distributed in LAMMPS's water example.
The structures are reproducible test geometries, not optimized clusters.
`--verify` checks all energy components, every force component, charges,
dipole vectors, per-atom bond-order sums, lone pairs, and bond counts against
`lmp_mpi`. It writes XYZ structures, Python results, reference inputs/outputs,
and a comparison summary. The output directory must not already exist.
Use `--ffield path/to/file` to select another compatible parameter set.

## Force derivative option

`evaluate(..., full_derivative=False)` is the default. It returns
**fixed-charge forces**, matching the LAMMPS convention: QEq is solved at the
input geometry, then those equilibrated charges are held constant while
differentiating the energy. Charges are recalculated for every new geometry;
"fixed-charge" describes only the derivative.

`evaluate(..., full_derivative=True)` returns **charge-response forces**, the
full negative derivative of the reported energy through the QEq solve.

Both modes return the selected array in `result.forces`, with
`result.full_derivative` and `result.force_convention` (`"fixed_charge"` or
`"charge_response"`) recording the choice. Only the selected derivative is
computed. Energies, charges, and other properties are identical between modes.
Both use automatic differentiation, not finite differences.

**API change in 0.3:** the separate `result.lammps_forces` field was removed.
Use default `result.forces` in its place. Code that previously used `forces`
for the full energy derivative must now pass `full_derivative=True`.

The two derivatives differ because LAMMPS uses `14.4` in the QEq coupling, `332.06371` in
the Coulomb energy, and `23.02` in the self-energy conversion;
`14.4 * 23.02 != 332.06371`. Thus its QEq solution is not precisely stationary
for its reported energy. For a ZnO dimer at 1.9 Å, the force difference is
approximately 0.038 kcal/mol/Å. A geometry stationary under charge-response
forces can have nonzero fixed-charge forces.

`relax()` always uses fixed-charge forces and returns an evaluation with
`full_derivative=False`. At each geometry it equilibrates charges, computes
fixed-charge forces, and updates positions with ASE FIRE (`backend="ase"`,
the default). `backend="native"` selects our retained FIRE implementation
without importing ASE. Both use damped fictitious dynamics, do not
differentiate through QEq, and do not perform an
energy line search. Convergence requires the largest absolute Cartesian force
component to be below `force_tolerance`; reaching `max_iterations` does not
imply convergence. Reported energies need not decrease at every step.

This uses the same force convention as LAMMPS, but does not promise the same
optimization trajectory or local minimum. The FIRE algorithm is described in
[Bitzek et al., Phys. Rev. Lett. 97, 170201 (2006)](https://doi.org/10.1103/PhysRevLett.97.170201).
See the [ASE relaxation results](validation/relaxation-ase-catorch3/summary.json)
and [native FIRE results](validation/relaxation-native-catorch3/summary.json).

The [single-point QEq audit](validation/qeq-audit/README.md) verifies that
`run 0` equilibrates charges and compares LAMMPS forces against finite
differences of LAMMPS's own energies, with fresh QEq at every displacement.
The force discrepancy persists after repeating QEq and changing the initial
charges; it is not caused by missing charge equilibration.

## Verify with lmp_mpi

```sh
mamba run -n catorch3 python -m pytest -q
```

Reference tests require `lmp_mpi` on PATH, or an explicit override:

```sh
XREAC_LAMMPS=/opt/homebrew/bin/lmp_mpi mamba run -n catorch3 python -m pytest -q -m reference
mamba run -n catorch3 python -m pytest -q -m 'not reference'  # checks without LAMMPS
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
mamba run -n catorch3 python scripts/validate.py
mamba run -n catorch3 python scripts/benchmark.py
# Optional timing of our retained FIRE optimizer:
mamba run -n catorch3 python scripts/benchmark.py --relax-iterations 20 --output validation/benchmark-native-fire-catorch3.json
```

Benchmarks compute only fixed-charge forces (`full_derivative=False`). Default
runs exclude relaxation; `--relax-iterations N` adds timings for **native FIRE**,
also with fixed-charge forces. Default results go to
`validation/benchmark-catorch3.json`. Environment metadata is recorded, and
previous benchmarks remain historical records. Iteration-capped relaxation
timings are not claims of convergence.

Acceptance targets are 1e-5 kcal/mol/atom for total/component energy,
1e-6 e for charges, and 1e-4 kcal/mol/Å for default fixed-charge `forces`.
Independent finite-difference tests check charge-response forces with QEq
re-solved on each displacement, and fixed-charge forces with charges held fixed.

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
