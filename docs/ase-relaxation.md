# ASE and geometry relaxation

## ASE calculator

Install `.[ase]` to use {py:class}`xreac.ase.ReaxFFCalculator` with ASE objects
and optimizers. The adapter reads `atoms.cell` and `atoms.pbc`; ASE requires
periodicity to be enabled explicitly on the atom object.

```{testcode} ase
from ase import Atoms
from ase.units import kcal, mol
from xreac.ase import ReaxFFCalculator

atoms = Atoms("OH2", positions=[[0, 0, 0], [0.97, 0, 0], [-0.243, 0.94, 0]],
              cell=[4.0]*3, pbc=True)
atoms.calc = ReaxFFCalculator(ForceField.bundled("qeq_ff.water"))
energy_ev = atoms.get_potential_energy()
forces_ev = atoms.get_forces()
assert forces_ev.shape == (3, 3)
assert abs(energy_ev - atoms.calc.evaluation.energy*kcal/mol) < 1e-10
```

ASE exposes energy in **eV** and forces in **eV/Å**. Charges use elementary
charge and dipoles use e Å. `atoms.calc.evaluation` retains the core result in
kcal/mol units, including bond properties and `cell_repetitions`.

ASE neighbor lists are the default (`neighbor_backend="ase"`). They retain
periodic image shifts and work directly on small input cells without replication.
To select the retained native method:

```python
atoms.calc.set(neighbor_backend="replicated")
# Switch back to the default:
atoms.calc.set(neighbor_backend="ase")
```

The last result records `evaluation.neighbor_backend`. The replication limit
`max_expanded_atoms` applies only to `"replicated"`; see
[neighbor handling](periodic.md#neighbor-backends).

The adapter caches results, invalidating them when positions, cell, periodicity,
or calculator settings change. Initial charge arrays must be finite and sum to
zero; their values do not freeze the QEq solution. Position constraints are
handled by ASE. Stress and net charge are unsupported.

## Relaxation convenience API

```python
relaxed = calc.relax(symbols, positions, cell=[6.0]*3,
                     force_tolerance=1e-4, max_iterations=500)
print(relaxed.converged, relaxed.iterations, relaxed.message)

native = calc.relax(symbols, positions, cell=[6.0]*3, backend="native")
```

The default `backend="ase"` uses ASE FIRE and ASE neighbor lists.
`backend="native"` keeps the original NumPy FIRE and replication method
available for benchmarks and installations without ASE.
Both re-equilibrate charges at every geometry and use **fixed-charge forces**.
The convenience method has no `full_derivative` option. It changes positions
while keeping the cell fixed.

Convergence requires the largest absolute Cartesian force component to be at
most `force_tolerance`, in kcal/mol/Å. Reaching `max_iterations` returns
`converged=False`; it does not imply that the final geometry is optimized.
Always inspect the convergence flag.

FIRE uses damped fictitious dynamics rather than an energy line search. The
reported energy need not decrease every step because the fixed-charge forces
are not exactly its total derivative through QEq. Weak periodic image torques
can make orientation relaxation much slower than bond relaxation.

## Direct ASE optimizers

```python
from ase.optimize import FIRE

with FIRE(atoms, logfile="relax.log", trajectory="relax.traj") as optimizer:
    converged = optimizer.run(fmax=1e-5, steps=500)
```

Direct ASE optimizers use the largest atomic **force-vector norm** in eV/Å,
which differs from the core convenience method's stopping criterion. Use the
adapter directly for optimizer selection, constraints, and trajectory output.
The default fixed-charge setting is the relaxation convention. The explicit
adapter option `full_derivative=True` is intended for single-point studies.

From a checkout, `examples/ase_water.py --verify` saves a FIRE trajectory and
checks the final geometry with LAMMPS. The script also supports `--optimizer BFGS`.
Independent LAMMPS minimizations are not expected to follow identical paths.
