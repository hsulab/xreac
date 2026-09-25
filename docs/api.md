# Python API

## Calculator

```{eval-rst}
.. autoclass:: xreac.Calculator
   :members: evaluate, relax
```

The constructor accepts a loaded `ForceField` and keyword-only
`max_expanded_atoms=512`. The latter limits automatic small-cell replication;
it must be a positive integer. The default core calculation and native FIRE do not require ASE.

`evaluate` takes atom labels and a finite `(N, 3)` coordinate array. Options are
`total_charge=0`, `cell=None`, `pbc=None`, `full_derivative=False`, and
`neighbors=None`. Pass `neighbors=(i, j, S)` to use externally built directed
neighbor arrays; omit it to use native replication. `total_charge` is a finite real
scalar in e, including fractional values, and applies to the input cell. Booleans
and arrays are rejected. The return value is an `Evaluation`.

`relax` adds `force_tolerance=1e-4`, `max_iterations=500`, and `backend="ase"`.
It always uses fixed-charge forces and returns a `Relaxation`.

## Evaluation

```{eval-rst}
.. autoclass:: xreac.Evaluation
   :members: force_convention
```

| Attribute | Meaning |
| --- | --- |
| `energy` | Total kcal/mol, per input cell under PBC |
| `forces` | `(N, 3)` array in kcal/mol/Å |
| `charges` | `(N,)` equilibrated charges in e |
| `components` | Dictionary of 14 energy components in kcal/mol |
| `full_derivative` | Whether the selected forces include charge response |
| `force_convention` | `"fixed_charge"` or `"charge_response"` |
| `bond_orders` | `(N, N)` corrected bond orders, summed over periodic images |
| `total_bond_orders` | Per-atom row sums of `bond_orders` |
| `lone_pairs` | Per-atom lone-pair values |
| `bond_counts` | Per-atom counts of image bonds with order greater than 0.3 |
| `dipole` | `(3,)` vector in e Å on the supplied coordinate branch |
| `cell_repetitions` | Three native neighbor-search replication factors; `(1, 1, 1)` for supplied lists |
| `neighbor_backend` | `"replicated"`, `"provided"` for explicit arrays, or `"ase"` when the adapter supplies them |

## Relaxation result

```{eval-rst}
.. autoclass:: xreac.Relaxation
```

The fields are `positions`, `evaluation`, `converged`, `iterations`, and
`message`. The final evaluation uses fixed-charge forces even when convergence
was not reached. The cell is fixed throughout.

## ForceField

```{eval-rst}
.. autoclass:: xreac.ForceField
   :members: bundled, from_file, elements, vdw_type, validate_model
   :undoc-members:
```

Use the loading class methods rather than constructing parameter dictionaries
manually. Metadata fields include `path`, `checksum` (SHA256), and `citation`.
Parsed data are available in `general`, `atoms`, `pairs`, `angles`, `torsions`,
and `hydrogen_bonds`. See [parameter format](force-fields.md) for supported files.

## ASE adapter

```{eval-rst}
.. autoclass:: xreac.ase.ReaxFFCalculator
   :members: set
```

Implemented ASE properties are `energy`, `forces`, `charges`, and `dipole`.
Supported calculator parameters are `full_derivative=False`, `total_charge=0`,
`max_expanded_atoms=512`, `neighbor_backend="ase"`, and `neighbor_skin=0.3` (Å).
Net charge must be specified explicitly with `total_charge`; it is not inferred
from ASE initial charges. Initial charges may be absent. When supplied, they must
be finite and sum to `total_charge` within `1e-8` e; their distribution does not
constrain the fresh QEq solution.

The skin controls ASE topology reuse; zero rebuilds at every evaluation.
`neighbor_list_builds` counts builds over the adapter's lifetime.
The expansion limit applies only when `neighbor_backend="replicated"`.
Update parameters with `atoms.calc.set(...)` to invalidate
the cache. ASE standard calculator initialization options are accepted through
`**kwargs`. `atoms.calc.evaluation` retains the last core `Evaluation`.

## LAMMPS reference

```{eval-rst}
.. autofunction:: xreac.reference.evaluate_lammps
```

```{eval-rst}
.. autoclass:: xreac.reference.ReferenceResult
```

The result includes `energy`, `forces`, `charges`, `components`,
`total_bond_orders`, `lone_pairs`, `bond_counts`, `dipole`, `cell_repetitions`,
`version`, and `directory`. It does not contain the full bond-order matrix.
Its numerical units match the core calculator, including input-cell
normalization. The executable is needed only when calling the function.

Pass `supplied_charges=charges` to `evaluate_lammps` to use a finite `(N,)`
array of charges in e, including nonzero net charge. This disables LAMMPS QEq
and verifies energies and fixed-charge forces with the supplied values.
Small-cell replication tiles the charges, and results are normalized back to
the input cell. Omit the argument for the existing neutral QEq reference.
