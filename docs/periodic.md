# Periodic cells and neighbor handling

## Cell conventions

The core API accepts three cell lengths or a `(3, 3)` matrix whose rows are
lattice vectors, in Å. Cells must be finite, nonsingular, and right-handed.
Supplying `cell` enables full periodicity unless `pbc` is explicitly specified.
`pbc` accepts one boolean or three flags selecting periodic lattice directions.

```python
result = calc.evaluate(symbols, positions, cell=[12.48]*3)
slab = calc.evaluate(symbols, positions, cell=[12.48, 12.48, 30.0],
                     pbc=[True, True, False])
isolated = calc.evaluate(symbols, positions, cell=[12.48]*3, pbc=False)
```

Positions may be wrapped or unwrapped. Adding lattice vectors to individual
atoms preserves energies, forces, charges, and bond properties. Dipoles follow
the supplied coordinate branch. Cells remain fixed during relaxation; stress
and variable-cell optimization are not implemented.

## Neighbor backends

The two backends evaluate the same image-resolved interactions:

| Backend | Default for | Neighbor construction |
| --- | --- | --- |
| Native replication | Core `Calculator.evaluate()` without `neighbors` and native FIRE | Dense image search, with tied copies for small cells |
| Supplied `(i, j, S)` | `ReaxFFCalculator` and ASE FIRE build these with ASE | Directed neighbor arrays on the input cell, retaining integer image shifts |

```{testcode} supplied-neighbors
from ase import Atoms
from ase.neighborlist import neighbor_list

ff = ForceField.bundled("qeq_ff.water")
calc = Calculator(ff)
atoms = Atoms("OH2", positions=[[0, 0, 0], [0.97, 0, 0], [-0.243, 0.94, 0]],
              cell=[4.0]*3, pbc=True)
i, j, S = neighbor_list("ijS", atoms, ff.general[12])
direct = calc.evaluate(atoms.get_chemical_symbols(), atoms.positions,
                       cell=atoms.cell.array, pbc=atoms.pbc, neighbors=(i, j, S))
native = calc.evaluate(atoms.get_chemical_symbols(), atoms.positions,
                       cell=atoms.cell.array, pbc=atoms.pbc)
assert direct.cell_repetitions == (1, 1, 1)
assert direct.neighbor_backend == "provided"
assert abs(direct.energy - native.energy) < 1e-9
assert np.max(abs(direct.forces - native.forces)) < 1e-9
```

The ASE adapter needs `.[ase]`. It calls
[ASE's neighbor list](https://docs.ase-lib.org/ase/neighborlist.html)
with the force field's nonbonded cutoff **before** calling the core evaluator,
then passes `neighbors=(i, j, S)` explicitly. The core has no ASE dependency and
performs no neighbor search when arrays are supplied. Other builders can supply
the same tuple. Every edge stores two input atom indices
and a lattice shift; its vector is `positions[j] - positions[i] + shift @ cell`.
Both edge directions and all images within the cutoff are retained. Only the
zero-shift self interaction is excluded. Bond-order, angle, torsion, and
hydrogen-bond screening then applies the same cutoffs as the native backend.
ASE's covalent-radius defaults are not used.

Angles and torsions track image identities along each bond chain. A hydrogen-bond
donor and acceptor can share an input atom index if they occupy different images.
Charges are solved on the input atoms; self-image couplings enter the QEq diagonal.
Bond counts apply their threshold per image before reduction to the input atoms.

The ASE adapter rebuilds the neighbor list for each new evaluation, including
changes to positions, cell, or PBC. During autodiff, the selected edges and lattice shifts
stay fixed while their vectors and distances remain differentiable. There is
no skin or list reuse between changed geometries; ASE caches unchanged results.
When calling the core with your own arrays, **you** must keep the list complete
for the current coordinates, cell, and PBC. Both edge directions are required;
half lists, duplicates, out-of-range indices, and zero-shift self edges are
rejected. Completeness cannot be checked without rebuilding the list, so it is
the caller's responsibility. Extra neighbors from a skin are supported and are
screened by the physical cutoffs. The core copies the arrays and does not mutate them.

Supplied arrays always give `cell_repetitions=(1, 1, 1)` and do not use
`max_expanded_atoms`. Pair parameters and distances are stored per edge. The QEq
matrix and returned bond-order matrix remain dense in the input atom count.

## Small cells with native replication

Without supplied arrays (or with `neighbor_backend="replicated"` on the ASE adapter),
small cells are replicated internally until every periodic face height exceeds
both the nonbonded cutoff and twice the bond cutoff—10 Å for the bundled files.
Nonperiodic directions are not repeated. For tilted cells, face heights rather
than lattice-vector lengths determine replication.

```{testcode} small-cell
ff = ForceField.bundled("qeq_ff.water")
calc = Calculator(ff)
symbols = ["O", "H", "H"]
positions = [[0, 0, 0], [0.97, 0, 0], [-0.243, 0.94, 0]]
result = calc.evaluate(symbols, positions, cell=[4.0]*3)
assert result.cell_repetitions == (3, 3, 3)
assert result.forces.shape == (3, 3)
assert result.components["hydrogen_bond"] < -0.3
```

This calculation has three input atoms and 81 internal atoms. Translated
copies have distinct interaction identities but share the input coordinates
and charges. QEq solves only for the input atoms. The energy is divided by the
number of copies before differentiation, folding all image forces back onto
the input coordinates. Returned energies and properties refer to the input cell.

The default limit is **512 internal atoms** for automatic replication. Exceeding
it raises an error before expanded pair arrays are allocated. Set
`Calculator(ff, max_expanded_atoms=...)` or the same option with
`ReaxFFCalculator(ff, neighbor_backend="replicated")` to change it.
This is a replication limit, not a limit on an already-large input
cell. Pair work and memory still grow with the expanded system.

## Native neighbor representation

The retained native backend uses dense pair arrays. At every geometry it:

1. Computes displacements between all internal atom pairs.
2. Examines shifts `-1`, `0`, and `+1` in each periodic direction, up to 27 image
   offsets after centering fractional displacements.
3. Uses the nearest bonded image with bond-order screening for bonded terms.
4. Sums all images within the nonbonded and hydrogen-bond cutoffs.

The internal cell-height condition guarantees a unique bonded image per pair
and sufficient image coverage. Distinct copies represent multiple neighbors
of the same primitive atom. Arrays are rebuilt for each evaluated geometry;
there is no neighbor skin or reuse between changed geometries. ASE caches a
result when its input has not changed.

Pair storage scales quadratically in the internal atom count, multiplied by
the number of image offsets. QEq's dense solve scales cubically in the input
atom count. The ASE backend avoids these expanded pair arrays.

## Reference convention

Small-cell verification uses equivalent **larger LAMMPS supercells**, normalized
to the input cell. Direct primitive-cell LAMMPS hydrogen-bond exclusions can
depend on original atom IDs. The reference harness expands small inputs by
default and records the mapping. Its `allow_small_cell=True` option is reserved
for diagnosing raw primitive-cell behavior.

Electrostatics follow the shielded, tapered finite-cutoff ReaxFF model. There is
no Ewald or PME summation. See [verification](validation.md) for retained results.
