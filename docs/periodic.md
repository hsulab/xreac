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

## Small cells

Small cells are replicated internally until every periodic face height exceeds
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
`Calculator(ff, max_expanded_atoms=...)` or the same ASE calculator option to
change it. This is a replication limit, not a limit on an already-large input
cell. Pair work and memory still grow with the expanded system.

## Current neighbor representation

The implementation uses dense pair arrays, not a spatially binned or Verlet
neighbor list. At every geometry it:

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
atom count. A direct image-aware sparse neighbor implementation could reduce
the replicated work; it is not the current backend.

## Reference convention

Small-cell verification uses equivalent **larger LAMMPS supercells**, normalized
to the input cell. Direct primitive-cell LAMMPS hydrogen-bond exclusions can
depend on original atom IDs. The reference harness expands small inputs by
default and records the mapping. Its `allow_small_cell=True` option is reserved
for diagnosing raw primitive-cell behavior.

Electrostatics follow the shielded, tapered finite-cutoff ReaxFF model. There is
no Ewald or PME summation. See [verification](validation.md) for retained results.
