# Energies, charges, and forces

## Charge equilibration

Every evaluated geometry receives a fresh, neutral QEq solution. A dense linear
system equalizes the charge chemical potentials while a Lagrange multiplier
enforces zero net charge. Previous or initial atomic charges are not frozen
across geometries.

For small cells, translated copies share the input atom's charge. QEq includes
coupling to nonzero images of the same atom and solves an `(N+1) × (N+1)` system
for the `N` input atoms. Global QEq can transfer charge between separated
fragments; it is not a guarantee of physically correct dissociation.

## Force convention

`evaluate(..., full_derivative=False)` is the default. These **fixed-charge
forces** hold the freshly equilibrated charges constant while differentiating
the energy. This follows the LAMMPS force convention. Charges are still
re-equilibrated at every new geometry; “fixed-charge” describes the derivative.

`evaluate(..., full_derivative=True)` selects **charge-response forces**:
the full negative derivative of the reported energy through QEq. Both choices
use automatic differentiation, not finite differences. Only the selected
derivative is computed and returned in `result.forces`.

```python
fixed = calc.evaluate(symbols, positions)
full = calc.evaluate(symbols, positions, full_derivative=True)
assert fixed.force_convention == "fixed_charge"
assert full.force_convention == "charge_response"
```

The two arrays can differ because the LAMMPS conventions use `14.4` in QEq,
`332.06371` for Coulomb energy, and `23.02` for QEq self-energy conversion.
Since `14.4 * 23.02 != 332.06371`, the QEq solution is not exactly stationary
for the reported energy. The retained [QEq audit](validation.md#qeq-audit)
checks this using fresh LAMMPS single-point calculations.

Relaxation always uses fixed-charge forces; see [optimizers](ase-relaxation.md).

## Properties

`total_bond_orders` contains per-atom row sums of `bond_orders`. `lone_pairs`
contains the per-atom lone-pair values used for comparison with LAMMPS.
`bond_counts` counts corrected bond orders strictly greater than `0.3`.

For periodic cells, `bond_orders[i, j]` sums contributions from all images of
atom `j`. The diagonal may be nonzero when an atom bonds to its own translated
copies. Bond counts apply the threshold to each individual image **before**
summing; counting entries of the aggregated matrix above `0.3` is not equivalent.

The neutral dipole is independent of the coordinate origin, but individual
atom wrapping changes its periodic coordinate branch. It is not a unique bulk
polarization.

## Branches and precision

Bond-order and interaction thresholds select discrete branches. Forces are
derivatives within the current branch; selection itself is not differentiated.
Cutoff crossings and singular geometries require care when comparing numerical
derivatives. An exactly collinear active torsion raises `ValueError`.

The calculator rejects coincident atoms, nonfinite input, unknown atom labels,
nonneutral charge, and singular QEq systems. The [API reference](api.md)
documents result fields and available options.
