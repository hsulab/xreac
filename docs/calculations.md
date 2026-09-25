# Energies, charges, and forces

## Charge equilibration

Every evaluated geometry receives a fresh QEq solution. A dense linear
system equalizes the charge chemical potentials while a Lagrange multiplier
enforces `sum(charges) = total_charge` (default zero, in e). Previous or initial atomic charges are not frozen
across geometries.

For small cells, translated copies share the input atom's charge. QEq includes
coupling to nonzero images of the same atom and solves an `(N+1) × (N+1)` system
for the `N` input atoms. Global QEq can transfer charge between separated
fragments; it is not a guarantee of physically correct dissociation.

## Charged systems

Use `calc.evaluate(symbols, positions, total_charge=-1)` or
`calc.relax(symbols, positions, total_charge=-1)`. ASE uses
`ReaxFFCalculator(ff, total_charge=-1)`. Fractional charges are accepted.
See `examples/charged_water.py` for hydroxide and hydronium demonstrations.

For periodic systems the charge belongs to the input cell. Internal neighbor
replication preserves it; explicitly repeating the physical cell requires
multiplying `total_charge` by the number of copies. Electrostatics retain the
existing shielded, finite-cutoff model, with no Ewald sum or compensating
background. Numerical support does not establish accurate ionic chemistry or
equivalence to charged periodic DFT.

Bučko's *Ab initio calculations of free-energy reaction barriers*
([2008, DOI](https://doi.org/10.1088/0953-8984/20/6/064211)) studies
Cl⁻ + CH₃Cl → CH₃Cl + Cl⁻, with total charge −1. This is a possible future
application. The bundled Hur `ffield.reax.CHOCl.2021` covers C/H/Cl, but using
it for this reaction still requires validation of fragment charges and reaction
energies, plus free-energy sampling.
No reproduction of its barriers is claimed. Only a global charge constraint
is supported; fragment constraints and spin states are not modeled.

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

Dipoles are evaluated relative to the force-field center of mass:
`sum(q[i] * (x[i] - center_of_mass))`. This convention is translation invariant
also for charged systems. A charged dipole depends on the choice of physical
origin; xreac fixes that choice to the center of mass. Individual atom wrapping
changes the periodic coordinate branch and the center of mass. It is not a unique bulk
polarization.

## Branches and precision

Bond-order and interaction thresholds select discrete branches. Forces are
derivatives within the current branch; selection itself is not differentiated.
Cutoff crossings and singular geometries require care when comparing numerical
derivatives. An exactly collinear active torsion raises `ValueError`.

The calculator rejects coincident atoms, nonfinite input, unknown atom labels,
invalid total charge, and singular QEq systems. The [API reference](api.md)
documents result fields and available options.
