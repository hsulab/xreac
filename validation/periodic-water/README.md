# Periodic water verification

All five cases pass against LAMMPS 22 Jul 2025, Update 4. Generated September
23, 2026 with `mamba run -n catorch3 python examples/periodic_water.py --verify
--output validation/periodic-water`. Reproduce using a new output directory.

See [summary.json](summary.json) for signed energy differences, maximum
absolute discrepancies, acceptance tolerances, timings, and environment.
Each case contains structures (JSON and extended XYZ), complete Python and
reference properties, and the raw LAMMPS run with the exact parameter file.
Energies use kcal/mol and forces use kcal/mol/Å. All comparisons use
`full_derivative=False`; charges are equilibrated at each single-point geometry.

| Case | Atoms | Python total energy (kcal/mol) | Maximum force difference (kcal/mol/Å) |
| --- | ---: | ---: | ---: |
| Boundary-crossing dimer | 6 | -498.827197208 | 1.6e-12 |
| Multiple nonbonded/H-bond images | 6 | -492.853710805 | 5.3e-13 |
| Rotated triclinic water | 24 | -1989.024212472 | 2.9e-10 |
| Partially periodic water | 24 | -1989.403415674 | 3.9e-10 |
| Bulk water box | 192 | -16017.189296171 | 8.0e-10 |

The bulk case has 64 water molecules in a cubic 12.48 Å box, approximately
1 g/cm³. These are seeded, nonoptimized test structures, not equilibrated
liquid snapshots. Atom wrapping intentionally splits molecules across faces.
Periodic cell heights exceed 10 Å; QEq, vdW, Coulomb, and hydrogen bonding sum
all contributing images. Dipoles refer to the supplied coordinate branch;
they are not unique bulk polarizations.

The full regression suite passes 136 tests, including both force derivatives,
supercell replication, lattice wrapping, symmetry checks, and fixed-charge
periodic relaxation with ASE and native FIRE.
