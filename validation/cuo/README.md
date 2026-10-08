# CuO surface validation

`surface_cuo_010.xyz` is a stoichiometric 96-atom tenorite (010) slab:
two conventional-cell layers, 2x3 repeat, 12 Angstrom vacuum on each side,
xy periodicity and seeded 0.01 Angstrom displacements. Bulk parameters are
a=4.6837, b=3.4226, c=5.1288 Angstrom, beta=99.54 degrees, C2/c,
Cu=(1/4,1/4,0), O=(0,0.416,1/4).

The unchanged `ffield.reax.CuOHCl.2010` is used; see [provenance](../../data/README.md).
ASE/native energies, forces, QEq charges and bonding properties passed
LAMMPS checks. The recorded force discrepancy was below 4.6e-10
kcal/mol/Angstrom. This fixture exposed the corrected distinction between
uncorrected `valency_val` and corrected `valency_boc` coordination.
It is unrelaxed and establishes numerical consistency, not surface stability.

```sh
python examples/cuo_surface.py
python scripts/validate.py --system cuo --verify
```

Only this input and summary are retained. Calculation/performance outputs
are written under ignored `validation/runs/`; historical timings remain in
Git history. LAMMPS verification is required for this example.
