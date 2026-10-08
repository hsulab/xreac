# Ag4/ZnO smoke test

Inputs are `bulk.initial.xyz` (128-atom perturbed wurtzite ZnO), `clean.xyz`
(192-atom nonpolar ZnO(10-10) slab), and `initial.xyz` (supported tetrahedral
Ag4). The slab fixes its bottom 48 atoms, uses xy periodicity and 15 Angstrom
vacuum per face. Ag4 starts with 2.85 Angstrom edges and 2.50 Angstrom minimum
surface separation. Only one starting arrangement was tested.

The unedited `data/ffield.AgZnO` is a public copy attributed to Lloyd et al.,
Surface Science 645 (2016) 67-73, DOI:10.1016/j.susc.2015.11.009. Identity with
the inaccessible author attachment and redistribution license remain
unverified. Its Ag mass entry is 63.546; LAMMPS data use ASE atomic masses.
See [parameter attribution](../../data/NOTICE).

Parsing, lossless export and original/exported LAMMPS comparisons passed.
Fixed-cell LAMMPS FIRE converged for both systems:

| System | Initial/final energy, eV | Max mobile force, eV/Angstrom | Steps |
| --- | --- | ---: | ---: |
| Bulk ZnO | -495.234033 / -495.505093 | 0.002483 | 191 |
| Ag4/ZnO | -715.704386 / -718.393882 | 0.002828 | 319 |

Bulk Zn retains four O neighbours. Ag4 remains a compact distorted 3D cluster,
with minimum Ag-O=2.29666 and Ag-Zn=2.99346 Angstrom; no penetration or
support collapse was observed and fixed atoms remained unchanged. This is a
local relaxation, not a deposition, global-minimum or adsorption-energy test.

```sh
python examples/ag_zno.py --lammps lmp_serial
python examples/ag_zno_build.py --output validation/runs/ag-zno-build
```

Only inputs and this summary are retained. Runners write generated outputs
to fresh directories under ignored `validation/runs/`.
