# Pd4/CeO2(111) input

**No relaxation performed: the published combined potential remains unresolved.**
The [parameter audit](../../data/PdCeO.parameters.md) documents conflicting
Ce/O rows in the public candidate. No missing parameters were inferred and
`data/ffield.PdCeO` has not been created.

`initial.xyz` contains Ce36O72Pd4: three complete O-Ce-O trilayers, 12
primitive surface cells, bottom 36 atoms fixed, xy PBC and 15 Angstrom vacuum
per face. The rectangular cell is 11.47846 x 13.25419 Angstrom; the initial
fluorite lattice constant is 5.411 Angstrom. Pd4 has 2.75 Angstrom edges and
its base is 2.2 Angstrom above surface O. Minimum initial Pd-O/Pd-Ce/Ce-O
distances are 2.22393/3.09674/2.34303 Angstrom. Geometry checks passed;
energies, forces and relaxed stability have not been evaluated.

```sh
python examples/pd_ceria_build.py --output validation/runs/pd-ceria/initial.xyz
# Only after the authentic combined parameter file is available:
python examples/pd_ceria.py --lammps lmp_serial
```

The runner targets 0.02 eV/Angstrom mobile forces and refuses to run without
the potential. Its physical relaxation path remains untested. Only this
input and summary are retained; fresh outputs belong in `validation/runs/`.
