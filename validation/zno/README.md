# Zn/O validation

Inputs cover an isolated Zn atom, ZnO dimer, O4 many-body terms, a 20-atom
cluster, a small periodic ZnO cell and a zinc self-image chain, using
`ffield.reax.ZnOH.2010`. All six routine cases passed fresh LAMMPS checks;
the recorded largest force difference was 1.505e-10 kcal/mol/Angstrom.
ASE/native neighbor builders agree. Dimer relaxation reuses the same input.

`periodic_bulk_zno_128.xyz` is the optional perturbed wurtzite fixture:
a=3.25, c=5.21 Angstrom, u=0.382, 4x2x2 orthorhombic repeats and seeded
0.01 Angstrom displacements. This fixed-cell input tests bulk periodicity
and performance, not equilibrium lattice parameters.

```sh
python scripts/validate.py --system zno --verify --include-bulk
python scripts/validate_relaxation.py --backend ase
python scripts/benchmark.py --compare-neighbors
```

Only XYZ inputs and this summary are retained. Timings, trajectories and
numerical output archives are generated under ignored `validation/runs/`;
historical performance records remain in Git history.
