# C/H/O and C/H/O/Cl validation

The four CHO2008 inputs test methane, CO, C2 lone-pair/triple-bond terms and
periodic carbon-chain torsions. All passed ASE/native/LAMMPS checks; the
recorded maximum force difference was 9.322e-12 kcal/mol/Angstrom.
The five `chocl_*.xyz` fixtures use Hur2021: three neutral molecules and
two charged SN2 geometries. All numerical comparisons passed. Nonneutral
LAMMPS checks supply xreac charges; QEq itself is validated independently.

```sh
python scripts/validate.py --system cho --include-chocl --verify
```

`sn2_path.xyz` retains only the 27 geometries used to initialise the charged
Cl- + CH3Cl NEB (total charge -1). It is a constrained scan path, not 27
independent equilibrium structures. The recorded constrained rise was
124.842 kJ/mol. A nine-image charge-response CI-NEB converged to 0.00878
eV/Angstrom with a 125.255 kJ/mol barrier; its peak passed LAMMPS checks.
A noisier path failed to converge and skipped an intervening high barrier.
No verified smaller barrier was found. Midpoint curvature and alternate-path
concerns remain, so these are exploratory potential energies, not validated
activation free energies or proof of a first-order saddle.

```sh
python scripts/neb_sn2.py --fmax 0.01 --output validation/runs/sn2-neb
python scripts/scan_sn2.py --output validation/runs/sn2-scan
```

The NEB defaults to the retained XYZ path; `--scan` also accepts a fresh
scan JSON. Only inputs and this summary are retained. Generated calculations
go under ignored `validation/runs/`; parameter provenance stays in `data/`.
