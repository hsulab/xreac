# Single-point charge-equilibration audit

Every geometry in this audit is evaluated in a fresh `lmp_mpi` process with
zero initial charges and `fix charges all qeq/reaxff 1 0 10 1e-12 reaxff
maxiter 2000`, followed by `run 0`. LAMMPS performs QEq during the run setup,
before evaluating forces. No nuclear relaxation is performed.

Reference version: **LAMMPS 22 Jul 2025, Update 4**. The water force field
and exact input/output files are retained in each run directory. The
[summary](summary.json) contains the numerical comparisons and parameter checksum.

Two independent controls check that the initial single point has equilibrated
charges: running five further steps with QEq but without a time-integration fix
(positions stay fixed), and starting from O = -0.8 e / H = +0.4 e instead of zero.
Neither changes the energy, charges, or forces beyond floating-point noise.
The maximum difference between atomic QEq chemical potentials and their mean is
8.9e-16 eV for the monomer and 5.0e-14 eV for the dimer. This checks the constrained
QEq stationarity condition directly with the charges returned by LAMMPS.

The force component below is the **y force on atom 1 (oxygen)** in each saved
geometry. Units are kcal/mol/Angstrom. Central finite differences use energies
from independent, freshly charge-equilibrated LAMMPS runs at y +/- h.

| Calculation | Monomer | Dimer |
| --- | ---: | ---: |
| LAMMPS reported force | 14.0073432764 | 14.4687110648 |
| LAMMPS energy finite difference, h = 1e-3 A | 14.1809861221 | 14.6524929715 |
| LAMMPS energy finite difference, h = 1e-4 A | 14.1812556819 | 14.6527624469 |
| LAMMPS energy finite difference, h = 1e-5 A | 14.1812583777 | 14.6527651566 |
| Python full energy-gradient force | 14.1812584046 | 14.6527651684 |

Thus the discrepancy persists with equilibrated charges, and is visible using
LAMMPS's own energies and forces. It is not caused by skipped QEq in `run 0`.
The finite-difference error relative to the Python full gradient falls to
2.7e-8 and 1.2e-8 kcal/mol/Angstrom, respectively.

For this reference version, QEq uses an off-diagonal coefficient of 14.4,
the self energy uses 23.02, and the Coulomb energy uses 332.06371.
Their mismatch is 332.06371 - 23.02 * 14.4 = 0.57571.
Writing S for the shielded, tapered Coulomb matrix, the omitted charge-response
force contribution is `-(332.06371 - 23.02*14.4) * (S @ q) dot (dq/dR)`.
Using charge finite differences returned by LAMMPS predicts +0.1739151282
(monomer) and +0.1840541040 (dimer), explaining the observed force differences
within the finite-difference errors.

Sources for the tested release:

- [QEq setup and solver](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/REAXFF/fix_qeq_reaxff.cpp)
- [ReaxFF constants](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/REAXFF/reaxff_defs.h)
- [Electrostatic energy and force evaluation](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/REAXFF/reaxff_nonbonded.cpp)

Reproduce in a new directory:

```sh
python scripts/audit_qeq.py --output validation/runs/qeq-audit-new
python -m pytest -q tests/test_reference.py -k reequilibrated
```
