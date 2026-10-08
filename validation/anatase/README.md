# Water on anatase (101): Monti2012

Retained inputs test 1, 2, 4 and 8 H2O on a p(1x4) anatase (101) slab:
192 slab atoms, four stoichiometric repeats, bottom 48 atoms fixed, xy
periodicity and fixed cell. Water starts 2.3 Angstrom above distinct Ti5c
sites. Four waters fill one row (0.5 ML); eight fill both rows (1 ML).
`slab_initial.xyz` and `water_initial.xyz` are the clean-slab and
isolated-water reference inputs. Files use the `.xyz` extension with extended XYZ contents (cell, periodicity
and constraints). Only inputs and this summary are retained.

Tests use the published `data/ffield.reax.TiOH.Monti2012` without numerical
edits. ASE FIRE uses dt=0.05, dtmax=0.3, maxstep=0.1 Angstrom and
fmax=0.02 eV/Angstrom; charges are re-equilibrated at each geometry and forces
follow the fixed-charge LAMMPS convention. Reference minima use fmax=0.005.

| Waters | Final Ti-O(water), Angstrom | Intact Ti-bound waters | Adsorption energy per water, eV |
| ---: | ---: | ---: | ---: |
| 1 | 2.3824 | 1/1 | — |
| 2 | 2.3779-2.3780 | 2/2 | — |
| 4 | 2.3759 | 4/4 | -1.0481 |
| 8 | 2.3923 | 8/8 | -0.9729 |

All four minimisations converged and passed initial/final LAMMPS checks of
energies, forces, charges and bonding properties. Fixed atoms did not move.
Adsorption energy is `(E(slab+nH2O)-E(clean slab)-n*E(gas H2O))/n`.
These are ordered molecular starting configurations, not a search over water
arrangements or a validation of dissociation barriers or slab-size convergence.

```sh
python examples/anatase_water.py --waters 1 2 4 8 --verify
python scripts/validate.py --system anatase --verify
```

Fresh results go to ignored `validation/runs/`; do not retain trajectories,
logs, calculation caches or output archives here. Geometry/reference tests
use Monti2012 only. Kim2013 and Ganeshan2020 remain available as parameter
files with provenance in [data/README.md](../../data/README.md).
