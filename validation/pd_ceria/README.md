# Pd4/CeO2(111) starting structure

Status: **relaxation not run; published potential identity unresolved**.
The [parameter audit](parameters.md) records a complete public candidate,
its agreement with the published cross terms, and its conflicting Ce/O rows.
No missing fields have been guessed or filled. There is deliberately no
`data/ffield.PdCeO`, `relaxed.xyz`, or `log.lammps` yet.

`build.py` uses ASE fluorite CeO2 at an initial lattice constant of 5.411 Å.
It cuts (111), transfers the terminal O plane to the bottom to obtain three
complete O–Ce–O trilayers, and builds 12 primitive surface cells in an
11.47846 × 13.25419 Å rectangular surface box. This modest enlargement from
the requested approximate 3×3 surface avoids the 9.94064 Å box height of a
literal 3×3 rhombus, below the 10 Å cutoff required by
[LAMMPS QEq](https://docs.lammps.org/fix_qeq_reaxff.html#restrictions).

The model is Ce36O72Pd4 (112 atoms), with the bottom 36 atoms fixed and the
upper two trilayers plus Pd4 movable. The tetrahedron has six 2.75000 Å
edges, a triangular base 2.2 Å above the upper oxygen plane, and an apex
2.24537 Å above its base. There is 15 Å vacuum on either side of the complete
slab/cluster system, with periodicity only in x/y. No adsorption-site search
was performed. `initial.xyz` is extended XYZ and includes cell, periodicity,
trilayer tags, and the fixed-atom constraint.

Initial geometry checks (minimum-image distances):

| Quantity | Initial value |
| --- | ---: |
| Pd–Pd distances (all six) | 2.75000 Å |
| Minimum Pd–O | 2.22393 Å |
| Minimum Pd–Ce | 3.09674 Å |
| Minimum Ce–O | 2.34303 Å |
| Initial/final energy | Not evaluated |
| Maximum force | Not evaluated |

There are nine atomic planes with 12 atoms each, arranged O–Ce–O three times;
Ce:O = 1:2. The initial tetrahedron is intact and entirely above the support.
There are no initial 1.89 Å Ce–O contacts. The Ce/O authors specifically
warn about a false minimum near that distance in reduced ceria; see the
retained `data/ffield.reax.CeO.Broqvist2015.comments`.
No conclusions about relaxation, flattening, penetration, or support
reconstruction can be drawn before an authenticated-potential calculation.

Rebuild from the repository root:

```sh
python validation/pd_ceria/build.py
```

Once a source-resolved `data/ffield.PdCeO` exists:

```sh
python validation/pd_ceria/relax.py --lammps lmp_serial
```

The runner uses LAMMPS ReaxFF, neutral-system QEq (including fixed atoms),
fixed-cell FIRE minimization, and a 0.02 eV/Å maximum mobile-atom force target.
It retains only the final extended XYZ and LAMMPS log, and prints energies,
force, pair distances, Pd height span, and vertical separation from oxygen.
It exits unsuccessfully if the force target is missed. The full simulation
path is **not yet exercised**; only syntax, build geometry, and the
missing-potential failure path have been checked. The build was executed
with ASE 3.27.0 in the existing `catorch3` environment.
