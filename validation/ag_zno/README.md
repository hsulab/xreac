# Ag4/ZnO ReaxFF smoke test

## Force field

Literature: A. Lloyd, D. Cornil, A. C. T. van Duin, D. van Duin, R. Smith,
S. D. Kenny, J. Cornil and D. Beljonne, *Development of a ReaxFF potential
for Ag/Zn/O and application to Ag deposition on ZnO*, Surface Science
**645** (2016) 67–73, [doi:10.1016/j.susc.2015.11.009](https://doi.org/10.1016/j.susc.2015.11.009).
The [public author manuscript](https://repository.lboro.ac.uk/articles/journal_contribution/9234722)
was read via <https://ndownloader.figshare.com/files/16816238>.

The authors' [public parameter release](https://www.researchgate.net/publication/309321932_Parameter_set_for_AgZnO_ReaxFF_potential)
lists `ffield`, uploaded by Adam Lloyd in October 2016. Its download returned
a ResearchGate “Temporarily Unavailable” HTML page during this task.
With the user's explicit approval, `data/ffield.AgZnO` instead preserves
the following publicly downloadable, attributed repository copy **byte for byte**:

<https://raw.githubusercontent.com/by-student-2017/lammps_education_reaxff_win/055804c838e0d18aad321b631e138ca2e2635d9a/potentials/ffield.reax.2016X19.CHONSiCuAgZn>

The collection's
[bibliography](https://github.com/by-student-2017/lammps_education_reaxff_win/blob/055804c838e0d18aad321b631e138ca2e2635d9a/potentials/README.txt)
attributes entry `2016X19` to the Lloyd paper. Its alternative
`potentials/ffield.reax.012.CHONSiCuAgZn` has identical whitespace-separated
tokens. This establishes repository attribution, **not byte identity with
the inaccessible author-hosted attachment**. No SCM/AMS installation was
used, and no parameter was reconstructed or edited. Downloaded 2026-09-26;
27,391 bytes; SHA256:
`bad189f45e07047512a2b243f15b86f6aedadfb4d07b4d671058588a3bfdbdfb`.
No explicit redistribution license was located for this parameter file.

| Source section | Count |
|---|---:|
| General parameters | 39 |
| Atom types | 8 |
| Bonds | 29 |
| Off-diagonal entries | 20 |
| Angles | 116 |
| Torsions | 49 |
| Hydrogen bonds | 4 |

Elements, in file order: **C H O N Si Cu Ag Zn**. Simulations contain only
Zn/O or Ag/Zn/O. All original types and interactions remain in the file.

- **LAMMPS load: PASS**, including QEq and finite energy/forces.
- **xreac parse: PASS**, using `ForceField.from_file`.
- **xreac round trip: PASS**, using `ForceField.to_file` and reading again.
  All 2,711 numeric tokens (including section counts and element indices),
  retained source records, and derived model tables agree exactly.
- **LAMMPS use of the xreac-written file: PASS**. Bulk single-point energy,
  forces and charges match the original file (zero observed discrepancy;
  force and charge checks use absolute tolerances 1e-10 eV/Å and 1e-12 e).

xreac previously had no writer. The small exporter retains explicit numeric
source records, including unused fields and wildcard entries, because model
preprocessing is not reversible. `section_counts` reports these original
records rather than expanded/mixed entries. Export requires an unmodified
parsed object; it rejects edits to derived model tables instead of silently
discarding them. The original source need not remain on disk.

Source oddity: the Ag atom record gives a mass of **63.546**, retained unchanged.
LAMMPS receives standard ASE masses from the structure data (Ag 107.8682).
Its ReaxFF potential-energy evaluation does not use this Ag mass field;
the data-file masses govern FIRE's fictitious dynamics. No physical MD was run.

## Reproduction and setup

From the repository root, with xreac dependencies, ASE and `lmp_serial` installed:

```sh
python validation/ag_zno/run.py
```

This checks parsing/export and LAMMPS equivalence in a temporary directory,
then minimizes bulk ZnO and, only if bulk passes, builds and minimizes Ag4/ZnO.
It regenerates the compact retained outputs here. The actual structure data
and LAMMPS inputs are retained, so the relaxations can also be reproduced
from this directory with `lmp_serial -in bulk.in.lammps -log bulk.log.lammps`
and `lmp_serial -in in.lammps -log log.lammps`.

LAMMPS version: **22 Jul 2025 – Update 4**. Both minimizations use the
**original** `../../data/ffield.AgZnO`, `units real`, ReaxFF, QEq tolerance
1e-10, and FIRE. Convergence uses zero energy tolerance and a global force
norm tolerance equivalent to 0.01 eV/Å; final per-atom force norms are reported
below. Cell vectors remain fixed. `run 0` warnings about no integration fix
are expected for these static evaluations; both minimizers stop on force
tolerance, without QEq warnings.

Bulk reuses `scripts/validate.py:zno_bulk_case`: 128 atoms, wurtzite
a = 3.25 Å, c = 5.21 Å, u = 0.382, with its existing seeded 0.01 Å perturbation.
This tests internal relaxation at fixed lattice dimensions, not equilibrium
cell parameters or elastic stability.

The nonpolar surface uses ASE's hexagonal three-index (100), equivalent to
(10-10), with four ZnO repeats normal to the surface (eight atomic planes,
9.382 Å thick), repeated 4×3 laterally. Its 13.0×15.63 Å lateral cell contains
96 Zn + 96 O. The lowest ZnO repeat (48 atoms) is fixed, and 15 Å vacuum
is added on each side of the complete supported cluster. `clean.xyz` retains
the clean unrelaxed slab. Ag4 starts tetrahedral with 2.85 Å edges and an
exact closest surface distance of 2.50 Å. Only one starting configuration is used.

## Results

| Quantity | Bulk ZnO | Ag4/ZnO |
|---|---:|---:|
| Relaxation | PASS | PASS |
| Initial total energy (eV) | -495.23403318 | -715.70438646 |
| Final total energy (eV) | -495.50509256 | -718.39388204 |
| Maximum mobile force (eV/Å) | 0.00248255 | 0.00282792 |
| FIRE iterations | 191 | 319 |
| Minimum Zn–O (Å) | 1.973965 | 1.905512 |

Bulk nearest Zn–O distances span 1.973965–1.974059 Å, and every Zn retains
four O neighbors within 2.5 Å. The wurtzite structure remains sensible.

Ag–Ag distances (pairs 1–2, 1–3, 1–4, 2–3, 2–4, 3–4):
**2.860068, 3.098789, 2.774024, 3.097857, 2.774579, 2.769024 Å**.
Minimum **Ag–O = 2.296662 Å**, **Ag–Zn = 2.993455 Å**.

Ag4 remains a compact, distorted **3D tetrahedron**: height span 2.553 Å and
tetrahedron volume 2.821 Å³. It neither flattens nor fragments/spreads over
the surface. Its lowest atom remains 1.523 Å above the highest oxide atom,
so there is no penetration. The support preserves its layered structure;
Zn has three or four O neighbors within 2.5 Å. The largest oxide displacement
is 0.503 Å using minimum-image displacements. Bottom atoms remain fixed.
Numerical inspection and two side views show surface relaxation without
obvious collapse or suspiciously short Zn–O contacts.

## Comparison and conclusion

Lloyd et al. report that Ag prefers clustering on the **nonpolar** surface,
whereas isolated adatoms are favored on the polar surface (author manuscript,
discussion of Figure 3). They also find less dimer splitting on the nonpolar
surface and penetration during sufficiently energetic impacts. Our intact,
surface-bound cluster is qualitatively consistent with the reported nonpolar
clustering preference. This local, zero-temperature optimization does not test
deposition kinetics, adsorption energies, or the global minimum. Their polar
interface spacings are not a quantitative target for this nonpolar Ag4 case.

The public attributed file is structurally complete; xreac round-trips it
numerically, and LAMMPS accepts both original and exported versions. Bulk ZnO
and this single supported cluster remain physically sensible. No obvious
misuse is indicated by these checks. The outstanding provenance limitation
is the lack of a direct comparison with the author-hosted attachment.
