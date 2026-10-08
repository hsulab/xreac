# Pd/Ce/O parameter provenance audit

Search date: 2026-09-26. No force-field values were invented, inferred from
defaults, fitted, or changed. A complete candidate was located, but it cannot
yet be identified unambiguously as the published combined potential.

## Published sources

1. Broqvist, Kullgren, Wolf, van Duin, and Hermansson, J. Phys. Chem. C
   **119**, 13598–13609 (2015), [DOI](https://doi.org/10.1021/acs.jpcc.5b01597).
   [NIST record](https://www.ctcms.nist.gov/potentials/entry/2015--Broqvist-P-Kullgren-J-Wolf-M-J-et-al--Ce-O/)
   says Kullgren supplied the file on 2016-12-19. The current `ipr2` version
   separates comments for LAMMPS compatibility (2020-03-15), without a reported
   reparameterization. Its internal header is dated 2016-02-09.
   [Parameter download](https://www.ctcms.nist.gov/potentials/Download/2015--Broqvist-P-Kullgren-J-Wolf-M-J-et-al--Ce-O/2/ffield_Ce-O)
   and [author comments](https://www.ctcms.nist.gov/potentials/Download/2015--Broqvist-P-Kullgren-J-Wolf-M-J-et-al--Ce-O/2/comments.txt)
   are preserved byte-for-byte in `data/ffield.reax.CeO.Broqvist2015*`.
2. Senftle, van Duin, and Janik, ACS Catalysis **7**, 327–332 (2017),
   [DOI](https://doi.org/10.1021/acscatal.6b02447), published online December 2016.
   The [paper on the author's website](https://senftle.blogs.rice.edu/files/2023/05/16_ACSCatal_PdCeria_ReaxFF.pdf)
   directs readers to its Pd/O/C/H supplement, to the Ce/O authors for that
   file, and to Penn State MCC for the complete parameter sets.
   [SI 1](https://ndownloader.figshare.com/files/7013171), pp. 5–6,
   gives cross-term Tables S2–S4. [SI 2](https://ndownloader.figshare.com/files/7013168)
   is a complete Pd/O/C/H file, with no Ce. Both PDFs are retained in `data/`.
   [SI 1 metadata](https://api.figshare.com/v2/articles/4300295) and
   [SI 2 metadata](https://api.figshare.com/v2/articles/4300292) specify CC BY-NC 4.0.
   The Pd/O lineage is Senftle et al., JCP **139**, 044109 (2013),
   [DOI](https://doi.org/10.1063/1.4815820); the **actual study-supplied SI 2**
   is the preferred Pd/O source over an independently merged older binary file.

The extracted SI 2 text was obtained with `pdftotext -raw`, removing the cover
and form-feed characters and folding numeric tokens into four-line atom,
two-line bond, and single-line other interaction records. Section counts are
39 general, 5 atom types, 10 bonds, 6 off-diagonals, 36 angles, 25 torsions,
1 hydrogen bond. Original numeric tokens and their order are unchanged.

## Comparison before any merge

| Shared block: NIST Ce/O versus study SI 2 | Result |
| --- | --- |
| General parameters | All 39 identical |
| O atomic parameters | All 32 identical |
| O–O bond | All 16 identical |
| O–O–O–O torsion | All 7 identical |
| Other duplicated interactions | None |

SI 2 additionally contains O–O–O angle parameters; NIST Ce/O has no such row.
This is an absence, not a numerical mismatch. Atomic masses and EEM fields
were included in the comparison. Interaction indices were mapped to element
names before comparing, including reversed-equivalent bond/angle/torsion keys.

Tables S2–S4 give only 6/16 fields for the Pd–Ce bond, 4/6 for the Pd–Ce
off-diagonal, and 5/7 per cross angle. In particular they do not specify
the bond-order correction switches or angle penalty/conjugation fields.
The standalone Ce/O and Pd/O files have no Pd–Ce rows to copy these from.

## Complete public candidate

`data/ffield.candidate.PdCHOCe.matsci2025` preserves
[`PdCHOCe.field`](https://matsci.org/uploads/short-url/t0eV52HsvTgbgah1W9wQRq9j8nz.field)
from [yfyin's forum post, March 25, 2025](https://matsci.org/t/first-nvt-then-nve-energy-is-not-constant/62251/3).
The same uploader reused the identical attachment URL in an
[April 20, 2025 post](https://matsci.org/t/compute-cn-pd-coord-atom-cutoff-3-5-group-pd/63098).
These are not two independent sources. Neither supplies an author provenance
chain, publication citation for the file, or explicit data license.
Header: `Reactive MD-force field: Oct15 2012 Ceria  Pd Added April23`.
That header does not establish the year of the Pd addition.

All 39 general parameters, all source atomic parameters, and all **78**
interaction rows of SI 2 match this candidate numerically. All **25** cross
values in Tables S2–S4 also match. Its full Pd–Ce bond is:

```text
68.3117 0.0000 0.0000 -0.3467 -0.2000 0.0000 16.0000 0.1073
11.8297 -0.2000 15.0000 1.0000 -0.1205 5.8218 0.0000 0.0000
```

Its off-diagonal is `0.2452 2.2178 12.8866 2.3042 -1.0000 -1.0000`.
Its complete cross angles (standard seven-column order) are:

```text
O-Pd-Ce  71.0692 9.9527 7.9620 20.0000 0.1000 0.0000 1.0102
O-Ce-Pd   1.8709 2.1607 7.9002 20.0000 0.1000 0.0000 1.0116
Ce-O-Pd  90.0000 6.7405 7.9834 90.0000 1.0000 0.0000 2.1474
```

The fourth column, `p_pen1`, is **20, 20, 90**, not zero. These are observations
from the candidate, not authenticated recoveries of the study's omitted fields.

The candidate's Ce and O atomic rows, Ce–Ce and O–O bonds, and OOOO torsion
match NIST. The following **26 values in six interaction rows differ**.
Indices are one-based positions in each numeric parameter row, excluding atom IDs:

| Row | Field indices | NIST Ce/O | Complete candidate |
| --- | --- | --- | --- |
| Ce–O bond | 1, 4, 8, 9, 13, 14 | 162.8952, -0.4477, 0.2661, 2.0330, -0.0610, 5.1946 | 163.2707, -0.4094, 0.2509, 1.8415, -0.0810, 5.3892 |
| Ce–O off-diagonal | 1, 2, 3, 4 | 0.3100, 1.8000, 11.6966, 1.6000 | 0.2752, 1.9083, 11.9233, 1.6089 |
| O–Ce–O angle | 1, 2, 3, 7 | 72.5707, 26.3851, 0.1355, 3.2470 | 71.9708, 24.9111, 0.0894, 3.2897 |
| Ce–O–Ce angle | 1, 2, 3, 7 | 79.7353, 5.3657, 7.0000, 1.0000 | 78.7111, 10.1184, 7.1795, 1.0050 |
| O–Ce–Ce angle | 1, 2, 3, 7 | 48.5535, 34.8873, 0.5643, 4.0000 | 50.2725, 23.7234, 0.0953, 3.8461 |
| O–O–Ce angle | 1, 2, 3, 7 | 64.6845, 16.7295, 1.0234, 1.0000 | 66.7976, 17.2897, 0.9855, 1.0381 |

The exact cross-value agreement is evidence of a relationship to the published
fit, but does not establish whether the Ce/O differences are the study's actual
branch, an earlier/later variant, or an undocumented modification.
Copying the candidate's cross rows into NIST Ce/O would create a hybrid that
neither complete source supplies. No such merge has been made.

## Additional searches and branch checks

- [Public ReaxFF collection](https://github.com/by-student-2017/lammps_education_reaxff_win),
  tree `055804c838e0d18aad321b631e138ca2e2635d9a`: its
  `potentials/ffield.reax.2017X2.CHOPd` is associated with the target DOI in
  `potentials/README.txt`. All 78 interaction rows match SI 2, but it has no Ce.
  Thus it is an independent mirror, not the missing combined potential.
- [ReaxCeriaEfield](https://github.com/warisa-r/ReaxCeriaEfield), tree
  `82baf7d2c7e9d18198735a1c3761a1e57543f417`, file
  `simulations/energy_vacancy_sim/ffield.reax`: header also mentions Oct15 2012
  ceria. All nine NIST Ce/O interaction rows agree exactly; it contains no Pd.
  This does not support replacing NIST with the candidate's differing Ce/O rows.
- [Senftle group website](https://senftle.blogs.rice.edu/) and its linked
  [GitHub account](https://github.com/tsenftle/): checked publication downloads
  and public repository inventory. `ReaxFF-GCMC` at
  `86334cb807e36a96f8b566825ecc7e2a5bbdac8c` contains Li/Si/C/N/H work;
  only its master branch and no tags were listed. `Coking-GCMC` at
  `8698075c6509cca0c096052dc142cb19613e96d8` contains Fe/C GCMC code.
  No Pd/Ce/O source found in these inventories.
- [Penn State MCC](https://www.mri.psu.edu/affiliated-institutes/computational-materials)
  provides ReaxFF project/contact information, not a located combined download.
  No authors or forum users have been contacted.
- [SCM catalogue](https://www.scm.com/doc/ReaxFF/Included_Forcefields.html)
  lists CeO and PdO separately, not the target combined field.
  [LAMMPS potentials](https://github.com/lammps/lammps/tree/develop/potentials)
  has no Pd/Ce/O ReaxFF file.
- Additional GitHub inventories checked: `annashchygol/ReaxFF_parameters_DataBase`
  (`879d646cca0ba9057c6b7289549b235da552ef76`),
  `mikepols/reaxff-force-fields` (`72d63e3fc9033102b33307490b8f1b33150b12bb`),
  and `fenggo/I-ReaxFF` (`d7cc28869fbe1c766f5e8bab0a78009284c35456`).
  No named Pd/Ce/O candidate located. Filename inventories are not exhaustive
  full-content searches of every repository history.
- Figshare API searches for `Pd CeO2 ReaxFF` returned only the two target
  supplements; `palladium ceria ReaxFF` returned none. A Zenodo search for
  `ReaxFF AND (ceria OR palladium)` found no relevant record.
- Earlier [2015 Pd/ceria study](https://doi.org/10.1021/acscatal.5b00741):
  downloaded and inspected [its supplement](https://ndownloader.figshare.com/files/3759709),
  which contains DFT structures/energies and no complete ReaxFF parameter file.
  Later [2018 CO oxidation study](https://doi.org/10.1021/jacs.7b13624)
  cites the target ReaxFF work but presents a DFT study, not a located reused
  force-field archive. Searches for subsequent Pd/ceria ReaxFF work and exact
  candidate names/numeric signatures did not identify an author-linked copy.
- Penn State dissertation catalogue search was blocked by its JavaScript
  traffic check; grep.app code search returned an access-check page.
  These unsuccessful accesses are not evidence that a file does not exist.

Outcome: **one complete but unauthenticated candidate, no independently
verified complete published parameterization**. The required next provenance
step is confirmation of the candidate's full cross rows and the Ce/O version
used by the study. A successful minimization alone would not resolve this.

## SHA256 of retained files

```text
527dfadb8b60f32814e582dfe8d90afb224ff5a7c0b2746b1a490fde0f5a338a  ffield.reax.CeO.Broqvist2015
938de3a5784b92b8cb74c293d1ad108469612309aac7c5c0b1e17083cd5567d0  ffield.reax.CeO.Broqvist2015.comments
1db3d1f9b81f8296d0acb3d2e242b556369b6ff2dedaac3b71aa1b04c993c083  ffield.reax.PdOCH.Senftle2017
2a056e3d494e18b6ab8139e71ddc5240a08640d3fcc8bf73ff6fb84e05ca8c28  PdCeO.Senftle2017.cross-parameters.pdf
440d1543743fda332d98abbc8e0efcd081ec2f0559a7006e821b0b7bd88233d5  PdOCH.Senftle2017.parameters.pdf
d61e89925accfcb81bc8e465c783cd4b26ff1f41da29571070682872ee5aba34  ffield.candidate.PdCHOCe.matsci2025
```
