# Parameter files and licenses

Files retain their original contents except the documented PDF extraction below. Each file has its own license;
the GPL license of the xreac code does not replace the license of the data.

| File | License | Source |
| --- | --- | --- |
| `ffield.reax.HO.2015` | GPL-2.0-or-later | LAMMPS `examples/reaxff/water/qeq_ff.water` |
| `ffield.reax.CHO.2008` | GPL-2.0-or-later | LAMMPS `potentials/ffield.reax.cho` |
| `ffield.reax.ZnOH.2010` | GPL-2.0-or-later | LAMMPS `potentials/ffield.reax.ZnOH` |
| `ffield.reax.CuOHCl.2010` | **CC-BY-NC-4.0** | ACS supplement `jp102272z_si_001.txt` |
| `ffield.reax.CHOCl.2021` | **CC-BY-NC-3.0** | Hur et al., RSC supplement, pp. S10–S12 |
| `ffield.reax.PtNiCHO.2016` | **CC-BY-NC-4.0** | Gai et al., ACS supplement `jp6b01064_si_001.pdf`, pp. S7–S10 |
| `ffield.reax.PtNiCHO.2026` | **CC-BY-NC-4.0** | Experimental derivative: O–H–Pt `p_val4` changed to 1.0250 |

The first three files come from LAMMPS `stable_22Jul2025_update4`.
See `NOTICE` for their original citations, source URLs, and checksums.

## Gai Pt/O/H parameters

L. Gai, Y. K. Shin, M. Raju, A. C. T. van Duin, and S. Raman,
“Atomistic Adsorption of Oxygen and Hydrogen on Platinum Catalysts by Hybrid
Grand Canonical Monte Carlo/Reactive Molecular Dynamics,” J. Phys. Chem. C
120, 9780–9793 (2016), [DOI:10.1021/acs.jpcc.6b01064](https://doi.org/10.1021/acs.jpcc.6b01064).

Extracted from [the publisher's supplement](https://ndownloader.figshare.com/files/5029120),
section 3, pp. S7–S10. The original header says “Pt/Ni/C/H/O force field 2015”;
the filename uses the publication year 2016. All C/H/O/Ni/Pt/X types are retained.
This is the original Pt/O/H publication's supplement, not the later Pt–Ni alloy
paper (DOI:10.1021/acs.jpca.6b06770).

[Publisher metadata](https://api.figshare.com/v2/articles/3204151) assigns
**CC BY-NC 4.0**. Attribution and noncommercial terms apply separately from
the code license; see `LICENSES/CC-BY-NC-4.0.txt` and `NOTICE`.
Extraction removes page numbers/whitespace and restores one missing separator
between `5.0000` and `9999.9999` in dummy X. No numerical values are changed.

```sh
python scripts/extract_gai_parameters.py jp6b01064_si_001.pdf data/ffield.reax.PtNiCHO.2016
```

Source PDF SHA256: `ab99e02eb876347637c084b1741028851513e640b2879920483fe42a4bb5fc71`.
Parameter SHA256: `a9769ca030d0a94b3a9d1d09acceb1d6e46b283e2fc701adb5aaaee3cc19588b`.
The source and metadata are retained in `validation/water/surfaces/gai2016_source.tar.gz`.
Pt(111)/water test results are described in `validation/water/README.md`.

### 2026 experimental revision

`ffield.reax.PtNiCHO.2026` is an xreac experimental derivative,
not the published parameter set. Only the O–H–Pt angle's bond-order exponent
`p_val4` changes from 1.0000 to 1.0250; all other numerical entries are identical.
The year 2026 identifies this local revision of Gai2016, not a new publication.
The change reduces the weak-bond cutoff force cusp that stalled water
relaxation on Pt(111). The original `ffield.reax.PtNiCHO.2016` is preserved.
Its header identifies the modification and source. It retains the attribution
and CC BY-NC 4.0 terms above. SHA256:
`4c4a54ee24565faf8fdfe0269da6b352e45ec7658bc2a64152a52e792be80209`.

The change permits the tested Pt(111)/water endpoints to reach 0.02 eV/Angstrom
and the seven-image CI-NEB to reach 0.05 eV/Angstrom. It was selected for
numerical convergence, **not fitted to DFT energies or barriers**. Residual
cutoff sensitivity remains; see the retained tuning report in
`validation/water/README.md`. Select it explicitly with
`ForceField.bundled("ffield.reax.PtNiCHO.2026")`.

The retained [Pt/Ni examples](../examples/README.md#water-dissociation-on-pt111-and-ni111)
use two-layer p(2x2) slabs, a fixed bottom layer at z=2.5 Angstrom, 15 Angstrom
total vacuum padding, and seven NEB images. Ni/O/H uses the unchanged Gai2016
entries; no Assowe2012 file is bundled.

| Example | NEB/CI steps | NEB runtime | Barrier |
| --- | ---: | ---: | ---: |
| Ni(111) | 134 | ~7 s | 0.52115 eV |
| Pt(111) | 1,810 | ~84 s | 0.98303 eV |

Both endpoints and bands meet their force targets. Ni passes all LAMMPS
comparisons; Pt retains a product equivalent-copy force-consistency failure.
These numerical tests do not establish DFT accuracy.

- [Validation notes and reproduction commands](../validation/water/README.md#two-layer-pt111-and-ni111-water-dissociation)
- [Ni results](../validation/water/surfaces/ni_2026_7images_lowered.json) and [structures](../validation/water/surfaces/ni_2026_7images_lowered_band.extxyz)
- [Pt results](../validation/water/surfaces/pt_2026_7images_lowered.json) and [structures](../validation/water/surfaces/pt_2026_7images_lowered_band.extxyz)
- [Energy curves](../validation/water/surfaces/metal_2026_7images_lowered.png)
- [Compact tuning provenance](../validation/water/surfaces/pt_2026_parameter_tuning.json)

## Cu/O/H/Cl parameters

Attribution: A. C. T. van Duin, V. S. Bryantsev, M. S. Diallo,
W. A. Goddard, O. Rahaman, D. J. Doren, D. Raymand, and K. Hermansson,
“Development and Validation of a ReaxFF Reactive Force Field for Cu Cation/Water
Interactions and Copper Metal/Metal Oxide/Metal Hydroxide Condensed Phases,”
Journal of Physical Chemistry A 114, 9507–9514 (2010).

- [Paper](https://doi.org/10.1021/jp102272z)
- [Supporting information](https://doi.org/10.1021/jp102272z.s001)
- [Publisher license metadata](https://api.figshare.com/v2/articles/2732188)
- [Original file](https://ndownloader.figshare.com/files/4408471)
- [CC BY-NC 4.0 terms](https://creativecommons.org/licenses/by-nc/4.0/)

**This parameter file carries a noncommercial restriction.** Redistribution
and use must follow its own license, including attribution and notices of
changes. Bundling it with xreac does not grant commercial-use permission.
The xreac code remains GPL-2.0-or-later.

The file was renamed to `ffield.reax.CuOHCl.2010`; its contents, including
the original header and line endings, are unchanged. SHA256:
`8b1a57a6945be8b1d9df3b69df8f329dc32561393ed1f02221a5ef0c21e8c07e`.

The full license and its disclaimer of warranties are included in
`LICENSES/CC-BY-NC-4.0.txt` in source distributions and in the wheel's
`.dist-info/licenses/LICENSES/` directory. Full attribution is also retained
in `NOTICE` in both distribution formats.

## Hur C/H/O/Cl parameters

`ffield.reax.CHOCl.2021` contains the C/H/O/Cl parameters (including dummy X)
from J. Hur, Y. N. Abousleiman, K. L. Hull, and M. J. Abdolhosseini Qomi,
“Reactive force fields for modeling oxidative degradation of organic matter in
geological formations,” RSC Advances 11, 29298–29307 (2021).

- [Paper](https://doi.org/10.1039/D1RA04397H)
- [Supplement, pp. S10–S12](https://www.rsc.org/suppdata/d1/ra/d1ra04397h/d1ra04397h2.pdf)
- [Europe PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9040638/supplementaryFiles)
- [Article license metadata](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9040638/fullTextXML)
- [CC BY-NC 3.0 legal terms](https://creativecommons.org/licenses/by-nc/3.0/)

Copyright 2021 The Royal Society of Chemistry. This data is distributed under
**CC BY-NC 3.0**, including a noncommercial restriction. It is not relicensed
under the code's GPL. The complete license is retained in
`LICENSES/CC-BY-NC-3.0.txt` and included in distributions.

Extracted from archive member `RA-011-D1RA04397H-s002.pdf`, SHA256
`867348b9395826206919eb3922ff5dba2f1a3301c307f0225e41b5a995fb1cfb`.
The bundled file adds an attribution header and removes PDF page furniture,
blank lines, and outer whitespace. Numerical entries and their order, including
repeated angle/torsion entries and dummy X, are unchanged. Both archived
supplements yield identical extracted parameter text. Bundled file SHA256:
`d468f31f71c1ca429aef9e576bf2c5fda3d8eb8ad0c7c4ae433197692ee5825b`.

Reproduce with Poppler's `pdftotext` and:

```sh
python scripts/extract_hur_parameters.py RA-011-D1RA04397H-s002.pdf data/ffield.reax.CHOCl.2021
```

The source parameter pages and provenance are archived in
`validation/cho/hur2021_source.tar.gz`. Numerical comparisons and raw LAMMPS
runs are retained alongside it. The original fit targets organic oxidation by
oxychlorine species; numerical agreement does not establish an accurate
chloride–chloromethane SN2 barrier.
