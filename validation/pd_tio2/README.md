# Pd4/rutile TiO2 smoke test (awaiting the force field)

`python validation/pd_tio2/build.py` uses ASE to produce `initial.xyz`:
48 Ti, 96 O and 4 Pd atoms in three stoichiometric O–TiO–O trilayers.
The bottom trilayer (48 atoms, tag 1) is fixed. The unreconstructed rutile
(110) surface is a deliberate simplification of the published reconstructed
(011)-(2×1) surface, not a reproduction of that experiment. Bulk construction
uses a = 4.594 Å, c = 2.959 Å and oxygen coordinate u = 0.305.
The lateral cell is 12.994 × 11.836 Å. The tetrahedral Pd4 cluster has
2.75 Å edges and its base starts 2.3 Å above the highest surface oxygen.
There is 15 Å vacuum on each side of the complete slab/cluster structure.
No hydroxyls, vacancies or adsorption-site search are included.

## Parameter source and current blocker

R. Addou et al., *Influence of Hydroxyls on Pd Atom Mobility and Clustering
on Rutile TiO2(011)-2 × 1*, ACS Nano **8**, 6321–6333 (2014),
[doi:10.1021/nn501817w](https://doi.org/10.1021/nn501817w).
[SCM's catalogue](https://www.scm.com/doc/ReaxFF/Included_Forcefields.html)
identifies the corresponding complete H/O/Ti/Pd potential as `HOTiPd.ff`.

The publisher's [SI metadata](https://api.figshare.com/v2/articles/2280571)
lists only `nn501817w_si_001.pdf`, available from
<https://ndownloader.figshare.com/files/3916969>. This PDF was downloaded
and inspected: its two tables contain H and Pd adsorption energies, not
force-field parameters. No separate parameter attachment is listed.
The authors' public publication pages and searched public ReaxFF collections
did not yield the complete file. A third-party
[archive listing](https://www.dssz.com/3665348.html) lists an 8013-byte
`HOTiPd.ff`, but downloading the archive requires login; its contents and
provenance have not been verified. Search checked 2026-09-26.

The actual `HOTiPd.ff` has **not** been obtained. No replacement or inferred
parameters have been added. The xreac parse/round-trip check and LAMMPS
relaxation have therefore **not run**. There are no energy, force or relaxed
structure results. The next prerequisite is the unchanged complete file
from SCM/AMS (normally `atomicdata/ForceFields/ReaxFF/HOTiPd.ff`) or another
verifiable distribution. The structure alone does not validate this potential.
