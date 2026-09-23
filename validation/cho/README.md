# C/H/O validation

All 4 retained cases pass fresh LAMMPS comparisons. The largest force
difference is 9.322e-12 kcal/mol/Å. Both Python neighbor builders agree.

| Case | Atoms | Checks |
| --- | ---: | --- |
| `methane` | 5 | C–H bonding and valence |
| `carbon_monoxide` | 2 | Heteronuclear C–O interactions |
| `carbon_dimer` | 2 | C2 lone-pair correction and triple-bond stabilization |
| `periodic_carbon_chain` | 3 | Active torsions across periodic images |

See [summary.json](summary.json), [structures.json](structures.json),
[full numerical results](results.json.gz), and [raw LAMMPS runs](reference.tar.gz).
The [independent baseline](baseline.json.gz) predates energy-model consolidation.

```sh
python scripts/validate.py --verify --system cho
```
