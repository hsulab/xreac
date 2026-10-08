# Validation inputs and results

Retained validation directories contain only input structures (`.xyz`, with
extended XYZ metadata) and brief READMEs describing tests and main results.
Calculation outputs, trajectories, plots, caches and archives belong in
ignored `validation/runs/`. Historical detailed records remain in Git history.
Parameter provenance and licenses belong in `data/`; reusable runners live in
`examples/` and `scripts/`.

| System | Tests and main status |
| --- | --- |
| [Water](water/README.md) | Molecular/periodic/charged checks pass; MD and metal-surface NEB qualifications recorded |
| [Zn/O](zno/README.md) | Six routine fixtures and optional 128-atom bulk; numerical checks pass |
| [C/H/O and Cl](cho/README.md) | Molecular/periodic checks pass; charged SN2 path is exploratory |
| [CuO](cuo/README.md) | 96-atom surface; ASE/native/LAMMPS agreement |
| [Ag/ZnO](ag_zno/README.md) | Bulk and supported Ag4 minimisations converge; source identity qualification retained |
| [Anatase](anatase/README.md) | Monti2012 with 1/2/4/8 waters; molecular Ti binding retained |
| [Pd/CeO2](pd_ceria/README.md) | Input only; potential provenance unresolved, no relaxation |
| [Pd/TiO2](pd_tio2/README.md) | Input only; complete potential unavailable, no relaxation |

```sh
python scripts/validate.py --verify
python scripts/validate.py --system water --include-charged --verify
python scripts/validate.py --system cho --include-chocl --verify
python scripts/validate.py --system cuo --verify
python scripts/validate.py --system anatase --verify
python scripts/benchmark.py --compare-neighbors
```

Shared fixtures are defined in `scripts/validate.py`; `--include-bulk` enables
larger systems. Fresh runs retain detailed checks for inspection. LAMMPS
agreement establishes implementation consistency, not physical accuracy.
See [verification conventions](../docs/validation.md).
