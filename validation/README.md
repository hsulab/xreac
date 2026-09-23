# Validation

Results are organized by chemical system. The default suite uses **16 structures**
for energy, forces, QEq charges, dipoles, bond properties, and neighbor-list
agreement. Each structure is stored once. The 192-atom water box is optional.

| System | Default cases | Purpose |
| --- | ---: | --- |
| [Water](water/README.md) | 6 | Molecules, hydrogen bonds, boundary crossings, small cells, triclinic/partial periodicity |
| [Zn/O](zno/README.md) | 6 | Atom, dimer, oxygen many-body terms, cluster, periodic/self-image bonds |
| [C/H/O](cho/README.md) | 4 | Methane, CO, C2 correction, periodic torsions |

```sh
python scripts/validate.py --verify
python scripts/validate.py --verify --system water --include-bulk
python -m pytest -q
```

Use `--output PATH` for a new output directory. Fresh runs go under ignored
`validation/runs/`; only selected records belong in the system folders.
Omit `--verify` for a comparison of the two Python neighbor backends alone.
Reference checks require `lmp_mpi` (LAMMPS 22 Jul 2025, Update 4).

Every system has the same compact files:

- `structures.json`: coordinates, atom labels, cell, and PBC for each case.
- `summary.json`: tolerances, discrepancies, timings, versions, and checksums.
- `results.json.gz`: full `ase`, `replicated`, and `reference` numerical results.
- `reference.tar.gz`: complete LAMMPS inputs and outputs for those cases.
- `baseline.json.gz`: selected independent results from before energy-model
  consolidation, used by regression tests. These values were not regenerated.

Read compressed results with `json.load(gzip.open(path, "rt"))`. The baseline
contains one native result per default case. The optional bulk case has fresh
LAMMPS verification but is outside the routine baseline suite.

Water and Zn/O also retain ASE/native FIRE relaxation on their shared monomer
and dimer. Water contains two targeted diagnostics and a [PDF report](water/report.pdf).
Large Zn/O sizes are confined to optional [performance benchmarks](zno/benchmarks.json).

Cutoff checks perturb the ZnO dimer around 5 and 10 Å. Symmetry, finite-difference,
QEq, cache, and optimizer checks reuse the same geometries. Distinct checks for
invalid inputs, inner-wall variants, and singular torsions remain in the tests.

The previous feature/version-based records and larger overlapping case lists
are available in Git history. Numerical tolerances are unchanged; see the
[verification guide](../docs/validation.md). These unoptimized fixtures verify
implementation consistency, not physical accuracy of a parameterization.
