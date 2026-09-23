# Small-cell implementation verification

Generated September 23, 2026 in `catorch3` with:

```sh
mamba run -n catorch3 python examples/small_cells.py --verify --output validation/small-cell-support
```

All nine cases pass against LAMMPS 22 Jul 2025, Update 4; **163 tests pass** in
the full regression suite. To reproduce, select a
new output directory or omit `--output` to use a timestamped directory.

| Input cell | Input atoms | Internal atoms | Python energy (kcal/mol/input-cell) |
| --- | ---: | ---: | ---: |
| Water, cubic 9 Å | 3 | 24 | -248.505693095 |
| Water, cubic 6 Å | 3 | 24 | -249.339384818 |
| Water, cubic 4 Å | 3 | 81 | -251.223114721 |
| Water, cubic 3.12 Å | 3 | 192 | -252.309567187 |
| Water dimer, cubic 6 Å | 6 | 48 | -501.221716306 |
| Water, rotated triclinic | 3 | 54 | -250.059751389 |
| Water, partial periodicity | 3 | 27 | -249.969350241 |
| Zn/O, cubic 4 Å | 2 | 54 | -120.377069619 |
| Zinc chain, periodic length 2.5 Å | 1 | 5 | -10.868673983 |

The largest absolute discrepancies over all cases are:

| Quantity | Difference |
| --- | ---: |
| Energy per input atom | 4.3e-12 kcal/mol |
| Energy component per input atom | 9.5e-12 kcal/mol |
| Force component | 4.2e-12 kcal/mol/Å |
| Charge | 6.3e-14 e |
| Dipole component on the input coordinate branch | 1.2e-13 e Å |
| Total bond order | 1.7e-15 |
| Lone-pair count | 9.8e-19 |
| Bond count | Exact match |

Each Python calculation uses the **primitive input**. Positions and charges of
internal copies remain tied, and QEq solves an `(N+1) × (N+1)` constrained
system for `N` input atoms. The energy is divided by the number of copies
before differentiation, which folds all image forces onto the input atoms.
All reference comparisons use fixed-charge forces (`full_derivative=False`).

The 4 Å water cell now includes the -0.370052 kcal/mol hydrogen-bond contribution
missing from a direct primitive-cell LAMMPS run. The one-atom zinc chain has
two bonded neighbors: translated copies of itself. Its aggregated bond-order
matrix has a nonzero diagonal, while its net atomic force is zero.

`python.json` and `reference.json` contain normalized input-cell results;
`structure.json` records the original geometry and PBC. **Raw LAMMPS files
describe the expanded supercell**. Its metadata records the original cell,
input atom count, repetitions, and energy divisor; `primitive.json` preserves
the input. Forces and atom properties are checked for consistency among copies
before averaging. Dipoles are reconstructed from the reference charges on the
input-coordinate branch.

The three-atom 3.12 Å water calculation took about 0.64 s locally, using 192
internal atoms. The 4 Å case took about 0.11 s with 81 internal atoms. These are
single-run timings, not averaged benchmarks. This implementation uses dense
replicated pair arrays; it is not a sparse primitive-cell algorithm. Automatic
replication is capped at 512 internal atoms by default and can be configured
with `max_expanded_atoms`.

These deterministic test geometries are not equilibrated liquids or crystals.
See [summary.json](summary.json) for full errors, tolerances, and timings, and
the [original diagnostic](../small-cells/README.md) for the primitive LAMMPS
hydrogen-bond exclusion that motivated supercell-based verification.
