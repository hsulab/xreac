# Force fields and cutoffs

## Parameter files

The parser supports the standard LAMMPS ReaxFF format with 39 global parameters
and four lines per atom type. It preserves atom labels and reads bonds,
off-diagonal overrides, repeated angle entries, hydrogen-bond parameters, and
explicit or terminal `0-i-j-0` wildcard torsions. Explicit torsions take
precedence over wildcard defaults. Fortran `D` exponents and an omitted final
hydrogen-bond block are accepted.

```python
from xreac import ForceField

ff = ForceField.bundled("qeq_ff.water")
print(ff.elements)
print(ff.checksum)
# ff = ForceField.from_file("my_force_field.ff")
```

| Bundled file | Example coverage |
| --- | --- |
| `qeq_ff.water` | Achtyl QEq water; H/O and dummy X type |
| `ffield.reax.cho` | Chenoweth C/H/O |
| `ffield.reax.ZnOH` | Raymand 2010 Zn/O/H |

Source files live under the repository's `data/` directory. Wheels include
them as `xreac.data`. Their original contents and citation headers are retained.
Dummy labels such as `X` are parameter metadata, not additional chemical
elements. The ASE adapter requires labels compatible with ASE chemical symbols.

Standard shielded vdW, inner-wall vdW, and their combination are implemented.
Five-line atom extensions, lgvdW, a nonzero lower taper radius, and charge models
such as ACKS2 are not supported. Successfully parsing a file does not establish
its physical suitability or implementation coverage.

## Cutoff values and their sources

| Setting | Source | Bundled/default value |
| --- | --- | --- |
| Nonbonded cutoff | Global parameter 13, `ff.general[12]`, upper taper radius | 10 Å |
| Bond candidate cutoff | Code constant `BOND_CUT` | 5 Å |
| Hydrogen-bond cutoff | Code constant `HBOND_CUT` | 7.5 Å |
| Raw bond-order screening threshold | `0.01 * ff.general[29]` | Parameter dependent |
| Corrected bond order for a donor–H bond | Code constant `HBOND_THRESHOLD` | 0.01 |
| Reported bond-count threshold | Code constant `BOND_GRAPH_CUT` | Strictly greater than 0.3 |

The bond and hydrogen-bond distances are capped by the nonbonded cutoff.
They follow the LAMMPS defaults and currently have no public override option.
The nonbonded cutoff applies to vdW, shielded Coulomb interactions, and QEq
coupling; the nonbonded terms taper to zero there. “Nonbonded” does not exclude
bonded atom pairs from these energy contributions.

Being within 5 Å does not automatically create a bond. Distance-dependent bond
orders and their corrections determine the surviving interactions. Corrected
bond order `0.001` is used for angle/torsion neighbor screening, with additional
product thresholds in those terms. This differs from the `0.3` threshold used
only for reported bond counts.

## Hydrogen bonds

The force-field atom role flags and donor–H–acceptor parameter entries determine
which triplets are eligible. The donor–H corrected bond order must be at least
`0.01`. All acceptor images within the H···acceptor cutoff are considered, and
the energy depends on that distance, donor–H bond order, and triplet angle.

In a small primitive cell, neighbors with different lattice shifts have distinct identities.
A donor oxygen may therefore hydrogen-bond to another image of the same
primitive oxygen. The actual donor is excluded as its own acceptor. This
preserves equivalence with a larger replicated cell; see the
[small-cell validation](validation.md#small-cell-verification).
