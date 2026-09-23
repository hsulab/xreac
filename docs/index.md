# xreac

xreac is a NumPy + Autograd implementation of ReaxFF for neutral molecules,
clusters, and periodic cells. It targets small-to-medium systems, roughly
1–200 atoms, and computes energies, equilibrated charges, forces, and bond
properties from standard parameter files.

Core calculations need only NumPy and HIPS Autograd. The optional ASE adapter
provides familiar atom objects and optimizers. LAMMPS is used for independent
verification and is not required to use the calculator or build these docs.

```python
from xreac import Calculator, ForceField

calc = Calculator(ForceField.bundled("qeq_ff.water"))
symbols = ["O", "H", "H"]
positions = [[0, 0, 0], [0.97, 0, 0], [-0.243, 0.94, 0]]
result = calc.evaluate(symbols, positions)
print(result.energy)   # kcal/mol
print(result.forces)   # kcal/mol/Angstrom; fixed-charge derivative
```

Start with [installation and a working calculation](getting-started.md), then
read about [force conventions](calculations.md) and [periodic cells](periodic.md).

```{toctree}
:maxdepth: 2
:caption: User guide

getting-started
calculations
force-fields
periodic
ase-relaxation
validation
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
```

## Scope

The implementation includes bond, coordination, lone-pair, angle, torsion,
conjugation, hydrogen-bond, vdW, Coulomb, and QEq self-energy terms. Small
periodic cells use tied internal replication, with results normalized to the
input cell. Fixed-cell relaxation uses fixed-charge forces with either ASE FIRE
or the retained native FIRE optimizer.

Net charge, stress, variable-cell relaxation, MD integration, external fields,
and alternative charge models such as ACKS2 are not implemented. Exactly
collinear active torsions are rejected because their dihedral derivative is
undefined. A parsed force field still needs validation for its intended use.

The code is GPL-2.0-or-later. Equations and conventions follow LAMMPS/PuReMD;
bundled parameter files retain their citations. Download the
{download}`license <../LICENSE>` and {download}`attribution notice <../NOTICE>`.
