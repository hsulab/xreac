# Water validation

Inputs cover isolated monomer/dimer, boundary-crossing water, rotated partial
periodicity, small cells with repeated images, and optional 192-atom bulk
water. Charged fixtures include hydroxide, hydronium and charged periodic
cells; `total_charge` is stored in their XYZ metadata. Parameters are
`ffield.reax.HO.2015`. Shared fixture builders remain in `scripts/validate.py`.

All seven stored neutral fixtures passed ASE/native/LAMMPS comparisons;
the recorded maximum force difference was 7.941e-10 kcal/mol/Angstrom.
Charged checks also passed: xreac QEq is checked independently, and LAMMPS
receives supplied charges for nonneutral systems. Monomer relaxation,
finite differences, symmetries, neighbor rebuilding and charge conservation
reuse these inputs. Small-cell hydrogen-bond comparisons use equivalent
supercells because LAMMPS primitive-cell acceptor exclusions differ.

A 1 ps Berendsen MD pilot (0.25 fs, target 300 K) passed geometry-matched
LAMMPS snapshot checks. Recorded xreac/LAMMPS mean times were 237.77/31.93
ms per step on one M1 Pro CPU thread; this is a historical throughput result,
not evidence of equilibrated water or identical integration trajectories.

```sh
python scripts/validate.py --system water --include-charged --verify
python scripts/validate_relaxation.py --backend ase
python examples/water_md.py --steps 4000 --warmup 100
```

## Two-layer Pt(111) and Ni(111) water dissociation

`pt_*_initial.xyz` and `ni_*_initial.xyz` retain the endpoint inputs: p(2x2)
two-layer slabs, eight metal atoms plus H2O, four fixed bottom atoms, xy PBC,
15 Angstrom total vacuum padding and bottom height 2.5 Angstrom. These use
the experimental `ffield.reax.PtNiCHO.2026` derivative, not Assowe2012.
Parameters and the single O-H-Pt edit are documented in [data](../../data/README.md).

Seven-image CI-NEB used charge-response forces; endpoints converged to
0.02 eV/Angstrom and bands to 0.05 eV/Angstrom. Recorded barriers were
0.52115 eV (Ni) and 0.98303 eV (Pt). Ni's final LAMMPS checks passed; Pt's
reactant and peak passed but its product failed equivalent-copy force
consistency. These are local exploratory paths, not validated DFT barriers,
global minimum paths, or slab-converged results.

```sh
python scripts/neb_pt_water.py --metal Ni --ffield data/ffield.reax.PtNiCHO.2026 \
  --layers 2 --vacuum 15 --bottom-height 2.5 --images 7 --neb-dtmax 0.1 --steps 2500 \
  --output validation/runs/ni-water-example
```

Use `--metal Pt` and a fresh output directory for Pt. Only inputs and this
summary are retained; outputs are generated in ignored `validation/runs/`.
Water on anatase has a separate [minimal record](../anatase/README.md).
