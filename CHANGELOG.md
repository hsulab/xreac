# Changelog

## 0.6.0

First public release of xreac, a NumPy and Autograd ReaxFF implementation for
small-to-medium neutral molecules, clusters, and periodic cells.

- Parameter-driven energies, QEq charges, forces, dipoles, and bond properties.
- ASE calculator and optimizers, with native FIRE retained for benchmarks.
- Fixed-cell periodic boundaries, including small cells and repeated images.
- One energy model accepting explicit `(i, j, S)` neighbor arrays; ASE and native
  neighbor builders produce matching results.
- Bundled `ffield.reax.HO.2015`, `ffield.reax.CHO.2008`, and
  `ffield.reax.ZnOH.2010` parameter files with original citations and contents.
- System-organized validation: 16 shared routine structures and an optional
  192-atom water box. All 163 regression tests and fresh comparisons with
  LAMMPS 22 Jul 2025, Update 4 pass locally.
- Markdown/Sphinx documentation and retained numerical records and water PDF.

Requires Python 3.10 or newer. Install the wheel attached to this release,
or install from the tag:

```sh
python -m pip install 'git+https://github.com/hsulab/xreac.git@v0.6.0'
```

The `ase` extra enables ASE integration; `dev` adds testing, documentation,
and formatting tools. Core dependencies are NumPy and HIPS Autograd.

Forces default to the LAMMPS fixed-charge convention, with fresh QEq at each
geometry. Full differentiation through QEq is an explicit single-point option.
Relaxation always uses fixed-charge forces. Net charge, stress, variable-cell
relaxation, and MD integration are not supported in this release.
