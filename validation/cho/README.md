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

## Hur C/H/O/Cl parameters

[Hur 2021 results](hur2021.json) retain CH3Cl, HCl, formyl chloride, and two
charged Cl- + CH3Cl geometries using the bundled `ffield.reax.CHOCl.2021`.
All five pass energy, component, fixed-charge force, and property comparisons
with LAMMPS. Neutral cases use fresh LAMMPS QEq; charged cases supply xreac
charges and validate QEq independently. The idealized SN2 geometries use
distances from Bucko (2008), but are not optimized stationary points and do
not constitute a barrier prediction.

[Raw reference runs](hur2021_reference.tar.gz) include force fields, structures,
scripts, outputs, and metadata. [Source parameter pages and provenance](hur2021_source.tar.gz)
record the PDF extraction. The data is CC BY-NC 3.0; see `data/README.md`.

```sh
python scripts/validate.py --system cho --include-chocl --verify
```

## Quick charged SN2 scan

[Scan data and geometries](sn2_scan.json), [plot](sn2_scan.png), and
[raw scan/reference archive](sn2_scan_raw.tar.gz) record an isolated
Cl− + CH3Cl exchange with total charge −1 and the Hur parameters. Fresh
xreac QEq is solved at every geometry; selected LAMMPS checks use those
supplied charges. Reproduce with:

```sh
python scripts/scan_sn2.py --output validation/runs/sn2-scan
```

The 27-point continuous inversion scan covers
xi = r(C–Cl1) − r(C–Cl2) from −1.3 to +1.3 Å. It constrains a collinear
Cl–C–Cl axis, threefold methyl symmetry, and hydrogen height interpolated
between the endpoint structures; mean C–Cl distance and H radius relax.
All these constrained optimizations meet the 1e-4 kcal/mol/Å gradient
tolerance. The midpoint lies **29.8379 kcal/mol (124.842 kJ/mol)** above
the constrained reactant complex. The C–Cl distances change from
1.5692/2.8692 Å at the endpoint to 2.4410/2.4410 Å at the midpoint.
Endpoints are grid-selected constrained structures, not fully optimized
reactant minima or separated reactants.

Two diagnostic scans illustrate the sensitivity to path constraints:

- A rigid interpolation of the paper-inspired structures rises by
  20.6937 kcal/mol (86.582 kJ/mol), relative to its own unrelaxed endpoint.
- Independently relaxing methyl height at each xi gives a lower envelope
  rising by 18.0158 kcal/mol (75.378 kJ/mol), but jumps between opposite
  pyramidal methyl orientations near xi = 0. This is not a continuous
  inversion pathway. Some envelope points fail the strict gradient check;
  their diagnostics are retained rather than treated as converged minima.

**An apparent barrier exists along the constrained continuous path, but a
minimum-energy barrier is not established.** Cartesian finite differences
at the planar midpoint give three negative curvatures. Two vary strongly
with displacement (1e-3, 1e-4, and 1e-5 Å), so this calculation does not
establish a smooth first-order saddle or trustworthy vibrational modes.
An unconstrained path/saddle search and investigation of this numerical
behavior would be needed before interpreting the height as an activation
energy. These are potential energies, not free energies.

All four retained LAMMPS comparisons pass: maximum energy discrepancy
5.46e-9 kcal/mol/atom and fixed-charge force discrepancy
4.21e-8 kcal/mol/Å. This verifies energy/force implementation at those
geometries, not the physical accuracy of the parameterization or charged QEq.

## Quick ASE NEB follow-up

[NEB results](sn2_neb.json), [energy plot](sn2_neb.png), and
[archived runs](sn2_neb_raw.tar.gz) record a nine-image unconstrained
[ASE climbing-image NEB](https://docs.ase-lib.org/ase/neb.html) calculation.
The full derivative through fresh QEq is used, with total charge −1.
The reactant was relaxed to a maximum force below 0.01 eV/Å; the product
is its symmetry-related chlorine exchange. Initial images came from the
continuous scan, corrected for endpoint relaxation and perturbed by seeded
Cartesian noise with amplitude 0.015 Å to break exact symmetry. NEB removes
overall rotation/translation but imposes no internal-coordinate constraints.
FIRE uses the improved tangent, spring constant 0.1 eV/Å², and maximum step
0.04 Å. The run converged after 80 regular and 22 climbing-image steps.

The resulting barrier is **125.255 kJ/mol (29.9376 kcal/mol)** with maximum
NEB force **0.00878 eV/Å**, below the 0.01 eV/Å stopping tolerance. The peak
C–Cl distances are 2.4397 and 2.4424 Å. Sampling 20 subdivisions of every
inter-image segment found no higher energy than the climbing image.
This is not a meaningful reduction from the constrained scan's 124.842
kJ/mol; endpoint relaxation slightly changes the energy zero. An initial
looser run (0.05 eV/Å tolerance) gave 125.265 kJ/mol.

A second initialization with 0.1 Å noise did not converge in 200 regular
steps. Its bent images all fell below the endpoints, but dense interpolation
between them exposed a **523.062 kJ/mol** peak: the sparse band skipped the
barrier. Continuing climbing-image optimization on that unresolved band
was unstable and is retained only as a failed diagnostic. The script now
skips climbing if the regular band has not converged. A separate 0.03 Å
endpoint perturbation relaxed to within 0.017 kJ/mol of the original endpoint.

Thus this quick NEB test found **no verified smaller barrier**. NEB force
convergence on the near-collinear branch does not resolve the midpoint
curvature concerns above or exclude other pathways. The converged peak
passes a supplied-charge LAMMPS comparison (energy error below 1e-12
kcal/mol/atom and force error below 1e-12 kcal/mol/Å).

```sh
python scripts/neb_sn2.py --output validation/runs/sn2-neb-tight --fmax 0.01
```
