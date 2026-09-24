# Molecular and periodic examples

`water_cluster.py` runs entirely from the repository with NumPy and Autograd:

```sh
python examples/water_cluster.py
```

To also verify against the pinned `lmp_mpi` executable:

```sh
python examples/water_cluster.py --verify --output validation/runs/my-water-check
```

Choose a new output directory. The script saves structures in XYZ format,
both sets of numerical results, raw LAMMPS inputs/outputs, and a summary with
explicit tolerances. Nonmatching reference results cause a nonzero exit status.
The default data source is the bundled LAMMPS QEq water parameter file; it
does not depend on the ZnOH example or on any Zn-specific calculator logic.

The example uses a monomer and a distorted dimer. The same geometries support
single-point energy, force, charge, property, symmetry, and derivative checks.
These are unoptimized validation fixtures.

Verified properties are charges, dipole vectors (e Å), per-atom total bond
orders, lone-pair counts, and bond counts at the LAMMPS default threshold of
0.3. Force comparisons use the default fixed-charge `forces`
(`full_derivative=False`), matching LAMMPS. The example computes no
charge-response forces. Those are available through
`calc.evaluate(symbols, positions, full_derivative=True)` and are checked
separately by finite differences in the tests.

See [saved results](../validation/water/summary.json) and the
[package README](../README.md) for supported formats and limitations.

`ase_water.py` uses the optional `ReaxFFCalculator` adapter with ASE optimizers:

```sh
python examples/ase_water.py --verify
python examples/ase_water.py --optimizer BFGS --verify
```

It writes the initial/final XYZ structures, an ASE trajectory, optimizer log,
and results in a fresh directory. ASE energies/forces and optimizer tolerance
use eV and eV/Angstrom; the saved core `python.json` and reference comparison
use kcal/mol and kcal/mol/Angstrom. Both optimizers use fixed-charge forces
and equilibrate charges at every evaluated geometry.

`Calculator.relax()` uses ASE FIRE by default. The original FIRE implementation
remains available through `backend="native"` and in the native relaxation benchmark.

`periodic_water.py` verifies fixed-cell periodic calculations:

```sh
python examples/periodic_water.py --verify
python examples/periodic_water.py --verify --include-bulk
```

It covers a boundary-crossing dimer and a rotated triclinic water slab.
The optional `--include-bulk` adds the 192-atom box. Structures are deterministic
fixtures, not equilibrated liquid snapshots. See the
[water results](../validation/water/README.md).

`water_md.py` runs the existing 192-atom water box with ASE's fixed-volume
Berendsen thermostat and compares moving-MD throughput with single-rank LAMMPS:

```sh
python examples/water_md.py --steps 1000 --warmup 100
# Longer example: 1 ps of timed dynamics at 0.25 fs/step.
python examples/water_md.py --steps 4000 --warmup 100
```

Defaults are 300 K and a 100 fs thermostat time constant. Both engines use the
same initial positions, seeded velocities, force-field masses, timestep, and
3N temperature degrees of freedom. ASE uses its neighbor list with a 0.3 Å
per-atom skin; LAMMPS uses a 0.6 Å pair skin. QEq is converged at every step.
The timer excludes warmup and reference checks. Three frames from each run
are checked against fresh LAMMPS calculations; mismatches stop the example.
Outputs include sampled coordinates, numerical checks, timings, and raw LAMMPS
inputs/outputs under `validation/runs/`.

ASE scales velocities before velocity Verlet; LAMMPS's
[`temp/berendsen`](https://docs.lammps.org/fix_temp_berendsen.html) scales after
integration with `fix nve`. Therefore trajectories are not required to match
point by point; energies, forces, charges, and properties are compared at the
same sampled geometries. The unrelaxed fixture heats during initial settling;
these runs demonstrate implementation and throughput, not equilibrated liquid
statistics. See [ASE's MD documentation](https://ase.gitlab.io/ase/ase/md.html)
for thermostat usage and ensemble limitations.

`small_cells.py` demonstrates repeated images and verification against normalized
larger LAMMPS supercells:

```sh
python examples/small_cells.py --verify
```

Four cases cover the 4 Å water monomer, a small triclinic water slab, periodic
ZnO, and a zinc atom bonded to its own images. Both neighbor builders feed one
energy model; QEq solves for input-cell charges. The native builder uses
replicated search cells; ASE builds image-resolved neighbors directly.

The shared system validation runner combines these structures with the C/H/O
and Zn/O checks, retaining each structure only once:

```sh
python scripts/validate.py --verify
python scripts/validate.py --verify --system cho
```

See [validation](../validation/README.md) for the coverage and retained records.

## Three performance pilots

```sh
python scripts/benchmark_lammps.py
```

The default benchmark runs exactly three shared cases: 192-atom bulk water,
128-atom wurtzite ZnO bulk, and the 96-atom CuO(010) surface. It uses one CPU
thread, builds **fresh ASE neighbor lists inside every timed evaluation**, and
checks all three against LAMMPS. ASE result caching is bypassed explicitly.
`scripts/benchmark.py` is an alias for the same driver and options.
Use `--skip-lammps-timing` to skip the LAMMPS timing loops while retaining
mandatory numerical verification. These fixtures are unrelaxed test geometries.

`cuo_surface.py` adds one optional 96-atom CuO(010) slab. It always checks
energies, forces, charges, and properties against LAMMPS, and measures
single-CPU performance. It uses the bundled Cu/O/H/Cl supplement, separately
licensed under **CC BY-NC 4.0**; see [parameter licenses](../data/README.md)
and [saved results](../validation/cuo/README.md).

```sh
python examples/cuo_surface.py
```

This example also times ASE neighbors by default. Its untimed native evaluation
checks backend agreement. `--neighbor-backend replicated` is available only
for reproducing the older native-neighbor benchmark records.

## Experimental QEq history reuse

```sh
python scripts/benchmark_water_qeq.py --steps 4000 --warmup 100
```

This experiment reruns the same bulk-water Berendsen example. It retains the
previous two QEq solutions and an LU factorization, predicts the next solution,
and refines against the **newly rebuilt** QEq matrix. It refactors when six
corrections fail to reach an absolute KKT residual of `1e-12`. Charges are
converged at every geometry; no force or charge update is skipped.

Every step is also solved directly on the identical matrix, with alternating
solver order. Charge differences must stay below `1e-10` e. Cached charges
drive the trajectory, and three sampled states are checked against fresh
LAMMPS calculations. The full experimental timer includes both solvers and
verification overhead. The report separates their costs and gives an estimated
fraction of MD time that caching could save; this is not an independently
measured production speedup. The production calculator keeps its direct solver.
