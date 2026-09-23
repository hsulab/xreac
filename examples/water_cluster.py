"""Calculate water clusters and optionally verify energies, forces, and properties.

Run: python examples/water_cluster.py --verify --output validation/water
The output directory must not already exist. LAMMPS is only needed with --verify.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def water_cases():
    """Deterministic nonperiodic geometries in Angstrom; not optimized structures."""
    monomer = np.array([[0., 0., 0.], [.97, 0., 0.], [-.243, .94, 0.]])
    half = np.deg2rad(104.5/2)
    acceptor = np.array([[0., 0., 0.], [.97*np.cos(half), .97*np.sin(half), 0.],
                         [.97*np.cos(half), -.97*np.sin(half), 0.]]) + [2.85, .15, .12]
    cases = {"monomer": (["O", "H", "H"], monomer),
             "dimer": (["O", "H", "H"]*2, np.vstack((monomer, acceptor)))}
    distorted = cases["dimer"][1].copy()
    distorted[1] += [.09, .04, -.02]
    distorted[4] += [-.04, .03, .06]
    cases["distorted_dimer"] = (cases["dimer"][0], distorted)
    for n, name in ((3, "trimer"), (6, "hexamer")):
        phi = np.arange(n)*2*np.pi/n
        radius = 2.85/(2*np.sin(np.pi/n))
        oxygens = np.column_stack((radius*np.cos(phi), radius*np.sin(phi), .1*np.cos(2*phi)))
        positions = []
        for i in range(n):
            direction = oxygens[(i+1) % n]-oxygens[i]
            direction /= np.linalg.norm(direction)
            transverse = np.cross(direction, [0., 0., 1.])
            transverse /= np.linalg.norm(transverse)
            # Tilt the donor O-H bond away from exact O-H...O collinearity:
            # weak intermolecular bonds otherwise activate singular dihedrals.
            tilt = np.deg2rad(4.)
            direction, transverse = (np.cos(tilt)*direction+np.sin(tilt)*transverse,
                                      -np.sin(tilt)*direction+np.cos(tilt)*transverse)
            angle = np.deg2rad(104.5)
            second = np.cos(angle)*direction+np.sin(angle)*transverse
            positions.extend((oxygens[i], oxygens[i]+.97*direction, oxygens[i]+.97*second))
        cases[name] = (["O", "H", "H"]*n, np.array(positions))
    return cases


def serialize(result):
    fields = ("energy", "forces", "charges", "components", "dipole",
              "total_bond_orders", "lone_pairs", "bond_counts")
    output = {name: getattr(result, name) for name in fields}
    for name in ("lammps_forces", "bond_orders"):
        if hasattr(result, name):
            output[name] = getattr(result, name)
    return {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in output.items()}


def comparison(actual, reference, n):
    errors = {
        "energy_per_atom": abs(actual.energy-reference.energy)/n,
        "components_per_atom": max(abs(actual.components[k]-reference.components[k]) for k in actual.components)/n,
        "lammps_forces": float(np.max(abs(actual.lammps_forces-reference.forces))),
    }
    for name in ("charges", "dipole", "total_bond_orders", "lone_pairs", "bond_counts"):
        errors[name] = float(np.max(abs(getattr(actual, name)-getattr(reference, name))))
    limits = {"energy_per_atom": 1e-5, "components_per_atom": 1e-5, "lammps_forces": 1e-4,
              "charges": 1e-6, "dipole": 1e-6, "total_bond_orders": 1e-8,
              "lone_pairs": 1e-8, "bond_counts": 0}
    return {"maximum_absolute_errors": errors, "limits": limits,
            "passed": all(errors[k] <= limits[k] for k in errors)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Run lmp_mpi and fail on reference discrepancies")
    parser.add_argument("--ffield", type=Path, help="Use another supported parameter file")
    parser.add_argument("--executable", help="Override the lmp_mpi executable")
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "runs" / ("water-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    args = parser.parse_args()
    ff = ForceField.from_file(args.ffield) if args.ffield else ForceField.bundled("qeq_ff.water")
    calc = Calculator(ff)
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {"force_field": ff.path.name, "sha256": ff.checksum, "citation": ff.citation,
               "units": {"energy": "kcal/mol", "forces": "kcal/mol/Angstrom", "charges": "e", "dipole": "e Angstrom"},
               "reference_verified": args.verify,
               "force_convention": "Full QEq response in forces; fixed-charge derivative in lammps_forces",
               "cases": {}}
    for name, (symbols, positions) in water_cases().items():
        directory = args.output / name
        directory.mkdir()
        xyz = [str(len(symbols)), f"xreac water {name}; coordinates in Angstrom"]
        xyz += [s+" "+" ".join(f"{v:.17g}" for v in x) for s, x in zip(symbols, positions)]
        (directory / "structure.xyz").write_text("\n".join(xyz)+"\n")
        result = calc.evaluate(symbols, positions)
        (directory / "python.json").write_text(json.dumps(serialize(result), indent=2)+"\n")
        row = {"atoms": len(symbols), "energy": result.energy, "hydrogen_bond_energy": result.components["hydrogen_bond"],
               "dipole": result.dipole.tolist(), "total_charge": float(result.charges.sum())}
        if args.verify:
            ref = evaluate_lammps(ff, symbols, positions, directory=directory/"lammps", executable=args.executable)
            (directory/"reference.json").write_text(json.dumps(serialize(ref), indent=2)+"\n")
            row.update(comparison(result, ref, len(symbols)))
            row["reference_version"] = ref.version
        summary["cases"][name] = row
        print(f"{name}: E={result.energy:.9f} kcal/mol" + (f", reference {'PASS' if row['passed'] else 'FAIL'}" if args.verify else ""), flush=True)
    if args.verify:
        summary["passed"] = all(v["passed"] for v in summary["cases"].values())
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(f"Results: {args.output}")
    if args.verify and not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
