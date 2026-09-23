"""Relax water with ASE FIRE or BFGS using the xreac calculator.

Run: python examples/ase_water.py --verify
"""

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS, FIRE

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from water_cluster import comparison, serialize, water_cases
from xreac import ForceField
from xreac.ase import ReaxFFCalculator
from xreac.reference import evaluate_lammps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optimizer", choices=("FIRE", "BFGS"), default="FIRE")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "validation/runs" / ("ase-water-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")),
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    ff = ForceField.bundled("qeq_ff.water")
    symbols, positions = water_cases()["monomer"]
    atoms = Atoms(symbols, positions=positions, calculator=ReaxFFCalculator(ff))
    write(args.output / "initial.xyz", atoms, format="xyz")
    initial_energy = atoms.get_potential_energy()
    optimizer_class = {"FIRE": FIRE, "BFGS": BFGS}[args.optimizer]
    with optimizer_class(
        atoms, logfile=str(args.output / "optimize.log"), trajectory=str(args.output / "relaxation.traj")
    ) as optimizer:
        converged = optimizer.run(fmax=1e-5, steps=500)
        steps = optimizer.nsteps
    forces = atoms.get_forces()
    result = atoms.calc.evaluation
    write(args.output / "relaxed.xyz", atoms, format="xyz")
    (args.output / "python.json").write_text(json.dumps(serialize(result), indent=2) + "\n")
    report = {
        "optimizer": f"ASE {args.optimizer}",
        "ase_version": version("ase"),
        "force_field_sha256": ff.checksum,
        "full_derivative": False,
        "force_convention": "fixed_charge",
        "energy_units": "eV",
        "force_units": "eV/Angstrom",
        "initial_energy": initial_energy,
        "energy": atoms.get_potential_energy(),
        "charges": atoms.get_charges().tolist(),
        "dipole_e_A": atoms.get_dipole_moment().tolist(),
        "max_atomic_force_norm": float(np.max(np.linalg.norm(forces, axis=1))),
        "fmax": 1e-5,
        "converged": bool(converged),
        "steps": steps,
    }
    if args.verify:
        ref = evaluate_lammps(ff, symbols, atoms.positions, directory=args.output / "lammps")
        (args.output / "reference.json").write_text(json.dumps(serialize(ref), indent=2) + "\n")
        report["reference"] = comparison(result, ref, len(atoms))
        report["reference"]["units"] = "Core units: kcal/mol, kcal/mol/Angstrom, e, e Angstrom"
    report["passed"] = bool(converged and report.get("reference", {"passed": True})["passed"])
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"ASE {args.optimizer}: {'PASS' if report['passed'] else 'FAIL'}, {steps} steps, E={report['energy']:.10f} eV"
    )
    print(f"Results: {args.output}")
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
