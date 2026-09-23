"""Audit single-point QEq using fresh LAMMPS energy finite differences.

Each displacement launches lmp_mpi with zero initial charges and active QEq.
Controls repeat QEq at fixed coordinates and start with different neutral charges.
"""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
from water_cluster import water_cases
from xreac import Calculator, ForceField
from xreac.energy import EnergyModel, QEQ_COULOMB, SELF_CONVERSION, C_ELE
from xreac.neighbors import replicated_neighbors
from xreac.reference import evaluate_lammps


def control(baseline, directory, *, repeat=False):
    directory.mkdir()
    for name in ("atoms.data", "ffield"):
        shutil.copy2(baseline.directory / name, directory / name)
    script = (baseline.directory / "in.lammps").read_text()
    if repeat:
        # No integration fix: coordinates stay fixed while QEq is repeated.
        script = script.replace("run 0\n", "run 0\nrun 5\n")
    else:
        script = script.replace("run 0\n", "set type 1 charge -0.8\nset type 2 charge 0.4\nrun 0\n")
    (directory / "in.lammps").write_text(script)
    metadata = json.loads((baseline.directory / "metadata.json").read_text())
    command = [metadata["executable"], "-in", "in.lammps", "-log", "log.lammps", "-nocite"]
    proc = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=120)
    (directory / "stdout.txt").write_text(proc.stdout)
    (directory / "stderr.txt").write_text(proc.stderr)
    if proc.returncode or "convergence failed" in proc.stdout.lower():
        raise RuntimeError(f"Control failed: {directory}")
    metadata["command"] = command
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    last_frame = (directory / "atoms.dump").read_text().rsplit("ITEM: ATOMS", 1)[1]
    atoms = np.loadtxt(last_frame.splitlines()[1:], ndmin=2)
    energy = float(np.loadtxt(directory / "energy.txt")[0])
    return {
        "energy_difference": energy - baseline.energy,
        "max_charge_difference": float(np.max(abs(atoms[:, 2] - baseline.charges))),
        "max_force_difference": float(np.max(abs(atoms[:, 6:9] - baseline.forces))),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "validation/runs/water/qeq-audit")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    ff = ForceField.bundled("ffield.reax.HO.2015")
    calc = Calculator(ff)
    report = {
        "force_field_sha256": ff.checksum,
        "energy_units": "kcal/mol",
        "force_units": "kcal/mol/Angstrom",
        "charge_units": "e",
        "qeq_residual_units": "eV",
        "qeq_tolerance": 1e-12,
        "cases": {},
    }
    for name in ("monomer",):
        symbols, x = water_cases()[name]
        directory = args.output / name
        actual = calc.evaluate(symbols, x, full_derivative=True)
        fixed = calc.evaluate(symbols, x)
        ref = evaluate_lammps(ff, symbols, x, directory=directory / "baseline")
        index = np.unravel_index(np.argmax(abs(actual.forces - ref.forces)), x.shape)
        neighbors, _ = replicated_neighbors(x, ff.general[12])
        model = EnergyModel(ff, symbols, neighbors)
        _, distances = model.edges.geometry(x)
        _, _, shield = model.electrostatics(distances)
        shield = model.edges.pair_sum(shield)
        hessian = QEQ_COULOMB * shield + np.diag(model.a["eta"])
        chemical_potential = hessian @ ref.charges + model.a["chi"]
        # Equilibration requires equal chemical potentials, not zero individually.
        residual = float(np.max(abs(chemical_potential - chemical_potential.mean())))
        row = {
            "reference_version": ref.version,
            "atom_index_one_based": int(index[0]) + 1,
            "axis": "xyz"[index[1]],
            "charges": ref.charges.tolist(),
            "net_charge": float(ref.charges.sum()),
            "qeq_stationarity_residual": residual,
            "lammps_force": float(ref.forces[index]),
            "python_matched_force": float(fixed.forces[index]),
            "python_full_gradient_force": float(actual.forces[index]),
            "repeated_qeq_control": control(ref, directory / "repeat", repeat=True),
            "different_initial_charges_control": control(ref, directory / "initial_charges"),
            "finite_differences": [],
        }
        for step in (1e-3, 1e-4, 1e-5):
            displaced = []
            for sign, label in ((1, "plus"), (-1, "minus")):
                positions = x.copy()
                positions[index] += sign * step
                displaced.append(evaluate_lammps(ff, symbols, positions, directory=directory / f"h-{step:g}-{label}"))
            plus, minus = displaced
            numerical_force = -(plus.energy - minus.energy) / (2 * step)
            # Charge-response term predicted specifically by the constant mismatch.
            dq = (plus.charges - minus.charges) / (2 * step)
            response = -(C_ELE - SELF_CONVERSION * QEQ_COULOMB) * np.dot(shield @ ref.charges, dq)
            row["finite_differences"].append(
                {
                    "step_A": step,
                    "energy_plus": plus.energy,
                    "energy_minus": minus.energy,
                    "max_charge_change_between_displacements": float(np.max(abs(plus.charges - minus.charges))),
                    "force_from_lammps_energies": float(numerical_force),
                    "minus_lammps_force": float(numerical_force - ref.forces[index]),
                    "minus_python_full_gradient": float(numerical_force - actual.forces[index]),
                    "predicted_charge_response_from_constants": float(response),
                }
            )
        final = row["finite_differences"][-1]
        assert residual < 1e-9
        assert abs(final["minus_python_full_gradient"]) < 1e-5
        assert abs(final["minus_lammps_force"] - final["predicted_charge_response_from_constants"]) < 1e-5
        for key in ("repeated_qeq_control", "different_initial_charges_control"):
            assert abs(row[key]["energy_difference"]) < 1e-9
            assert row[key]["max_charge_difference"] < 1e-10
            assert row[key]["max_force_difference"] < 1e-8
        report["cases"][name] = row
        print(
            f"{name}: LAMMPS force={ref.forces[index]:.10f}; LAMMPS energy derivative={final['force_from_lammps_energies']:.10f}; QEq residual={residual:.3e}",
            flush=True,
        )
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
