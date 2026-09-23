"""Audit primitive LAMMPS cells against equivalent supercells.

Primitive-cell LAMMPS runs below the documented QEq size limit are diagnostics.
The trusted comparison uses xreac and LAMMPS on an exactly replicated cell with
all face heights above the cutoff. Every run is a fixed-geometry single point.
"""

import argparse
from datetime import datetime, timezone
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
from water_cluster import comparison, serialize, water_cases
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps


def image_qeq(ff, symbols, x, cell):
    """Independent dense QEq sum, including nonzero images of the same atom.

    This diagnostic deliberately enumerates a generous, geometry-dependent
    integer shell; it does not reuse xreac's large-cell image generator/masks.
    Only the zero-shift i=j term is excluded from the periodic coupling.
    """
    n = len(x)
    cutoff = ff.general[12]
    inverse = np.linalg.inv(cell)
    x = (x @ inverse % 1) @ cell
    heights = 1 / np.linalg.norm(inverse, axis=0)
    extent = np.ceil(cutoff / heights).astype(int) + 1
    shifts = np.array(list(product(*(range(-v, v + 1) for v in extent))))
    delta = x[None, :, None, :] - x[None, None, :, :] + (shifts @ cell)[:, None, None, :]
    r = np.linalg.norm(delta, axis=-1)
    mask = r < cutoff
    mask &= ~((shifts == 0).all(axis=1)[:, None, None] & np.eye(n, dtype=bool)[None, :, :])
    u = np.minimum(r / cutoff, 1.0)
    taper = (1 + u**4 * (-35 + u * (84 + u * (-70 + 20 * u)))) * mask
    gamma = np.array([[ff.pairs[s, t]["gamma"] for t in symbols] for s in symbols])
    coupling = np.sum(taper / (r**3 + gamma) ** (1 / 3), axis=0)
    chi = np.array([ff.atoms[s]["chi"] for s in symbols])
    eta = np.array([ff.atoms[s]["eta"] for s in symbols])
    hessian = 14.4 * coupling + np.diag(eta)
    kkt = np.block([[hessian, np.ones((n, 1))], [np.ones((1, n)), np.zeros((1, 1))]])
    charges = np.linalg.solve(kkt, np.r_[-chi, 0])[:n]
    bond_candidates = mask & (r <= min(5.0, cutoff))
    counts = mask.sum(axis=0)
    raw = np.zeros_like(r)
    for i, s in enumerate(symbols):
        for j, t in enumerate(symbols):
            pair = ff.pairs[s, t]
            for atom_radius, pair_radius, power1, power2, factor in (
                ("r_s", "r_s", "p_bo1", "p_bo2", 1 + 0.01 * ff.general[29]),
                ("r_pi", "r_p", "p_bo3", "p_bo4", 1.0),
                ("r_pi_pi", "r_pp", "p_bo5", "p_bo6", 1.0),
            ):
                if ff.atoms[s][atom_radius] > 0 and ff.atoms[t][atom_radius] > 0:
                    raw[:, i, j] += factor * np.exp(pair[power1] * (r[:, i, j] / pair[pair_radius]) ** pair[power2])
    active_bonds = bond_candidates & (raw >= 0.01 * ff.general[29])
    return {
        "charges": charges.tolist(),
        "coupling": coupling.tolist(),
        "qeq_self_image_diagonal": np.diag(14.4 * coupling).tolist(),
        "self_images_within_nonbonded_cutoff": np.diag(counts).tolist(),
        "maximum_images_per_atom_pair": int(counts.max()),
        "maximum_bond_candidate_images_per_atom_pair": int(bond_candidates.sum(axis=0).max()),
        "self_images_within_bond_cutoff": np.diag(bond_candidates.sum(axis=0)).tolist(),
        "maximum_active_bond_images_per_atom_pair": int(active_bonds.sum(axis=0).max()),
        "maximum_integer_shift_within_nonbonded_cutoff": int(np.max(abs(shifts[mask.any(axis=(1, 2))]), initial=0)),
    }


def repeat_structure(symbols, x, cell, repeats):
    shifts = np.array(list(product(*(range(n) for n in repeats)))) @ cell
    return symbols * len(shifts), np.concatenate([x + shift for shift in shifts]), cell * np.array(repeats)[:, None]


def differences(primitive, expanded, copies, n):
    """Compare equal geometries per primitive cell; do not compare branch dipoles."""
    result = {
        "energy_difference_per_primitive_cell": primitive.energy - expanded.energy / copies,
        "component_differences_per_primitive_cell": {
            key: primitive.components[key] - expanded.components[key] / copies for key in primitive.components
        },
    }
    for key in ("forces", "charges", "total_bond_orders", "lone_pairs", "bond_counts"):
        expected = np.tile(getattr(primitive, key), (copies, 1) if key == "forces" else copies)
        result["max_" + key + "_difference"] = float(np.max(abs(expected - getattr(expanded, key))))
    result["consistent"] = (
        abs(result["energy_difference_per_primitive_cell"]) / n <= 1e-5
        and max(abs(v) for v in result["component_differences_per_primitive_cell"].values()) / n <= 1e-5
        and result["max_forces_difference"] <= 1e-4
        and result["max_charges_difference"] <= 1e-6
        and result["max_total_bond_orders_difference"] <= 1e-8
        and result["max_lone_pairs_difference"] <= 1e-8
        and result["max_bond_counts_difference"] == 0
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "validation"
        / "runs"
        / ("small-cells-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")),
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    water = ForceField.bundled("qeq_ff.water")
    zno = ForceField.bundled("ffield.reax.ZnOH")
    symbols, x = water_cases()["monomer"]
    cases = [(f"water_{side:g}A", water, symbols, x, np.diag([side] * 3)) for side in (12.0, 9.0, 6.0, 4.0, 3.12)]
    # Avoid a bond aligned exactly with a lattice vector, which activates an
    # undefined collinear dihedral in the densest replicated geometry.
    rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
    name, ff, labels, positions, cell = cases[-1]
    cases[-1] = (name, ff, labels, positions @ rotation, cell)
    ds, dx = water_cases()["distorted_dimer"]
    cases.append(("water_dimer_6A", water, ds, dx, np.diag([6.0] * 3)))
    cases.append(("zno_4A", zno, ["Zn", "O"], np.array([[0.0, 0, 0], [1.9, 0.1, 0.2]]), np.diag([4.0] * 3)))
    report = {
        "scope": "Primitive LAMMPS representation diagnostic; supercells are the reference",
        "energy_units": "kcal/mol",
        "force_units": "kcal/mol/Angstrom",
        "charge_units": "e",
        "lammps_qeq_restriction": "https://docs.lammps.org/fix_qeq_reaxff.html#restrictions",
        "python": sys.version,
        "cases": {},
    }
    for name, ff, symbols, x, cell in cases:
        work = args.output / name
        work.mkdir()
        calc = Calculator(ff)
        heights = 1 / np.linalg.norm(np.linalg.inv(cell), axis=0)
        minimum = max(ff.general[12], 2 * min(5.0, ff.general[12]))
        repeats = (np.floor(minimum / heights) + 1).astype(int)
        ss, sx, sc = repeat_structure(symbols, x, cell, repeats)
        copies, n = int(np.prod(repeats)), len(x)
        (work / "structure.json").write_text(
            json.dumps(
                dict(
                    symbols=symbols, positions=x.tolist(), cell=cell.tolist(), pbc=[True] * 3, repeats=repeats.tolist()
                ),
                indent=2,
            )
            + "\n"
        )
        row = {
            "primitive_atoms": n,
            "cell": cell.tolist(),
            "repeats": repeats.tolist(),
            "expanded_atoms": len(sx),
            "force_field_sha256": ff.checksum,
        }
        try:
            calc.evaluate(symbols, x, cell=cell)
            row["python_primitive_status"] = "supported"
        except ValueError as exc:
            row["python_primitive_status"] = str(exc)
        independent = image_qeq(ff, symbols, x, cell)
        row["explicit_image_qeq"] = independent
        actual = calc.evaluate(ss, sx, cell=sc)
        expanded = evaluate_lammps(ff, ss, sx, cell=sc, directory=work / "lammps_supercell")
        (work / "python_supercell.json").write_text(json.dumps(serialize(actual), indent=2) + "\n")
        (work / "reference_supercell.json").write_text(json.dumps(serialize(expanded), indent=2) + "\n")
        row["supercell_python_lammps"] = comparison(actual, expanded, len(sx))
        row["explicit_qeq_supercell_charge_difference"] = float(
            np.max(abs(np.tile(independent["charges"], copies) - expanded.charges))
        )
        row["reference_version"] = expanded.version
        try:
            small = evaluate_lammps(
                ff, symbols, x, cell=cell, allow_small_cell=True, directory=work / "lammps_primitive"
            )
            (work / "reference_primitive.json").write_text(json.dumps(serialize(small), indent=2) + "\n")
            row["primitive_lammps_status"] = "completed"
            row["primitive_vs_supercell"] = differences(small, expanded, copies, n)
            row["explicit_qeq_primitive_charge_difference"] = float(
                np.max(abs(np.asarray(independent["charges"]) - small.charges))
            )
        except RuntimeError as exc:
            row["primitive_lammps_status"] = str(exc)
        report["cases"][name] = row
        print(
            name
            + ": "
            + json.dumps(
                {
                    "supercell_verified": row["supercell_python_lammps"]["passed"],
                    "primitive_vs_supercell": row.get("primitive_vs_supercell", row["primitive_lammps_status"]),
                }
            ),
            flush=True,
        )
        (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    report["supported_supercells_passed"] = all(
        row["supercell_python_lammps"]["passed"] for row in report["cases"].values()
    )
    report["explicit_qeq_supercells_passed"] = all(
        row["explicit_qeq_supercell_charge_difference"] <= 1e-6 for row in report["cases"].values()
    )
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    if not report["supported_supercells_passed"] or not report["explicit_qeq_supercells_passed"]:
        raise SystemExit("A supported supercell or explicit image QEq sum failed verification")


if __name__ == "__main__":
    main()
