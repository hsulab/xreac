"""Screen a one-parameter regularization of Gai2016's O-H-Pt angle cutoff.

This is a numerical experiment, not a DFT fit. The native ReaxFF equations and
cutoffs are unchanged. Only O-H-Pt p_val4 is varied from its original 1.0000.
"""

import argparse
import json
from pathlib import Path

from ase.optimize import FIRE
from ase.io import write

from screen_pt_water import closest_pt_h_cutoff, describe
from neb_pt_water import ForceField, ReaxFFCalculator, geometry, np


def variant_text(original, exponent):
    if not 1 < exponent <= 1.2:
        raise ValueError("This experimental screen only allows 1 < p_val4 <= 1.2")
    rows = original.splitlines()
    start = next(i for i, row in enumerate(rows) if "Nr of angles" in row)
    stop = start + 1 + int(rows[start].split()[0])
    hits = []
    for i in range(start + 1, stop):
        tokens = rows[i].split()
        if tokens[:3] == ["3", "2", "5"]:
            if len(tokens) != 10 or tokens[-1] != "1.0000":
                raise ValueError("Unexpected source O-H-Pt angle record")
            tokens[-1] = f"{exponent:.4f}"
            rows[i] = " ".join(tokens)
            hits.append(i)
    if len(hits) != 1:
        raise ValueError("Expected exactly one O-H-Pt angle record")
    rows[0] = (
        f"EXPERIMENTAL Gai2016 derivative; O-H-Pt p_val4 1.0000 -> {exponent:.4f}; "
        "not DFT fitted; DOI:10.1021/acs.jpcc.6b01064; CC-BY-NC-4.0; see NOTICE"
    )
    return "\n".join(rows) + "\n"


def minimize(a, ff, output, steps=800):
    a.calc = ReaxFFCalculator(ff, full_derivative=True)
    opt = FIRE(a, dt=0.02, maxstep=0.02, logfile=str(output))
    history = []
    for _ in opt.irun(fmax=0.02, steps=steps):
        history.append(a.positions.copy())
        if len(history) > 30 and np.max(abs(history[-1] - history[-30])) < 1e-10:
            break
    fmax = float(np.linalg.norm(a.get_forces(), axis=1).max())
    return dict(
        converged=fmax <= 0.02,
        fmax_ev_angstrom=fmax,
        energy_ev=a.get_potential_energy(),
        steps=opt.nsteps,
        positions=a.positions.tolist(),
        geometry=describe(a),
        closest_pt_h_cutoff=closest_pt_h_cutoff(a, ff),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    original = ForceField.bundled("ffield.reax.PtNiCHO.2016")
    report = dict(
        source_sha256=original.checksum,
        parameter="O-H-Pt p_val4",
        source_value=1.0,
        fmax_target_ev_angstrom=0.02,
        full_derivative=True,
        cases=[],
    )
    for value in (1.0, 1.01, 1.025, 1.05, 1.1, 1.2):
        folder = args.output / f"pval4-{value}"
        folder.mkdir()
        if value == 1:
            ff = original
        else:
            path = folder / "ffield"
            path.write_text(variant_text(original.path.read_text(), value))
            ff = ForceField.from_file(path)
        row = dict(value=value, sha256=ff.checksum, endpoints={})
        for name, product in (("reactant", False), ("product", True)):
            a = geometry(product=product)
            result = minimize(a, ff, folder / f"{name}.log")
            row["endpoints"][name] = result
            write(folder / f"{name}.traj", a)
            write(folder / f"{name}.extxyz", a)
            print(value, name, result["converged"], result["fmax_ev_angstrom"], result["energy_ev"], flush=True)
        report["cases"].append(row)
        (args.output / "screen.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
