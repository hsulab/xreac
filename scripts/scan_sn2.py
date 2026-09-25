"""Quick isolated Cl- + CH3Cl SN2 potential-energy scan with the Hur force field.

Use a collinear Cl-C-Cl axis and a threefold symmetric methyl group. At fixed
xi = r(C-Cl_1) - r(C-Cl_2), relax the mean C-Cl distance and hydrogen radius
and height. Also scan a continuous inversion with linearly constrained H
height and relaxed mean C-Cl distance and H radius. These are restricted
potential-energy scans, not free-energy profiles or transition-state/IRC
verification. Requires scipy and matplotlib.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
for path in ("src", "examples", "tests"):
    sys.path.insert(0, str(ROOT / path))

from validate import CHOCL_FORCE_FIELD, chocl_validation_cases
from xreac import Calculator, ForceField
from xreac.reference import evaluate_lammps
from water_cluster import comparison

PHI = np.arange(3) * 2 * np.pi / 3
RADIAL = np.column_stack((np.cos(PHI), np.sin(PHI)))
SYMBOLS = chocl_validation_cases()["chocl_sn2_reactant"][1]


def geometry(xi, parameters):
    mean, radius, height = parameters
    x = np.zeros((6, 3))
    x[1:4, :2] = radius * RADIAL
    x[1:4, 2] = height
    x[4, 2] = mean + xi / 2
    x[5, 2] = -mean + xi / 2
    return x


def objective(calc, xi, parameters):
    result = calc.evaluate(SYMBOLS, geometry(xi, parameters), total_charge=-1, full_derivative=True)
    f = result.forces
    gradient = -np.array([f[4, 2] - f[5, 2], np.sum(f[1:4, :2] * RADIAL), f[1:4, 2].sum()])
    return result.energy, gradient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=27)
    parser.add_argument("--extent", type=float, default=1.3)
    parser.add_argument("--output", type=Path, default=ROOT / "validation/runs/sn2-scan")
    args = parser.parse_args()
    if args.points < 3 or args.points % 2 != 1 or not 0 < args.extent < 7:
        parser.error("points must be odd and >=3; extent must be in (0, 7) Angstrom")
    args.output.mkdir(parents=True, exist_ok=False)
    ff = ForceField.bundled(CHOCL_FORCE_FIELD)
    calc = Calculator(ff)
    rows = []
    previous = None
    for xi in np.linspace(-args.extent, args.extent, args.points):
        bounds = [(abs(xi) / 2 + 1.2, 5.0), (0.6, 1.4), (-1.0, 1.0)]
        starts = [np.array([max(2.35, 1.85 + abs(xi) / 2), 1.04, z]) for z in (-0.36, 0.0, 0.36)]
        if previous is not None:
            starts.append(previous)
        candidates = [
            minimize(
                lambda p: objective(calc, xi, p),
                np.clip(p, *np.array(bounds).T),
                jac=True,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 300, "ftol": 1e-13, "gtol": 1e-6, "maxls": 40},
            )
            for p in starts
        ]
        best = min(candidates, key=lambda result: result.fun)
        previous = best.x.copy()
        x = geometry(xi, best.x)
        result = calc.evaluate(SYMBOLS, x, total_charge=-1, full_derivative=True)
        energy, gradient = objective(calc, xi, best.x)
        row = dict(
            xi=float(xi),
            energy=energy,
            parameters=best.x.tolist(),
            positions=x.tolist(),
            charges=result.charges.tolist(),
            components=result.components,
            c_cl_distances=[float(best.x[0] + xi / 2), float(best.x[0] - xi / 2)],
            success=bool(best.success),
            message=str(best.message),
            max_parameter_gradient=float(np.max(abs(gradient))),
            active_bounds=bool(
                np.any(np.isclose(best.x, np.array(bounds)[:, 0], atol=1e-5))
                or np.any(np.isclose(best.x, np.array(bounds)[:, 1], atol=1e-5))
            ),
        )
        rows.append(row)
        print(f"xi={xi:+.3f} E={energy:.6f} qCl={result.charges[-2:]} converged={best.success}", flush=True)

    # Paper-inspired rigid interpolation through the shared reactant and symmetric fixtures.
    rigid = []
    cases = chocl_validation_cases()
    reactant = cases["chocl_sn2_reactant"][2]
    symmetric = cases["chocl_sn2_symmetric"][2]
    for xi in np.linspace(-1.3, 1.3, 27):
        endpoint = reactant.copy()
        if xi > 0:
            endpoint[:, 2] *= -1
            endpoint[[4, 5]] = endpoint[[5, 4]]
        x = symmetric + abs(xi) / 1.3 * (endpoint - symmetric)
        rigid.append(dict(xi=float(xi), energy=calc.evaluate(SYMBOLS, x, total_charge=-1).energy))

    energies = np.array([row["energy"] for row in rows])
    left = int(np.argmin(energies[: len(rows) // 2 + 1]))
    right = len(rows) // 2 + int(np.argmin(energies[len(rows) // 2 :]))
    peak = left + int(np.argmax(energies[left : right + 1]))
    # Keep methyl inversion continuous: interpolate its height between the two
    # minima, optimizing only the mean C-Cl distance and H radius at each point.
    continuous = []
    for xi in np.linspace(rows[left]["xi"], rows[right]["xi"], args.points):
        height = np.interp(
            xi,
            [rows[left]["xi"], rows[right]["xi"]],
            [rows[left]["parameters"][2], rows[right]["parameters"][2]],
        )

        def restricted(p):
            energy, gradient = objective(calc, xi, [*p, height])
            return energy, gradient[:2]

        fits = [
            minimize(
                restricted,
                [mean, 1.08],
                jac=True,
                method="L-BFGS-B",
                bounds=[(abs(xi) / 2 + 1.2, 5), (0.6, 1.4)],
                options={"maxiter": 300, "ftol": 1e-13, "gtol": 1e-6, "maxls": 40},
            )
            for mean in (2.1, 2.35, 2.6)
        ]
        fit = min(fits, key=lambda fit: fit.fun)
        continuous.append(
            dict(
                xi=float(xi),
                energy=float(fit.fun),
                parameters=[*fit.x.tolist(), float(height)],
                positions=geometry(xi, [*fit.x, height]).tolist(),
                max_parameter_gradient=float(np.max(abs(fit.jac))),
                success=bool(fit.success),
            )
        )
    midpoint = geometry(0, [*continuous[len(continuous) // 2]["parameters"][:2], 0])
    # A Cartesian curvature check distinguishes a first-order saddle from an
    # unstable symmetry-constrained midpoint. Units: kcal/mol/Angstrom^2.
    curvature_checks = []
    for step in (1e-3, 1e-4, 1e-5):
        hessian = np.zeros((midpoint.size, midpoint.size))
        for i in range(midpoint.size):
            plus, minus = midpoint.copy(), midpoint.copy()
            plus.flat[i] += step
            minus.flat[i] -= step
            hessian[:, i] = -(
                calc.evaluate(SYMBOLS, plus, total_charge=-1, full_derivative=True).forces
                - calc.evaluate(SYMBOLS, minus, total_charge=-1, full_derivative=True).forces
            ).ravel() / (2 * step)
        curvature = np.linalg.eigvalsh((hessian + hessian.T) / 2)
        curvature_checks.append(
            dict(step=step, eigenvalues=curvature.tolist(), negative=int(np.sum(curvature < -1e-2)))
        )
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        force_field=CHOCL_FORCE_FIELD,
        force_field_sha256=ff.checksum,
        total_charge=-1,
        periodic=False,
        energy_units="kcal/mol",
        distance_units="Angstrom",
        symbols=SYMBOLS,
        method=__doc__,
        derivative="Full energy derivative through fresh QEq for constrained optimization",
        bounds="mean C-Cl in [abs(xi)/2+1.2,5]; H radius [0.6,1.4]; H height [-1,1] Angstrom",
        rows=rows,
        rigid=rigid,
        continuous=continuous,
        continuous_barrier=float(max(r["energy"] for r in continuous) - continuous[0]["energy"]),
        continuous_converged=all(r["success"] and r["max_parameter_gradient"] < 1e-4 for r in continuous),
        midpoint_curvature_checks=curvature_checks,
        left_minimum_index=left,
        right_minimum_index=right,
        peak_index=peak,
        barrier_from_left_minimum=float(energies[peak] - energies[left]),
        symmetry_error=float(np.max(abs(energies - energies[::-1]))),
        all_converged=all(
            row["success"] and row["max_parameter_gradient"] < 1e-4 and not row["active_bounds"] for row in rows
        ),
    )
    references = {}
    actual = calc.evaluate(SYMBOLS, midpoint, total_charge=-1)
    ref = evaluate_lammps(
        ff, SYMBOLS, midpoint, supplied_charges=actual.charges, directory=args.output / "reference" / "planar-midpoint"
    )
    references["planar-midpoint"] = comparison(actual, ref, len(midpoint))
    for index in sorted({0, left, len(rows) // 2, peak, right, len(rows) - 1}):
        x = np.asarray(rows[index]["positions"])
        actual = calc.evaluate(SYMBOLS, x, total_charge=-1)
        ref = evaluate_lammps(
            ff, SYMBOLS, x, supplied_charges=actual.charges, directory=args.output / "reference" / f"point-{index:03}"
        )
        references[str(index)] = comparison(actual, ref, len(x))
    report["lammps_supplied_charge_checks"] = references
    (args.output / "scan.json").write_text(json.dumps(report, indent=2) + "\n")
    np.savetxt(
        args.output / "scan.csv",
        np.array([[r["xi"], r["energy"], *r["parameters"], *r["charges"]] for r in rows]),
        delimiter=",",
        header="xi_A,energy_kcal_mol,mean_CCl_A,H_radius_A,H_height_A,qC,qH1,qH2,qH3,qCl1,qCl2",
    )
    xyz = []
    for row in rows:
        xyz.extend(["6", f"xi={row['xi']:.6f} E={row['energy']:.9f} kcal/mol Q=-1"])
        xyz.extend(s + " " + " ".join(f"{v:.12f}" for v in pos) for s, pos in zip(SYMBOLS, row["positions"]))
    (args.output / "scan.xyz").write_text("\n".join(xyz) + "\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True, layout="constrained")
    xi = np.array([r["xi"] for r in rows])
    axes[0].plot(xi, (energies - energies[left]) * 4.184, ":", label="Independent minima (inversion discontinuity)")
    axes[0].plot(
        [r["xi"] for r in continuous],
        [(r["energy"] - continuous[0]["energy"]) * 4.184 for r in continuous],
        label="Continuous methyl inversion; mean C–Cl and H radius relaxed",
    )
    axes[0].plot(
        [r["xi"] for r in rigid],
        [(r["energy"] - rigid[0]["energy"]) * 4.184 for r in rigid],
        "--",
        label="Rigid interpolation (own endpoint zero)",
    )
    axes[0].set_ylabel("Relative potential energy (kJ/mol)")
    axes[0].legend(fontsize=8)
    axes[0].set_title("Cl⁻ + CH₃Cl: Hur ReaxFF, isolated system, Q = −1")
    for i, label in ((4, "Cl₁"), (5, "Cl₂")):
        axes[1].plot(xi, [r["charges"][i] for r in rows], label=label)
    axes[1].set_ylabel("QEq charge (e)")
    axes[1].set_xlabel("r(C–Cl₁) − r(C–Cl₂) (Å)")
    axes[1].legend()
    fig.savefig(args.output / "scan.png", dpi=180)
    fig.savefig(args.output / "scan.pdf")
    print(
        json.dumps(
            {
                key: report[key]
                for key in (
                    "barrier_from_left_minimum",
                    "continuous_barrier",
                    "continuous_converged",
                    "symmetry_error",
                    "all_converged",
                )
            }
        ),
        flush=True,
    )
    if not report["continuous_converged"] or not all(r["passed"] for r in references.values()):
        raise SystemExit("Scan needs review: optimization or reference tolerance not met")


if __name__ == "__main__":
    main()
