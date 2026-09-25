"""Render a PDF and CSV from saved water validation results; no calculations rerun.

Install the optional ASE extra, then run python scripts/water_report.py.
"""

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import tempfile

# Keep plotting caches out of the user's home directory.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "xreac-mpl"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INK, BLUE, RED = "#193247", "#167c91", "#d54c46"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    }
)


def read_case(name, structure_data, result_data):
    symbols = structure_data["symbols"]
    xyz = np.asarray(structure_data["positions"])
    actual, reference = result_data["ase"], result_data["reference"]
    if actual.get("full_derivative") is not False:
        raise ValueError("Reference report requires fixed-charge forces")
    delta = np.asarray(actual["forces"]) - reference["forces"]
    assert delta.shape == xyz.shape
    values = {
        "case": name,
        "atoms": len(symbols),
        "python_energy_kcal_mol": actual["energy"],
        "lammps_energy_kcal_mol": reference["energy"],
        "delta_energy_kcal_mol": actual["energy"] - reference["energy"],
        "max_force_component_error_kcal_mol_A": float(np.max(abs(delta))),
        "rms_force_component_error_kcal_mol_A": float(np.sqrt(np.mean(delta**2))),
    }
    assert all(np.isfinite(v) for k, v in values.items() if k != "case" and v is not None)
    return dict(
        name=name,
        symbols=symbols,
        xyz=xyz,
        actual=actual,
        reference=reference,
        delta=delta,
        values=values,
    )


def page(title, subtitle, number):
    fig = plt.figure(figsize=(11.7, 8.3), facecolor="white")
    fig.text(0.055, 0.944, title, fontsize=23, weight="bold")
    fig.text(0.055, 0.906, subtitle, fontsize=10, color=BLUE)
    fig.text(0.055, 0.026, "XREAC  /  Water-cluster reference validation", fontsize=8, color="#617383")
    fig.text(0.945, 0.026, str(number), ha="right", fontsize=8, color="#617383")
    return fig


def table(fig, bounds, headings, rows, widths=None, size=9):
    ax = fig.add_axes(bounds)
    ax.axis("off")
    tab = ax.table(
        cellText=rows, colLabels=headings, colWidths=widths, cellLoc="right", colLoc="right", bbox=[0, 0, 1, 1]
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(size)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_facecolor(INK if r == 0 else ("#eef4f7" if r % 2 else "#f8fafb"))
        cell.set_text_props(color="white" if r == 0 else INK, weight="bold" if r == 0 else "normal")
        if c == 0:
            cell.set_text_props(ha="left")
    return tab


def structure(ax, case, side=False, labels=True):
    xyz = case["xyz"] - np.mean(case["xyz"], axis=0)
    # Orthographic views of the actual saved coordinates; no geometry optimization.
    if side:
        angle = np.deg2rad(55.0)
        points = np.column_stack((xyz[:, 0], xyz[:, 1] * np.cos(angle) + xyz[:, 2] * np.sin(angle)))
    else:
        points = xyz[:, :2]
    bo = np.asarray(case["actual"]["bond_orders"])
    for i in range(len(points)):
        for j in range(i):
            if bo[i, j] > 0.3:
                ax.plot(*points[[i, j]].T, color="#a0aeb7", lw=3, zorder=1)
    for i, (symbol, point) in enumerate(zip(case["symbols"], points)):
        ax.scatter(
            *point,
            s=190 if symbol == "O" else 75,
            c=RED if symbol == "O" else "#e5edf2",
            edgecolors=INK,
            lw=0.6,
            zorder=3,
        )
        if labels:
            ax.annotate(str(i + 1), point, xytext=(7, 7), textcoords="offset points", fontsize=7)
    span = max(float(np.ptp(points, axis=0).max()), 2.0)
    center = (points.min(axis=0) + points.max(axis=0)) / 2
    ax.set_xlim(center[0] - 0.66 * span, center[0] + 0.66 * span)
    ax.set_ylim(center[1] - 0.66 * span, center[1] + 0.66 * span)
    ax.set_aspect("equal")
    ax.axis("off")
    if labels:
        x, y = center - 0.55 * span
        ax.plot([x, x + 1], [y, y], color=INK, lw=1.5)
        ax.text(x + 0.5, y - 0.08 * span, "1 Å", ha="center", fontsize=8)


def overview(pdf, cases, summary):
    fig = page(
        "Water clusters | Python vs LAMMPS",
        f"{len(cases)} isolated, neutral, unoptimized geometries · single-point energies and forces",
        1,
    )
    for i, case in enumerate(cases):
        ax = fig.add_axes([0.055 + i * 0.89 / len(cases), 0.65, 0.85 / len(cases), 0.225])
        structure(ax, case, labels=False)
        ax.set_title(case["name"].replace("_", " ").title(), fontsize=10)
    fig.text(0.055, 0.627, "TOTAL ENERGIES   /   kcal/mol; ΔE = Python − LAMMPS", fontsize=10, weight="bold")
    table(
        fig,
        [0.055, 0.424, 0.89, 0.187],
        ["Structure", "Atoms", "Python total", "LAMMPS total", "ΔE"],
        [
            [
                c["name"].replace("_", " "),
                c["values"]["atoms"],
                f"{c['actual']['energy']:.9f}",
                f"{c['reference']['energy']:.9f}",
                f"{c['values']['delta_energy_kcal_mol']:+.3e}",
            ]
            for c in cases
        ],
        [0.23, 0.07, 0.25, 0.25, 0.20],
    )
    fig.text(0.055, 0.388, "FORCE DIFFERENCES   /   kcal/mol/Å; Cartesian components", fontsize=10, weight="bold")
    table(
        fig,
        [0.055, 0.187, 0.89, 0.185],
        ["Structure", "Fixed-charge max |ΔF|", "Fixed-charge RMS ΔF"],
        [
            [
                c["name"].replace("_", " "),
                f"{c['values']['max_force_component_error_kcal_mol_A']:.3e}",
                f"{c['values']['rms_force_component_error_kcal_mol_A']:.3e}",
            ]
            for c in cases
        ],
        [0.34, 0.33, 0.33],
    )
    fig.text(
        0.055,
        0.15,
        "Fixed-charge forces: QEq at each structure; charges held constant only during differentiation (default).\n"
        "All differences subtract LAMMPS forces. RMS averages all 3N Cartesian components.",
        fontsize=8.5,
        linespacing=1.5,
        va="top",
    )
    version = next(iter(summary["cases"].values()))["reference_version"]
    fig.text(
        0.055,
        0.064,
        f"Reference: {version} via lmp_mpi. Parameters: {summary['force_field']}; Achtyl et al. (2015).",
        fontsize=8,
    )
    fig.text(0.055, 0.046, f"Parameter SHA-256: {summary['sha256']}", fontsize=7)
    pdf.savefig(fig)
    plt.close(fig)


def detail(pdf, case, number):
    v, actual, reference = case["values"], case["actual"], case["reference"]
    fig = page(
        case["name"].replace("_", " ").title(),
        f"{v['atoms']} atoms · O in red, H in pale gray · labels follow XYZ atom order · lines: bond order > 0.3",
        number,
    )
    for side, left in ((False, 0.05), (True, 0.285)):
        ax = fig.add_axes([left, 0.555, 0.225, 0.30])
        structure(ax, case, side=side)
        ax.set_title("XY projection" if not side else "Tilted projection (55° about x)", fontsize=9)
    table(
        fig,
        [0.55, 0.552, 0.395, 0.307],
        ["Quantity", "Value"],
        [
            ["Python E (kcal/mol)", f"{actual['energy']:.12f}"],
            ["LAMMPS E (kcal/mol)", f"{reference['energy']:.12f}"],
            ["ΔE (kcal/mol)", f"{v['delta_energy_kcal_mol']:+.6e}"],
            ["Fixed-charge max |ΔF|", f"{v['max_force_component_error_kcal_mol_A']:.6e}"],
            ["Fixed-charge RMS ΔF", f"{v['rms_force_component_error_kcal_mol_A']:.6e}"],
        ],
        [0.54, 0.46],
        size=8.5,
    )
    fig.text(0.055, 0.495, "ENERGY COMPONENTS   /   kcal/mol", fontsize=10, weight="bold")
    rows = [
        [key.replace("_", " "), f"{value:.8f}", f"{value - reference['components'][key]:+.2e}"]
        for key, value in actual["components"].items()
    ]
    table(fig, [0.055, 0.105, 0.43, 0.37], ["Component", "Python", "Δ vs LAMMPS"], rows, [0.38, 0.34, 0.28], size=8)
    ax = fig.add_axes([0.59, 0.19, 0.35, 0.275])
    indices = np.arange(1, v["atoms"] + 1)
    floor = 1e-16
    ax.semilogy(
        indices,
        np.maximum(np.max(abs(case["delta"]), axis=1), floor),
        "o-",
        ms=4,
        color=BLUE,
        label="Fixed-charge forces",
    )
    ax.set_xticks(indices)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("Atom index (XYZ order)")
    ax.set_ylabel("Max component |ΔF| (kcal/mol/Å)")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=8, loc="best")
    ax.set_title("Force differences from LAMMPS", fontsize=10, pad=12)
    fig.text(
        0.55,
        0.10,
        "Each point: max over x, y, z for one atom.\nLog scale; values below 10⁻¹⁶ shown at 10⁻¹⁶.\nForce summary above uses kcal/mol/Å throughout.",
        fontsize=8,
        linespacing=1.4,
    )
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "validation/water")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/water/diagnostics/report.pdf")
    args = parser.parse_args()
    summary = json.loads((args.input / "summary.json").read_text())
    if not summary.get("reference_verified"):
        parser.error("Input must contain verified LAMMPS reference results")
    structures = json.loads((args.input / "structures.json").read_text())
    with gzip.open(args.input / "results.json.gz", "rt") as stream:
        results = json.load(stream)
    cases = [
        read_case(name.removeprefix("cluster_water_"), structures[name], results[name])
        for name in summary["cases"]
        if structures[name]["cell"] is None
    ]
    if not cases:
        parser.error("Input contains no isolated water structures")
    first = next(iter(summary["cases"].values()))
    summary["force_field"], summary["sha256"] = first["force_field"], first["force_field_sha256"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(
        args.output,
        metadata={
            "Title": "Water clusters: Python ReaxFF vs LAMMPS",
            "Author": "xreac",
            "Subject": "Structures, total energies and force differences",
            "Keywords": "ReaxFF, water, LAMMPS, validation",
            "CreationDate": None,
            "ModDate": None,
        },
    ) as pdf:
        overview(pdf, cases, summary)
        for number, case in enumerate(cases, 2):
            detail(pdf, case, number)
    csv_path = args.output.with_suffix(".csv")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]["values"]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(case["values"] for case in cases)
    print(f"Created {args.output} ({len(cases) + 1} pages) and {csv_path}")


if __name__ == "__main__":
    main()
