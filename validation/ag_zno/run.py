"""Run bulk first, then one Ag4/ZnO minimization. Requires ASE and lmp_serial.

From the repository root: python validation/ag_zno/run.py
The original ffield is used for both minimizations; a temporary xreac-written
file is also checked with a LAMMPS bulk single point before either relaxation.
"""

import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.io import read, write
from ase.units import kcal, mol

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "scripts")]
from xreac import ForceField
from validate import zno_bulk_case
from build import build

EV = kcal / mol


def same(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            same(a[k], b[k])
    else:
        np.testing.assert_array_equal(a, b)


def numerical_tokens(path):
    """Independent audit of every number in this source file (including section counts)."""
    values = []
    for line in path.read_text().splitlines()[1:]:
        for token in line.split("!")[0].split("#")[0].split():
            try:
                values.append(float(token.replace("D", "E").replace("d", "e")))
            except ValueError:
                pass  # Atom labels and column headings are not numeric fields.
    return values


def calculate(atoms, prefix, potential, minimize=True, directory=HERE):
    data, inp, log, dump = (prefix + suffix for suffix in ("initial.data", "in.lammps", "log.lammps", "relaxed.dump"))
    order = [s for s in ("Zn", "O", "Ag") if s in atoms.get_chemical_symbols()]
    write(directory / data, atoms, format="lammps-data", atom_style="charge", specorder=order, masses=True)
    fixed = atoms.constraints[0].get_indices() if atoms.constraints else np.array([], dtype=int)
    commands = (
        "units real\natom_style charge\n"
        f"boundary {' '.join('p' if p else 'f' for p in atoms.pbc)}\nread_data {data}\n"
        "pair_style reaxff NULL safezone 3 mincap 100\n"
        f"pair_coeff * * {potential} {' '.join(order)}\n"
        "neighbor 2.0 bin\nneigh_modify every 1 delay 0 check yes\n"
        "fix charge all qeq/reaxff 1 0.0 10.0 1.0e-10 reaxff maxiter 1000\n"
    )
    if len(fixed):
        commands += f"group bottom id {' '.join(str(i + 1) for i in fixed)}\nfix freeze bottom setforce 0 0 0\n"
    commands += (
        "thermo 100\nthermo_style custom step pe fmax fnorm\nthermo_modify lost error flush yes\n"
        'run 0\nprint "INITIAL_ENERGY $(pe:%.16g)"\n'
    )
    if minimize:
        commands += f"timestep 0.1\nmin_style fire\nmin_modify dmax 0.05\nminimize 0.0 {0.01 / EV:.14g} 20000 200000\n"
    commands += (
        'run 0\nprint "FINAL_ENERGY $(pe:%.16g)"\n'
        f"write_dump all custom {dump} id type x y z q fx fy fz modify sort id format float %.16g\n"
    )
    (directory / inp).write_text(commands)
    process = subprocess.run(["lmp_serial", "-in", inp, "-log", log, "-screen", "none"], cwd=directory)
    if process.returncode:
        raise RuntimeError((directory / log).read_text()[-4000:])
    relaxed = read(directory / dump, format="lammps-dump-text", specorder=order, units="real")
    forces = relaxed.get_forces()
    mobile = np.ones(len(atoms), dtype=bool)
    mobile[fixed] = False
    np.testing.assert_allclose(relaxed.positions[fixed], atoms.positions[fixed], atol=1e-7, rtol=0)
    report = {}
    for line in (directory / log).read_text().splitlines():
        for label in ("INITIAL_ENERGY", "FINAL_ENERGY"):
            if line.startswith(label + " "):
                report[label.lower() + "_eV"] = float(line.split()[1]) * EV
    report["maximum_mobile_force_eV_A"] = float(np.linalg.norm(forces[mobile], axis=1).max())
    relaxed.set_tags(atoms.get_tags())
    relaxed.set_constraint(atoms.constraints)
    relaxed.info.update(report)
    write(directory / (prefix + "relaxed.xyz"), relaxed)
    (directory / dump).unlink()
    return relaxed, report


def main():
    potential = ROOT / "data" / "ffield.AgZnO"
    ff = ForceField.from_file(potential)
    summary = {"elements": ff.elements, "counts": ff.section_counts, "sha256": ff.checksum}
    _, symbols, positions, cell, pbc = zno_bulk_case()
    bulk = Atoms(symbols, positions=positions, cell=cell, pbc=pbc)
    write(HERE / "bulk.initial.xyz", bulk)
    with tempfile.TemporaryDirectory(prefix="ag_zno_roundtrip_") as temporary:
        work = Path(temporary)
        exported = work / "ffield.roundtrip"
        ff.to_file(exported)
        again = ForceField.from_file(exported)
        np.testing.assert_array_equal(numerical_tokens(potential), numerical_tokens(exported))
        summary["roundtrip_numerical_tokens"] = len(numerical_tokens(potential))
        assert ff._sections == again._sections
        for name in ("general", "atoms", "pairs", "angles", "torsions", "hydrogen_bonds"):
            same(getattr(ff, name), getattr(again, name))
        original, original_report = calculate(bulk, "original.", potential, False, work)
        written, written_report = calculate(bulk, "written.", exported, False, work)
        np.testing.assert_allclose(original.get_forces(), written.get_forces(), atol=1e-10, rtol=0)
        np.testing.assert_allclose(original.get_initial_charges(), written.get_initial_charges(), atol=1e-12, rtol=0)
        assert original_report == written_report
        summary["roundtrip"] = "PASS: all source records and model tables exact; LAMMPS energy/forces/charges agree"
    print(json.dumps(summary, indent=2), flush=True)
    relaxed, report = calculate(bulk, "bulk.", "../../data/ffield.AgZnO")
    symbols = np.array(relaxed.get_chemical_symbols())
    zn_o = relaxed.get_all_distances(mic=True)[np.ix_(symbols == "Zn", symbols == "O")]
    report["Zn_O_nearest_range_A"] = [float(zn_o.min(axis=1).min()), float(zn_o.min(axis=1).max())]
    report["Zn_coordination_range_2.5_A"] = [int((zn_o < 2.5).sum(axis=1).min()), int((zn_o < 2.5).sum(axis=1).max())]
    summary["bulk"] = report
    (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Bulk:", report, flush=True)
    if report["maximum_mobile_force_eV_A"] > 0.01 or not (zn_o.min() > 1.5 and np.all((zn_o < 2.5).sum(axis=1) == 4)):
        raise SystemExit("Bulk failed the convergence/structure check; stopped before building Ag4/ZnO.")
    slab = build()
    relaxed, report = calculate(slab, "", "../../data/ffield.AgZnO")
    symbols = np.array(relaxed.get_chemical_symbols())
    distances = relaxed.get_all_distances(mic=True)
    report["Ag_Ag_A"] = distances[np.ix_(symbols == "Ag", symbols == "Ag")][np.triu_indices(4, 1)].tolist()
    for a, b in (("Ag", "O"), ("Ag", "Zn"), ("Zn", "O")):
        report[f"minimum_{a}_{b}_A"] = float(distances[np.ix_(symbols == a, symbols == b)].min())
    report["Ag_height_span_A"] = float(np.ptp(relaxed.positions[symbols == "Ag", 2]))
    report["lowest_Ag_minus_highest_oxide_A"] = float(
        relaxed.positions[symbols == "Ag", 2].min() - relaxed.positions[symbols != "Ag", 2].max()
    )
    _, displacement = find_mic(relaxed.positions - slab.positions, slab.cell, slab.pbc)
    report["maximum_oxide_displacement_A"] = float(displacement[symbols != "Ag"].max())
    summary["surface"] = report
    (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Surface:", report, flush=True)
    if report["maximum_mobile_force_eV_A"] > 0.01:
        raise SystemExit("Surface force target not reached; inspect log.lammps.")


if __name__ == "__main__":
    main()
