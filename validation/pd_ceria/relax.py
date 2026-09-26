"""Relax initial.xyz with LAMMPS once the complete published potential is available.

Run from any directory: python validation/pd_ceria/relax.py [--lammps lmp_serial]
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from ase.io import read, write
from ase.units import kcal, mol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", default="lmp_serial")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    potential = here.parents[1] / "data" / "ffield.PdCeO"
    if not potential.is_file():
        parser.error(
            "data/ffield.PdCeO is unavailable: the complete public candidate has unresolved provenance "
            "and Ce/O conflicts. Resolve these before running; see parameters.md."
        )
    initial = read(here / "initial.xyz")
    fixed = initial.constraints[0].get_indices()
    mobile = np.ones(len(initial), dtype=bool)
    mobile[fixed] = False
    with tempfile.TemporaryDirectory(prefix="pd_ceria_") as temporary:
        work = Path(temporary)
        shutil.copyfile(potential, work / "ffield.PdCeO")
        write(work / "initial.data", initial, format="lammps-data", atom_style="charge", specorder=["Ce", "O", "Pd"])
        (work / "in.relax").write_text(
            "units real\natom_style charge\nboundary p p f\nread_data initial.data\n"
            "pair_style reaxff NULL safezone 3 mincap 100\npair_coeff * * ffield.PdCeO Ce O Pd\n"
            "neighbor 2.0 bin\nneigh_modify every 1 delay 0 check yes\n"
            "fix charge all qeq/reaxff 1 0.0 10.0 1.0e-10 reaxff maxiter 1000\n"
            f"group bottom id {' '.join(str(i + 1) for i in fixed)}\n"
            "fix freeze bottom setforce 0.0 0.0 0.0\n"
            "thermo 100\nthermo_style custom step pe fmax fnorm\nthermo_modify lost error flush yes\n"
            'run 0\nvariable energy equal pe\nprint "INITIAL_ENERGY $(v_energy:%.16g)"\n'
            "timestep 0.1\nmin_style fire\nmin_modify dmax 0.05\n"
            f"minimize 0.0 {0.02 / (kcal / mol):.12g} 20000 200000\n"
            'run 0\nprint "FINAL_ENERGY $(v_energy:%.16g)"\n'
            "write_dump all custom relaxed.dump id type x y z q fx fy fz modify sort id format float %.16g\n"
        )
        try:
            subprocess.run([args.lammps, "-in", "in.relax", "-screen", "none"], cwd=work, check=True)
        finally:
            if (work / "log.lammps").exists():
                shutil.copyfile(work / "log.lammps", here / "log.lammps")
        relaxed = read(work / "relaxed.dump", format="lammps-dump-text", specorder=["Ce", "O", "Pd"], units="real")
        forces = relaxed.get_forces()
        energy = {}
        for line in (work / "log.lammps").read_text().splitlines():
            for label in ("INITIAL_ENERGY", "FINAL_ENERGY"):
                if line.startswith(label + " "):
                    energy[label] = float(line.split()[1]) * kcal / mol
        # Carry constraints and trilayer labels through the LAMMPS dump.
        relaxed.set_tags(initial.get_tags())
        relaxed.set_constraint(initial.constraints)
        assert np.allclose(relaxed.positions[fixed], initial.positions[fixed], atol=1e-5, rtol=0)
        maximum_force = float(np.linalg.norm(forces[mobile], axis=1).max())
        relaxed.info.update(energy)
        relaxed.info["maximum_mobile_force_eV_per_A"] = maximum_force
        relaxed.info["converged"] = maximum_force <= 0.02
        write(here / "relaxed.xyz", relaxed, format="extxyz")

    print(f"Initial energy: {energy['INITIAL_ENERGY']:.8f} eV")
    print(f"Final energy: {energy['FINAL_ENERGY']:.8f} eV")
    print(f"Maximum mobile force: {maximum_force:.6f} eV/Angstrom")
    symbols = np.array(relaxed.get_chemical_symbols())
    distances = relaxed.get_all_distances(mic=True)
    pd_distances = distances[np.ix_(symbols == "Pd", symbols == "Pd")][np.triu_indices(4, 1)]
    print("Pd-Pd distances (Angstrom):", np.round(pd_distances, 5).tolist())
    for a, b in [("Pd", "O"), ("Pd", "Ce"), ("Ce", "O")]:
        print(f"Minimum {a}-{b}: {distances[np.ix_(symbols == a, symbols == b)].min():.5f} Angstrom")
    print(f"Pd height span: {np.ptp(relaxed.positions[symbols == 'Pd', 2]):.5f} Angstrom")
    print(
        f"Lowest Pd minus highest O: {relaxed.positions[symbols == 'Pd', 2].min() - relaxed.positions[symbols == 'O', 2].max():.5f} Angstrom"
    )
    if maximum_force > 0.02:
        raise SystemExit("Relaxation did not meet the 0.02 eV/Angstrom mobile-force target; inspect log.lammps.")


if __name__ == "__main__":
    main()
