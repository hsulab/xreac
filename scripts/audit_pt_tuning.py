"""Audit a Pt/O/H experimental variant on retained common geometries.

Compare fixed-geometry energies rather than confounding parameter changes
with relaxation. Check charge-response forces against energy differences.
"""

import argparse
import json
from pathlib import Path

from ase.io import read
from ase.mep import NEB

from neb_pt_water import ForceField, ReaxFFCalculator, np
from screen_pt_water import closest_pt_h_cutoff


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ffield", type=Path, required=True)
    parser.add_argument("--band", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    original = ForceField.bundled("ffield.reax.PtNiCHO.2016")
    tuned = ForceField.from_file(args.ffield)
    # Numerical tables must be identical except for the chosen exponent and
    # its symmetry-related reverse entry. Retain a machine-readable diff.
    changes = []
    np.testing.assert_array_equal(original.general, tuned.general)
    assert original.atoms == tuned.atoms and original.pairs == tuned.pairs
    for table in ("angles", "torsions", "hydrogen_bonds"):
        a, b = getattr(original, table), getattr(tuned, table)
        assert a.keys() == b.keys()
        for key in a:
            x, y = np.array(a[key]), np.array(b[key])
            for index in zip(*np.where(x != y)):
                changes.append(
                    dict(
                        table=table,
                        key=key,
                        index=tuple(map(int, index)),
                        before=float(x[index]),
                        after=float(y[index]),
                    )
                )
    assert len(changes) == 2 and all(
        c["key"] in (("O", "H", "Pt"), ("Pt", "H", "O")) and c["index"] == (0, 6) for c in changes
    )
    rows = []
    images = read(args.band, index=":")
    metal = images[0][0].symbol
    for i, a in enumerate(images):
        a.calc = ReaxFFCalculator(original, full_derivative=True)
        old = a.get_potential_energy()
        a.calc = ReaxFFCalculator(tuned, full_derivative=True)
        new = a.get_potential_energy()
        rows.append(dict(image=i, original_energy_ev=old, tuned_energy_ev=new, difference_ev=new - old))
    audits = []
    rng = np.random.default_rng(20260925)
    for i in (0, int(np.argmax([r["tuned_energy_ev"] for r in rows])), len(images) - 1):
        a = images[i].copy()
        a.calc = ReaxFFCalculator(tuned, full_derivative=True)
        x = a.positions.copy()
        direction = rng.normal(size=x.shape)
        direction[:4] = 0  # same fixed bottom layer
        direction /= np.linalg.norm(direction)
        analytic = -float(np.sum(a.get_forces() * direction))
        fd = []
        for h in (1e-5, 1e-6, 1e-7):
            a.positions = x + h * direction
            plus = a.get_potential_energy()
            a.positions = x - h * direction
            minus = a.get_potential_energy()
            derivative = (plus - minus) / (2 * h)
            fd.append(
                dict(step_angstrom=h, derivative_ev_angstrom=derivative, error_ev_angstrom=abs(derivative - analytic))
            )
        a.positions = x
        audits.append(
            dict(
                image=i,
                analytic_ev_angstrom=analytic,
                finite_differences=fd,
                **{f"closest_{metal.lower()}_h_cutoff": closest_pt_h_cutoff(a, tuned, metal=metal)},
            )
        )
    controls = []
    n_pt = len(images[0]) - 3
    assert images[0].get_chemical_symbols() == [metal] * n_pt + ["O", "H", "H"]
    for name, selection in (
        ("bare_slab", slice(0, n_pt)),
        ("water", slice(n_pt, n_pt + 3)),
        ("slab_plus_H", list(range(n_pt)) + [n_pt + 2]),
    ):
        a = images[0][selection]
        a.set_constraint()
        energies = []
        forces = []
        for ff in (original, tuned):
            a.calc = ReaxFFCalculator(ff, full_derivative=True)
            energies.append(a.get_potential_energy())
            forces.append(a.get_forces())
        controls.append(
            dict(
                name=name,
                energy_change_ev=energies[1] - energies[0],
                max_force_change_ev_angstrom=float(np.max(abs(forces[1] - forces[0]))),
            )
        )
    neb = NEB(images, k=0.1, method="improvedtangent", climb=True)
    neb_fmax = float(np.linalg.norm(neb.get_forces(), axis=1).max())
    peak = int(np.argmax([r["tuned_energy_ev"] for r in rows]))
    report = dict(
        original_sha256=original.checksum,
        tuned_sha256=tuned.checksum,
        band=str(args.band),
        parameter_changes=changes,
        common_geometry_energies=rows,
        force_audits=audits,
        controls=controls,
        independently_recomputed_ci_neb_fmax_ev_angstrom=neb_fmax,
        peak_physical_fmax_ev_angstrom=float(np.linalg.norm(images[peak].get_forces(), axis=1).max()),
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
