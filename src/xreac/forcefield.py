"""Reader for the standard four-line-atom ReaxFF parameter format.

Parameter conventions follow LAMMPS stable_22Jul2025_update4 (GPL-2.0+).
Only the bundled Zn/O model has been validated against the reference engine.
"""
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np

ATOM_NAMES = (
    "r_s valency mass r_vdw epsilon gamma r_pi valency_e "
    "alpha gamma_w valency_boc p_ovun5 unused12 chi eta p_hbond "
    "r_pi_pi p_lp2 heat b_o_131 b_o_132 b_o_133 bcut_acks2 unused23 "
    "p_ovun2 p_val3 unused26 valency_val p_val5 rcore ecore acore"
).split()
BOND_NAMES = (
    "De_s De_p De_pp p_be1 p_bo5 v13cor p_bo6 p_ovun1 "
    "p_be2 p_bo3 p_bo4 unused11 p_bo1 p_bo2 ovc unused15"
).split()


@dataclass(frozen=True)
class ForceField:
    path: Path
    checksum: str
    citation: str
    general: np.ndarray
    atoms: dict
    pairs: dict
    angles: dict
    torsions: dict
    hydrogen_bonds: dict

    @classmethod
    def zno(cls):
        """Load the bundled Raymand 2010 ZnOH set (Zn/O calculations only)."""
        # Wheels bundle the top-level data directory as xreac.data. Source
        # checkouts (including editable installs) read it directly from the repo.
        module = Path(__file__).resolve()
        bundled = module.with_name("data") / "ffield.reax.ZnOH"
        if bundled.is_file():
            return cls.from_file(bundled)
        return cls.from_file(module.parents[2] / "data" / "ffield.reax.ZnOH")

    @classmethod
    def from_file(cls, path):
        path = Path(path).resolve()
        raw = path.read_bytes()
        lines = raw.decode().splitlines()
        citation = lines.pop(0)
        rows = iter(lines)

        def row():
            while True:
                values = next(rows).split("!")[0].split("#")[0].split()
                if values:
                    return values

        def count():
            n = int(row()[0])
            if n < 0:
                raise ValueError("Negative section count")
            return n

        def floats(values, n):
            if len(values) != n:
                raise ValueError(f"Expected {n} parameters, received {len(values)}")
            arr = np.array([float(x.replace("D", "E").replace("d", "e")) for x in values])
            if not np.isfinite(arr).all():
                raise ValueError("Non-finite force-field parameters")
            return arr

        try:
            ng = count()
            if ng != 39:
                raise ValueError("Only the standard 39-global-parameter format is supported")
            g = np.array([float(row()[0]) for _ in range(ng)])
            na = count()
            for _ in range(3):
                next(rows)
            atoms, names = {}, []
            for _ in range(na):
                first = row()
                name = first[0].capitalize()
                if name in atoms:
                    raise ValueError(f"Duplicate element {name}")
                vals = floats(first[1:] + row() + row() + row(), 32)
                a = dict(zip(ATOM_NAMES, vals))
                a["eta"] *= 2
                if a["mass"] < 21:
                    a["valency_val"] = a["valency_boc"]
                atoms[name] = a
                names.append(name)

            def indices(tokens, wildcard=False):
                result = []
                for token in tokens:
                    i = int(token)
                    if i == 0 and wildcard:
                        result.append("*")
                    elif 1 <= i <= len(names):
                        result.append(names[i-1])
                    else:
                        raise ValueError("Invalid element index")
                return tuple(result)

            nb = count()
            next(rows)
            pairs = {}
            for _ in range(nb):
                first = row()
                key = indices(first[:2])
                vals = floats(first[2:] + row(), 16)
                p = dict(zip(BOND_NAMES, vals))
                a, b = [atoms[s] for s in key]
                for out, field in [("r_s", "r_s"), ("r_p", "r_pi"), ("r_pp", "r_pi_pi")]:
                    p[out] = (a[field] + b[field]) / 2
                for out, field in [("p_boc3", "b_o_132"), ("p_boc4", "b_o_131"),
                                   ("p_boc5", "b_o_133"), ("D", "epsilon"),
                                   ("alpha", "alpha"), ("gamma_w", "gamma_w"),
                                   ("r_vdW", "r_vdw")]:
                    p[out] = np.sqrt(a[field] * b[field])
                p["r_vdW"] *= 2
                p["gamma"] = (a["gamma"] * b["gamma"]) ** -1.5
                pairs[key] = pairs[key[::-1]] = p
            for _ in range(count()):
                values = row()
                key = indices(values[:2])
                vals = floats(values[2:], 6)
                if key not in pairs:
                    raise ValueError("Off-diagonal entry without bond parameters")
                for field, value in zip(("D", "r_vdW", "alpha", "r_s", "r_p", "r_pp"), vals):
                    if value > 0:
                        pairs[key][field] = value * (2 if field == "r_vdW" else 1)
            angles = {}
            for _ in range(count()):
                values = row()
                key = indices(values[:3])
                vals = floats(values[3:], 7)
                for k in {key, key[::-1]}:
                    angles.setdefault(k, []).append(vals)
            torsions = {}
            for _ in range(count()):
                values = row()
                key = indices(values[:4], wildcard=True)
                vals = floats(values[4:], 7)
                if "*" in key:
                    raise ValueError("Wildcard torsion parameters are not supported")
                torsions[key] = torsions[key[::-1]] = vals[:5]
            hydrogen = {}
            for _ in range(count()):
                values = row()
                hydrogen[indices(values[:3])] = floats(values[3:], 4)
            if any(line.split("!")[0].split("#")[0].strip() for line in rows):
                raise ValueError("Unsupported trailing parameter section")
            if not np.isfinite(g).all():
                raise ValueError("Non-finite general parameters")
        except (StopIteration, IndexError, ValueError) as exc:
            raise ValueError(f"Invalid or unsupported ReaxFF file {path}: {exc}") from exc
        return cls(path, sha256(raw).hexdigest(), citation, g, atoms, pairs, angles, torsions, hydrogen)

    def validate_model(self):
        """Conservatively reject models not covered by the implemented equations."""
        if not {"Zn", "O"} <= self.atoms.keys():
            raise ValueError("The force field must contain Zn and O")
        for s in ("Zn", "O"):
            a = self.atoms[s]
            if a["gamma_w"] <= .5 or a["rcore"] != 0 or a["ecore"] != 0:
                raise ValueError("Only shielded vdW without an inner wall is supported")
            if a["eta"] <= 0 or a["gamma"] <= 0:
                raise ValueError("Invalid QEq hardness or shielding")
            if a["r_s"] <= 0 or a["valency"] <= 0:
                raise ValueError("Positive sigma radius and valency are required")
            for t in ("Zn", "O"):
                if (s, t) not in self.pairs:
                    raise ValueError(f"Missing bond parameters for {s}/{t}")
        if self.general[11] != 0 or self.general[12] <= 5 or self.general[37] != 0:
            raise ValueError("Unsupported taper or triple-bond stabilization variant")
