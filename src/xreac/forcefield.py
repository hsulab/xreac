"""Reader for the standard four-line-atom ReaxFF parameter format.

Parameter conventions follow LAMMPS stable_22Jul2025_update4 (GPL-2.0+).
Element labels, interactions, and mixing rules are read from the parameter file.
"""
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations_with_replacement
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
        """Backward-compatible shortcut for the bundled Raymand 2010 ZnOH set."""
        return cls.bundled("ffield.reax.ZnOH")

    @classmethod
    def bundled(cls, name):
        """Load a named parameter file shipped in the repository's data directory."""
        if not isinstance(name, str) or Path(name).name != name or name in ("", ".", ".."):
            raise ValueError("A bundled parameter name must be a plain filename")
        # Wheels bundle the top-level data directory as xreac.data. Source
        # checkouts (including editable installs) read it directly from the repo.
        module = Path(__file__).resolve()
        bundled = module.with_name("data") / name
        if bundled.is_file():
            return cls.from_file(bundled)
        return cls.from_file(module.parents[2] / "data" / name)

    @property
    def elements(self):
        """Parameter-file atom labels, in their original order (including dummy types)."""
        return tuple(self.atoms)

    @property
    def vdw_type(self):
        """Standard LAMMPS vdW variant: 1 shielded, 2 inner wall, 3 both."""
        types = set()
        for a in self.atoms.values():
            shield = a["gamma_w"] > .5
            core = a["rcore"] > .01 and a["acore"] > .01
            types.add(3 if shield and core else 1 if shield else 2 if core else 0)
        if len(types) != 1 or 0 in types:
            raise ValueError("Inconsistent or unsupported van der Waals method across atom types")
        return types.pop()

    @classmethod
    def from_file(cls, path):
        path = Path(path).resolve()
        raw = path.read_bytes()
        lines = raw.decode().splitlines()
        if not lines:
            raise ValueError(f"Empty ReaxFF parameter file: {path}")
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
            g = np.array([floats(row()[:1], 1)[0] for _ in range(ng)])
            na = count()
            if na == 0:
                raise ValueError("A force field must define at least one atom type")
            for _ in range(3):
                next(rows)
            atoms, names = {}, []
            for _ in range(na):
                first = row()
                name = first[0]
                if any(name.casefold() == existing.casefold() for existing in atoms):
                    raise ValueError(f"Duplicate element {name}")
                vals = floats(first[1:] + row() + row() + row(), 32)
                a = dict(zip(ATOM_NAMES, vals))
                a["eta"] *= 2
                a["p_hbond"] = int(a["p_hbond"])
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
                second = row()
                # Standard files may omit the last unused bond parameter.
                if len(second) == 7:
                    second.append("0")
                vals = floats(first[2:] + second, 16)
                p = dict(zip(BOND_NAMES, vals))
                pairs[key] = pairs[key[::-1]] = p
            # Nonbonded mixing is defined for every pair, independently of
            # whether an explicit bond record exists (including dummy types).
            for s, t in combinations_with_replacement(names, 2):
                key = (s, t)
                p = pairs.get(key, dict.fromkeys(BOND_NAMES, 0.0))
                a, b = atoms[s], atoms[t]
                for out, field in [("r_s", "r_s"), ("r_p", "r_pi"), ("r_pp", "r_pi_pi")]:
                    p[out] = (a[field] + b[field]) / 2
                for out, field in [("p_boc3", "b_o_132"), ("p_boc4", "b_o_131"),
                                   ("p_boc5", "b_o_133"), ("D", "epsilon"),
                                   ("alpha", "alpha"), ("gamma_w", "gamma_w"),
                                   ("r_vdW", "r_vdw"), ("rcore", "rcore"),
                                   ("ecore", "ecore"), ("acore", "acore")]:
                    if a[field]*b[field] < 0:
                        raise ValueError(f"Invalid mixing parameters for {s}/{t}: {field}")
                    p[out] = np.sqrt(a[field] * b[field])
                p["r_vdW"] *= 2
                if a["gamma"]*b["gamma"] <= 0:
                    raise ValueError("Charge shielding parameters must be positive")
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
            explicit_torsions = set()
            for _ in range(count()):
                values = row()
                key = indices(values[:4], wildcard=True)
                if len(values[4:]) not in (5, 6, 7):
                    raise ValueError("Expected 5 to 7 torsion parameters")
                vals = floats(values[4:], len(values[4:]))[:5]
                if "*" not in key:
                    torsions[key] = torsions[key[::-1]] = vals
                    explicit_torsions.update((key, key[::-1]))
                else:
                    if key[0] != "*" or key[3] != "*" or "*" in key[1:3]:
                        raise ValueError("Only terminal 0-i-j-0 torsion wildcards are supported")
                    # LAMMPS: explicit entries always win; later wildcard rows
                    # replace earlier wildcard defaults for the same pair.
                    for s in names:
                        for t in names:
                            expanded = (s, key[1], key[2], t)
                            for k in (expanded, expanded[::-1]):
                                if k not in explicit_torsions:
                                    torsions[k] = vals
            hydrogen = {}
            try:
                nh = count()
            except StopIteration:
                nh = 0  # LAMMPS permits omission of the final hydrogen-bond block.
            for _ in range(nh):
                values = row()
                hydrogen[indices(values[:3])] = floats(values[3:], 4)
            if any(line.split("!")[0].split("#")[0].strip() for line in rows):
                raise ValueError("Unsupported trailing parameter section")
            if not np.isfinite(g).all():
                raise ValueError("Non-finite general parameters")
        except (StopIteration, IndexError, ValueError) as exc:
            raise ValueError(f"Invalid or unsupported ReaxFF file {path}: {exc}") from exc
        return cls(path, sha256(raw).hexdigest(), citation, g, atoms, pairs, angles, torsions, hydrogen)

    def validate_model(self, elements=None):
        """Conservatively reject models not covered by the implemented equations."""
        self.vdw_type
        if self.general[11] != 0 or self.general[12] <= 0:
            raise ValueError("Only zero lower taper radius and a positive upper cutoff are supported")
        if self.general[1] <= 0 or self.general[28] <= 0 or self.general[29] <= 0:
            raise ValueError("Invalid bond-order correction, vdW exponent, or bond-order cutoff")
        if elements is None:
            return
        missing = set(elements)-self.atoms.keys()
        if missing:
            raise ValueError(f"Atom types absent from the force field: {', '.join(sorted(missing))}")
        for s in set(elements):
            a = self.atoms[s]
            if a["eta"] <= 0 or a["gamma"] <= 0:
                raise ValueError("Invalid QEq hardness or shielding")
            if a["mass"] <= 0 or a["valency"] < 0:
                raise ValueError("Positive atomic masses and nonnegative valencies are required")
            for t in set(elements):
                p = self.pairs[s, t]
                if p["r_vdW"] <= 0:
                    raise ValueError(f"Invalid van der Waals radius for {s}/{t}")
                for atom_radius, pair_radius in (("r_s", "r_s"), ("r_pi", "r_p"), ("r_pi_pi", "r_pp")):
                    if a[atom_radius] > 0 and self.atoms[t][atom_radius] > 0 and p[pair_radius] <= 0:
                        raise ValueError(f"Invalid {pair_radius} for {s}/{t}")
