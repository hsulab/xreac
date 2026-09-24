"""Optional ASE adapter. Core xreac imports do not require ASE."""

from dataclasses import replace
from numbers import Real

import numpy as np

try:
    from ase.calculators.calculator import Calculator as ASECalculator, all_changes
    from ase.neighborlist import PrimitiveNeighborList
    from ase.units import kcal, mol
except ImportError as exc:
    raise ImportError("Install xreac[ase] to use the ASE calculator") from exc

from .calculator import Calculator
from .geometry import Boundary, validate_expansion_limit

KCAL_MOL_TO_EV = kcal / mol


def _directed_neighbors(neighbor_list):
    """Flatten ASE's half list and add reverse edges without Python pair loops."""
    counts = [len(row) for row in neighbor_list.neighbors]
    i = np.repeat(np.arange(len(counts), dtype=int), counts)
    j = np.concatenate(neighbor_list.neighbors) if counts else np.empty(0, dtype=int)
    shifts = np.concatenate(neighbor_list.displacements).reshape(-1, 3) if counts else np.empty((0, 3), dtype=int)
    return np.concatenate((i, j)), np.concatenate((j, i)), np.concatenate((shifts, -shifts))


class ReaxFFCalculator(ASECalculator):
    """ReaxFF for neutral ASE Atoms, with optional fixed-cell periodicity.

    Energies and forces use ASE units (eV, eV/A); charges and dipoles use e
    and e A. ``evaluation`` retains the last core result in kcal/mol units.
    Fixed-charge forces are the default; full_derivative=True is an explicit
    option for single-point charge-response calculations.
    neighbor_backend="ase" uses ASE's image-resolved neighbor list by default.
    Set it to "replicated" to use the native dense image/replication method.
    neighbor_skin is ASE's per-atom displacement allowance in Angstrom;
    zero forces a fresh list on every evaluation. Only topology is reused.
    """

    implemented_properties = ["energy", "forces", "charges", "dipole"]
    default_parameters = {
        "full_derivative": False,
        "total_charge": 0,
        "max_expanded_atoms": 512,
        "neighbor_backend": "ase",
        "neighbor_skin": 0.3,
    }

    def __init__(
        self,
        force_field,
        *,
        full_derivative=False,
        total_charge=0,
        max_expanded_atoms=512,
        neighbor_backend="ase",
        neighbor_skin=0.3,
        **kwargs,
    ):
        self.core = Calculator(force_field, max_expanded_atoms=max_expanded_atoms)
        self.evaluation = None
        self._neighbor_list = None
        self._neighbor_key = None
        self._neighbors = None
        self.neighbor_list_builds = 0
        super().__init__(
            full_derivative=full_derivative,
            total_charge=total_charge,
            max_expanded_atoms=max_expanded_atoms,
            neighbor_backend=neighbor_backend,
            neighbor_skin=neighbor_skin,
            **kwargs,
        )

    def set(self, **kwargs):
        unknown = set(kwargs) - set(self.default_parameters)
        if unknown:
            raise ValueError(f"Unknown ReaxFF calculator parameters: {sorted(unknown)}")
        if "full_derivative" in kwargs and not isinstance(kwargs["full_derivative"], bool):
            raise ValueError("full_derivative must be a boolean")
        if kwargs.get("total_charge", 0) != 0:
            raise ValueError("Only neutral systems are supported")
        if "max_expanded_atoms" in kwargs:
            validate_expansion_limit(kwargs["max_expanded_atoms"])
        if kwargs.get("neighbor_backend", "ase") not in ("ase", "replicated"):
            raise ValueError("neighbor_backend must be 'replicated' or 'ase'")
        if "neighbor_skin" in kwargs:
            skin = kwargs["neighbor_skin"]
            if isinstance(skin, bool) or not isinstance(skin, Real) or not np.isfinite(skin) or skin < 0:
                raise ValueError("neighbor_skin must be a finite nonnegative number")
        changed = super().set(**kwargs)
        if changed:
            self.core.max_expanded_atoms = self.parameters.max_expanded_atoms
            self.reset()
        return changed

    def reset(self):
        super().reset()
        self.evaluation = None
        self._neighbor_list = None
        self._neighbor_key = None
        self._neighbors = None

    def _ase_neighbors(self):
        """Update ASE topology before differentiation; recompute physics in core."""
        atoms = self.atoms
        Boundary(atoms.cell.array, atoms.pbc)
        cutoff = np.nextafter(float(self.core.force_field.general[12]), np.inf)
        skin = self.parameters.neighbor_skin
        key = (cutoff, skin, tuple(atoms.numbers))
        if self._neighbor_list is None or key != self._neighbor_key or skin == 0:
            self._neighbor_list = PrimitiveNeighborList(
                np.full(len(atoms), cutoff / 2), skin=skin, self_interaction=False, bothways=False
            )
            self._neighbor_key = key
        try:
            # ASE adds skin to each radius: pair cutoff grows by 2*skin.
            # It rebuilds after any atom moves more than skin from its build
            # position, or cell/PBC changes. Wrapped jumps rebuild as well.
            rebuilt = self._neighbor_list.update(atoms.pbc, atoms.cell, atoms.positions)
        except Exception:
            self._neighbor_list = None
            self._neighbors = None
            raise
        if rebuilt:
            self._neighbors = _directed_neighbors(self._neighbor_list)
            self.neighbor_list_builds += 1
        return self._neighbors

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {}
        self.evaluation = None
        if self.atoms is None:
            raise ValueError("An ASE Atoms object is required")
        initial_charges = self.atoms.get_initial_charges()
        if not np.isfinite(initial_charges).all() or abs(initial_charges.sum()) > 1e-8:
            raise ValueError("Initial charges must be finite and sum to zero; only neutral systems are supported")
        neighbors = None
        if self.parameters.neighbor_backend == "ase":
            # Build integer topology here, before entering the core evaluator.
            # The core only differentiates x[j] - x[i] + S @ cell, never this call.
            neighbors = self._ase_neighbors()
        result = self.core.evaluate(
            self.atoms.get_chemical_symbols(),
            self.atoms.positions,
            cell=self.atoms.cell.array,
            pbc=self.atoms.pbc,
            total_charge=self.parameters.total_charge,
            full_derivative=self.parameters.full_derivative,
            neighbors=neighbors,
        )
        # Only the adapter knows which builder supplied the core's arrays.
        result = replace(result, neighbor_backend=self.parameters.neighbor_backend)
        self.evaluation = result
        self.results = {
            "energy": result.energy * KCAL_MOL_TO_EV,
            "forces": result.forces * KCAL_MOL_TO_EV,
            "charges": result.charges.copy(),
            "dipole": result.dipole.copy(),
        }


def relax_with_ase(core, symbols, positions, force_tolerance, max_iterations, *, cell=None, pbc=False):
    """Internal bridge to ASE FIRE; retain xreac's Cartesian-component stop rule."""
    from ase import Atoms
    from ase.optimize import FIRE

    atoms = Atoms(symbols, positions=positions, cell=cell, pbc=pbc)
    adapter = ReaxFFCalculator(core.force_field, full_derivative=False, max_expanded_atoms=core.max_expanded_atoms)
    # Reuse the caller's calculator, including any instrumentation/subclass.
    adapter.core = core
    atoms.calc = adapter
    # ASE uses max atomic vector norm for fmax. Observe each iteration ourselves
    # to preserve xreac's existing max Cartesian component criterion.
    with FIRE(atoms, logfile=None, downhill_check=False) as optimizer:
        for _ in optimizer.irun(fmax=force_tolerance * KCAL_MOL_TO_EV, steps=max_iterations):
            result = adapter.evaluation
            if np.max(np.abs(result.forces)) <= force_tolerance:
                return atoms.positions.copy(), result, True, optimizer.nsteps
        return atoms.positions.copy(), adapter.evaluation, False, optimizer.nsteps
