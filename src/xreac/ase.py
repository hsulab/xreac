"""Optional ASE adapter. Core xreac imports do not require ASE."""

from dataclasses import replace

import numpy as np

try:
    from ase.calculators.calculator import Calculator as ASECalculator, all_changes
    from ase.neighborlist import neighbor_list
    from ase.units import kcal, mol
except ImportError as exc:
    raise ImportError("Install xreac[ase] to use the ASE calculator") from exc

from .calculator import Calculator
from .geometry import validate_expansion_limit

KCAL_MOL_TO_EV = kcal / mol


class ReaxFFCalculator(ASECalculator):
    """ReaxFF for neutral ASE Atoms, with optional fixed-cell periodicity.

    Energies and forces use ASE units (eV, eV/A); charges and dipoles use e
    and e A. ``evaluation`` retains the last core result in kcal/mol units.
    Fixed-charge forces are the default; full_derivative=True is an explicit
    option for single-point charge-response calculations.
    neighbor_backend="ase" uses ASE's image-resolved neighbor list by default.
    Set it to "replicated" to use the native dense image/replication method.
    """

    implemented_properties = ["energy", "forces", "charges", "dipole"]
    default_parameters = {
        "full_derivative": False,
        "total_charge": 0,
        "max_expanded_atoms": 512,
        "neighbor_backend": "ase",
    }

    def __init__(
        self,
        force_field,
        *,
        full_derivative=False,
        total_charge=0,
        max_expanded_atoms=512,
        neighbor_backend="ase",
        **kwargs,
    ):
        self.core = Calculator(force_field, max_expanded_atoms=max_expanded_atoms)
        self.evaluation = None
        super().__init__(
            full_derivative=full_derivative,
            total_charge=total_charge,
            max_expanded_atoms=max_expanded_atoms,
            neighbor_backend=neighbor_backend,
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
        changed = super().set(**kwargs)
        if changed:
            self.core.max_expanded_atoms = self.parameters.max_expanded_atoms
            self.reset()
        return changed

    def reset(self):
        super().reset()
        self.evaluation = None

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
            cutoff = np.nextafter(float(self.core.force_field.general[12]), np.inf)
            i, j, S = neighbor_list("ijS", self.atoms, cutoff, self_interaction=False)
            neighbors = (i, j, S)
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
