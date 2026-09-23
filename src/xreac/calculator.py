"""Public calculator and force-based geometry relaxation."""
from dataclasses import dataclass

import autograd.numpy as anp
from autograd import grad
import numpy as np

from .energy import COMPONENTS, EnergyModel
from .geometry import Boundary


@dataclass(frozen=True)
class Evaluation:
    energy: float
    forces: np.ndarray
    charges: np.ndarray
    components: dict[str, float]
    full_derivative: bool
    bond_orders: np.ndarray
    total_bond_orders: np.ndarray
    lone_pairs: np.ndarray
    bond_counts: np.ndarray
    dipole: np.ndarray

    @property
    def force_convention(self):
        """Name of the derivative used for ``forces``."""
        return "charge_response" if self.full_derivative else "fixed_charge"


@dataclass(frozen=True)
class Relaxation:
    positions: np.ndarray
    evaluation: Evaluation
    converged: bool
    iterations: int
    message: str


def validate_input(symbols, positions, total_charge=0, cell=None, pbc=None):
    symbols = tuple(symbols)
    if not symbols or any(not isinstance(s, str) or not s for s in symbols):
        raise ValueError("A nonempty sequence of parameter-file atom labels is required")
    if total_charge != 0:
        raise ValueError("Only neutral systems are validated")
    x = np.array(positions, dtype=np.float64, copy=True)
    if x.shape != (len(symbols), 3) or not np.isfinite(x).all():
        raise ValueError("positions must be a finite (N, 3) array in Angstrom")
    distances = np.linalg.norm(Boundary(cell, pbc).minimum_displacements(x), axis=-1)
    if np.any(distances[np.triu_indices(len(x), 1)] < 1e-6):
        raise ValueError("Coincident atoms are not supported")
    return symbols, x


class Calculator:
    def __init__(self, force_field):
        force_field.validate_model()
        self.force_field = force_field

    def evaluate(self, symbols, positions, *, total_charge=0, cell=None,
                 pbc=None, full_derivative=False):
        """Return energy, equilibrated charges, and the selected forces.

        By default, fixed-charge forces hold the freshly equilibrated charges
        constant during differentiation, matching LAMMPS. Set full_derivative
        to True for charge-response forces: the full reported energy derivative
        through QEq. Only the selected derivative is evaluated.

        cell contains three lengths or three row vectors in Angstrom. Supplying
        it enables all periodic directions unless pbc is explicitly set to a
        boolean or three flags. Periodic cell heights must exceed both the
        nonbonded cutoff and twice the bond cutoff (10 A for bundled files).
        Dipoles use the supplied coordinate branch and change upon wrapping.
        """
        if not isinstance(full_derivative, bool):
            raise ValueError("full_derivative must be a boolean")
        symbols, x = validate_input(symbols, positions, total_charge, cell, pbc)
        model = EnergyModel(self.force_field, symbols, cell, pbc)
        try:
            components, charges = model.components(x)
            fixed_charges = None if full_derivative else charges
            force = -grad(lambda y: anp.sum(model.components(y, fixed_charges)[0]))(x)
        except np.linalg.LinAlgError as exc:
            raise ValueError("Singular QEq system; check the geometry and parameters") from exc
        if not all(np.isfinite(v).all() for v in (components, charges, force)):
            raise ValueError("Non-finite energy, charges, or forces; check the geometry")
        if abs(charges.sum()) > 1e-8:
            raise ValueError("QEq charge constraint failed")
        properties = model.properties(x, charges)
        return Evaluation(float(components.sum()), force, charges,
                          dict(zip(COMPONENTS, map(float, components))), full_derivative, **properties)

    def relax(self, symbols, positions, *, force_tolerance=1e-4, max_iterations=500,
              total_charge=0, cell=None, pbc=None, backend="ase"):
        """Relax with fixed-charge forces, using ASE FIRE by default.

        QEq is solved at each geometry, without differentiation through the
        charge solve. Convergence requires the largest absolute Cartesian force
        component to be at most force_tolerance (kcal/mol/A). FIRE uses damped
        fictitious dynamics, not an energy line search or physical time evolution.
        Set backend="native" to use the original NumPy FIRE implementation
        without ASE, including for benchmarks. cell and pbc follow evaluate();
        the cell is fixed during relaxation.
        """
        symbols, x = validate_input(symbols, positions, total_charge, cell, pbc)
        boundary = Boundary(cell, pbc)
        boundary.validate_cutoff(self.force_field.general[12], min(5., self.force_field.general[12]))
        if not np.isfinite(force_tolerance) or force_tolerance <= 0:
            raise ValueError("force_tolerance must be positive and finite")
        if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
            raise ValueError("max_iterations must be a positive integer")
        if backend not in ("ase", "native"):
            raise ValueError("backend must be 'ase' or 'native'")
        if backend == "ase":
            from .ase import relax_with_ase

            x, result, converged, iterations = relax_with_ase(
                self, symbols, x, force_tolerance, max_iterations, cell=cell, pbc=boundary.pbc)
            message = "Fixed-charge force tolerance reached" if converged else "Maximum relaxation iterations reached"
            return Relaxation(x, result, converged, iterations, message)
        # FIRE (Bitzek et al., Phys. Rev. Lett. 97, 170201, 2006).
        # Unit fictitious masses; step parameters are optimizer scales, not fs.
        velocity = np.zeros_like(x)
        dt, dt_max, alpha = .02, .2, .1
        positive_steps = 0
        result = self.evaluate(symbols, x, cell=cell, pbc=pbc, full_derivative=False)
        for iteration in range(max_iterations+1):
            force = result.forces
            if np.max(np.abs(force)) <= force_tolerance:
                return Relaxation(x, result, True, iteration, "Fixed-charge force tolerance reached")
            if iteration == max_iterations:
                break
            power = float(np.sum(velocity*force))
            if power > 0:
                velocity = (1-alpha)*velocity + alpha*np.linalg.norm(velocity)/np.linalg.norm(force)*force
                positive_steps += 1
                if positive_steps > 5:
                    dt = min(dt*1.1, dt_max)
                    alpha *= .99
            else:
                velocity.fill(0.)
                positive_steps = 0
                alpha = .1
                if iteration > 0:
                    dt *= .5
            velocity += dt*force
            displacement = dt*velocity
            # Cap the largest atomic displacement at 0.1 A per iteration.
            largest_step = np.max(np.linalg.norm(displacement, axis=1))
            if largest_step > .1:
                displacement *= .1/largest_step
            x = x + displacement
            result = self.evaluate(symbols, x, cell=cell, pbc=pbc, full_derivative=False)
        return Relaxation(x, result, False, max_iterations, "Maximum relaxation iterations reached")
