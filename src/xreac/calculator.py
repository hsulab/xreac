"""Public calculator and optional SciPy geometry relaxation."""
from dataclasses import dataclass

import autograd.numpy as anp
from autograd import grad, value_and_grad
import numpy as np

from .energy import COMPONENTS, EnergyModel


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


def validate_input(symbols, positions, total_charge=0, cell=None):
    symbols = tuple(symbols)
    if not symbols or any(not isinstance(s, str) or not s for s in symbols):
        raise ValueError("A nonempty sequence of parameter-file atom labels is required")
    if total_charge != 0:
        raise ValueError("Only neutral systems are validated")
    if cell is not None:
        raise ValueError("Periodic cells are not supported")
    x = np.array(positions, dtype=np.float64, copy=True)
    if x.shape != (len(symbols), 3) or not np.isfinite(x).all():
        raise ValueError("positions must be a finite (N, 3) array in Angstrom")
    distances = np.linalg.norm(x[:, None]-x[None, :], axis=-1)
    if np.any(distances[np.triu_indices(len(x), 1)] < 1e-6):
        raise ValueError("Coincident atoms are not supported")
    return symbols, x


class Calculator:
    def __init__(self, force_field):
        force_field.validate_model()
        self.force_field = force_field

    def evaluate(self, symbols, positions, *, total_charge=0, cell=None,
                 full_derivative=False):
        """Return energy, equilibrated charges, and the selected forces.

        By default, fixed-charge forces hold the freshly equilibrated charges
        constant during differentiation, matching LAMMPS. Set full_derivative
        to True for charge-response forces: the full reported energy derivative
        through QEq. Only the selected derivative is evaluated.
        """
        if not isinstance(full_derivative, bool):
            raise ValueError("full_derivative must be a boolean")
        symbols, x = validate_input(symbols, positions, total_charge, cell)
        model = EnergyModel(self.force_field, symbols)
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
              total_charge=0, cell=None):
        """Minimize energy with charge-response forces (full_derivative=True).

        Unlike evaluate()'s default, energy minimization requires the full
        derivative. The returned evaluation records this mode. Force tolerance
        is in kcal/mol/A.
        """
        try:
            from scipy.optimize import minimize
        except ImportError as exc:
            raise ImportError("Install xreac[relax] to use geometry relaxation") from exc
        symbols, x = validate_input(symbols, positions, total_charge, cell)
        if not np.isfinite(force_tolerance) or force_tolerance <= 0:
            raise ValueError("force_tolerance must be positive and finite")
        if not isinstance(max_iterations, int) or max_iterations < 1:
            raise ValueError("max_iterations must be a positive integer")
        model = EnergyModel(self.force_field, symbols)
        vg = value_and_grad(lambda flat: anp.sum(model.components(flat.reshape(x.shape))[0]))

        def objective(flat):
            validate_input(symbols, flat.reshape(x.shape))
            value, derivative = vg(flat)
            if not np.isfinite(value) or not np.isfinite(derivative).all():
                raise ValueError("Non-finite energy or gradient during relaxation")
            return float(value), np.asarray(derivative)

        opt = minimize(objective, x.ravel(), method="L-BFGS-B", jac=True,
                       options={"gtol": force_tolerance, "ftol": 0.0,
                                "maxiter": max_iterations, "maxls": 40})
        positions = opt.x.reshape(x.shape)
        result = self.evaluate(symbols, positions, full_derivative=True)
        converged = bool(np.max(np.abs(result.forces)) <= force_tolerance)
        return Relaxation(positions, result, converged, int(opt.nit), str(opt.message))
