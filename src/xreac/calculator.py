"""Public calculator and force-based geometry relaxation."""

from dataclasses import dataclass
from numbers import Real

import autograd.numpy as anp
from autograd import make_vjp
import numpy as np

from .energy import COMPONENTS, EnergyModel
from .geometry import Boundary, validate_expansion_limit
from .neighbors import replicated_neighbors


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
    cell_repetitions: tuple[int, int, int] = (1, 1, 1)
    neighbor_backend: str = "replicated"

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


def validate_total_charge(total_charge):
    """Validate a net charge in elementary-charge units, including fractions."""
    if (
        isinstance(total_charge, (bool, np.bool_))
        or not isinstance(total_charge, Real)
        or not np.isfinite(total_charge)
    ):
        raise ValueError("total_charge must be a finite real scalar")


def validate_input(symbols, positions, total_charge=0, cell=None, pbc=None):
    symbols = tuple(symbols)
    if not symbols or any(not isinstance(s, str) or not s for s in symbols):
        raise ValueError("A nonempty sequence of parameter-file atom labels is required")
    validate_total_charge(total_charge)
    x = np.array(positions, dtype=np.float64, copy=True)
    if x.shape != (len(symbols), 3) or not np.isfinite(x).all():
        raise ValueError("positions must be a finite (N, 3) array in Angstrom")
    distances = np.linalg.norm(Boundary(cell, pbc).minimum_displacements(x), axis=-1)
    if np.any(distances[np.triu_indices(len(x), 1)] < 1e-6):
        raise ValueError("Coincident atoms are not supported")
    return symbols, x


class Calculator:
    def __init__(self, force_field, *, max_expanded_atoms=512):
        force_field.validate_model()
        validate_expansion_limit(max_expanded_atoms)
        self.force_field = force_field
        self.max_expanded_atoms = max_expanded_atoms

    def evaluate(
        self, symbols, positions, *, total_charge=0, cell=None, pbc=None, full_derivative=False, neighbors=None
    ):
        """Return energy, equilibrated charges, and the selected forces.

        By default, fixed-charge forces hold the freshly equilibrated charges
        constant during differentiation, matching LAMMPS. Set full_derivative
        to True for charge-response forces: the full reported energy derivative
        through QEq. Only the selected derivative is evaluated.

        total_charge is a finite real scalar in elementary-charge units and
        constrains the sum of charges in the input cell (default zero).
        Fractional charges are allowed. Periodic electrostatics use the existing
        finite cutoff, without a compensating background or Ewald summation.

        cell contains three lengths or three row vectors in Angstrom. Supplying
        it enables all periodic directions unless pbc is explicitly set to a
        boolean or three flags. Pass neighbors=(i, j, S) for a full directed
        image list built by the caller, e.g. ASE neighbor_list("ijS", atoms,
        cutoff). Each edge has vector x[j]-x[i]+S@cell. Include both directions
        and all images within the nonbonded cutoff; exclude zero-shift self
        edges. Extra neighbors beyond the cutoff (a skin) are allowed. The
        caller is responsible for rebuilding the list when needed. No neighbor
        search or ASE import occurs when arrays are supplied. Autograd treats
        i, j, S as constants and differentiates the vectors through x only.
        Without neighbors, use native dense image search and small-cell replication.
        Energies and properties refer to the input cell. cell_repetitions records
        neighbor-search replication. Both paths use the same input-cell energy
        model. bond_orders sums over images; bond_counts counts each image.
        Dipoles use the force-field center of mass as origin and the supplied
        coordinate branch; they change upon individual atom wrapping.
        """
        if not isinstance(full_derivative, bool):
            raise ValueError("full_derivative must be a boolean")
        symbols, x = validate_input(symbols, positions, total_charge, cell, pbc)
        self.force_field.validate_model(symbols)
        if neighbors is not None:
            repetitions = (1, 1, 1)
            neighbor_backend = "provided"
        else:
            neighbors, repetitions = replicated_neighbors(
                x, self.force_field.general[12], cell, pbc, max_expanded_atoms=self.max_expanded_atoms
            )
            neighbor_backend = "replicated"
        model = EnergyModel(self.force_field, symbols, neighbors, cell, pbc, total_charge=total_charge)
        try:
            # Solve QEq outside the trace for fixed-charge forces. The traced
            # energy pass supplies both values and forces, without replaying
            # all energy terms merely to obtain the reported components.
            fixed_charges = None
            if not full_derivative:
                _, distances = model.edges.geometry(x)
                fixed_charges = model.electrostatics(distances)[0]

            bond_state = None

            def values(y):
                nonlocal bond_state
                components, charges, bond_state = model.components(y, fixed_charges, return_bond_state=True)
                return anp.concatenate((components, charges))

            backward, values_at_x = make_vjp(values)(x)
            components, charges = values_at_x[: len(COMPONENTS)], values_at_x[len(COMPONENTS) :]
            force = -backward(np.concatenate((np.ones(len(COMPONENTS)), np.zeros(len(x)))))
        except np.linalg.LinAlgError as exc:
            raise ValueError("Singular QEq system; check the geometry and parameters") from exc
        if not all(np.isfinite(v).all() for v in (components, charges, force)):
            raise ValueError("Non-finite energy, charges, or forces; check the geometry")
        if abs(charges.sum() - total_charge) > 1e-8:
            raise ValueError("QEq charge constraint failed")
        properties = model.properties(x, charges, bond_state)
        return Evaluation(
            float(components.sum()),
            force,
            charges,
            dict(zip(COMPONENTS, map(float, components))),
            full_derivative,
            cell_repetitions=repetitions,
            neighbor_backend=neighbor_backend,
            **properties,
        )

    def relax(
        self,
        symbols,
        positions,
        *,
        force_tolerance=1e-4,
        max_iterations=500,
        total_charge=0,
        cell=None,
        pbc=None,
        backend="ase",
    ):
        """Relax with fixed-charge forces, using ASE FIRE by default.

        QEq is solved at each geometry, without differentiation through the
        charge solve. Convergence requires the largest absolute Cartesian force
        component to be at most force_tolerance (kcal/mol/A). FIRE uses damped
        fictitious dynamics, not an energy line search or physical time evolution.
        Set backend="native" to use the original NumPy FIRE implementation
        with native replication, without ASE, including for benchmarks. The
        default ASE FIRE path uses ASE neighbor lists. cell and pbc follow evaluate();
        the cell is fixed during relaxation.
        """
        symbols, x = validate_input(symbols, positions, total_charge, cell, pbc)
        boundary = Boundary(cell, pbc)
        if not np.isfinite(force_tolerance) or force_tolerance <= 0:
            raise ValueError("force_tolerance must be positive and finite")
        if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
            raise ValueError("max_iterations must be a positive integer")
        if backend not in ("ase", "native"):
            raise ValueError("backend must be 'ase' or 'native'")
        if backend == "ase":
            from .ase import relax_with_ase

            x, result, converged, iterations = relax_with_ase(
                self,
                symbols,
                x,
                force_tolerance,
                max_iterations,
                cell=cell,
                pbc=boundary.pbc,
                total_charge=total_charge,
            )
            message = "Fixed-charge force tolerance reached" if converged else "Maximum relaxation iterations reached"
            return Relaxation(x, result, converged, iterations, message)
        # FIRE (Bitzek et al., Phys. Rev. Lett. 97, 170201, 2006).
        # Unit fictitious masses; step parameters are optimizer scales, not fs.
        velocity = np.zeros_like(x)
        dt, dt_max, alpha = 0.02, 0.2, 0.1
        positive_steps = 0
        result = self.evaluate(symbols, x, total_charge=total_charge, cell=cell, pbc=pbc, full_derivative=False)
        for iteration in range(max_iterations + 1):
            force = result.forces
            if np.max(np.abs(force)) <= force_tolerance:
                return Relaxation(x, result, True, iteration, "Fixed-charge force tolerance reached")
            if iteration == max_iterations:
                break
            power = float(np.sum(velocity * force))
            if power > 0:
                velocity = (1 - alpha) * velocity + alpha * np.linalg.norm(velocity) / np.linalg.norm(force) * force
                positive_steps += 1
                if positive_steps > 5:
                    dt = min(dt * 1.1, dt_max)
                    alpha *= 0.99
            else:
                velocity.fill(0.0)
                positive_steps = 0
                alpha = 0.1
                if iteration > 0:
                    dt *= 0.5
            velocity += dt * force
            displacement = dt * velocity
            # Cap the largest atomic displacement at 0.1 A per iteration.
            largest_step = np.max(np.linalg.norm(displacement, axis=1))
            if largest_step > 0.1:
                displacement *= 0.1 / largest_step
            x = x + displacement
            result = self.evaluate(symbols, x, total_charge=total_charge, cell=cell, pbc=pbc, full_derivative=False)
        return Relaxation(x, result, False, max_iterations, "Maximum relaxation iterations reached")
