"""Physical relations and dimensionless equations for the drift-only model."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.constants import Boltzmann, elementary_charge, epsilon_0, h

from charge_injection_sim.config import BarrierCaseSI, ResolvedInputsSI

FloatArray = NDArray[np.float64]


def effective_density_of_states_m3(effective_mass_kg: float, temperature_k: float) -> float:
    """Return the parabolic-band effective density of states in m^-3."""
    return float(2.0 * (2.0 * np.pi * effective_mass_kg * Boltzmann * temperature_k / h**2) ** 1.5)


def contact_density_m3(
    density_of_states_m3: float, barrier_j: float, temperature_k: float
) -> float:
    """Return the thermally activated carrier density at a contact."""
    return float(density_of_states_m3 * np.exp(-barrier_j / (Boltzmann * temperature_k)))


def vacancy_conductivity_s_per_m(
    oxygen_vacancy_m3: float, oxygen_vacancy_mobility_m2_per_v_s: float
) -> float:
    """Return the doubly charged oxygen-vacancy conductivity."""
    return float(2.0 * elementary_charge * oxygen_vacancy_mobility_m2_per_v_s * oxygen_vacancy_m3)


def conductivity_components_s_per_m(
    electron_m3: FloatArray,
    hole_m3: FloatArray,
    electron_mobility_m2_per_v_s: float,
    hole_mobility_m2_per_v_s: float,
    vacancy_conductivity: float,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return electron, hole, vacancy, and total conductivity arrays."""
    electron = elementary_charge * electron_mobility_m2_per_v_s * electron_m3
    hole = elementary_charge * hole_mobility_m2_per_v_s * hole_m3
    vacancy = np.full_like(electron, vacancy_conductivity)
    return electron, hole, vacancy, electron + hole + vacancy


@dataclass(frozen=True)
class CasePhysics:
    """Precomputed dimensional constants and dimensionless scales for one case."""

    name: str
    length_m: float
    voltage_v: float
    electric_field_scale_v_per_m: float
    permittivity_f_per_m: float
    electron_mobility_m2_per_v_s: float
    hole_mobility_m2_per_v_s: float
    recombination_m3_per_s: float
    vacancy_conductivity_s_per_m: float
    electron_contact_m3: float
    hole_contact_m3: float
    equilibrium_electron_m3: float
    equilibrium_hole_m3: float
    electron_scale_m3: float
    hole_scale_m3: float
    conductivity_scale_s_per_m: float
    current_scale_a_per_m2: float
    gamma_vacancy: float
    gamma_electron: float
    gamma_hole: float
    lambda_electron: float
    lambda_hole: float
    rho_hole: float
    rho_electron: float
    equilibrium_product_scaled: float
    solve_for_electron: bool

    @property
    def electron_contact_scaled(self) -> float:
        return self.electron_contact_m3 / self.electron_scale_m3

    @property
    def hole_contact_scaled(self) -> float:
        return self.hole_contact_m3 / self.hole_scale_m3

    @property
    def equilibrium_electron_scaled(self) -> float:
        return self.equilibrium_electron_m3 / self.electron_scale_m3

    @property
    def equilibrium_hole_scaled(self) -> float:
        return self.equilibrium_hole_m3 / self.hole_scale_m3


def prepare_case(
    inputs: ResolvedInputsSI,
    case: BarrierCaseSI,
    *,
    electron_contact_m3: float | None = None,
    hole_contact_m3: float | None = None,
) -> CasePhysics:
    """Precompute constants and well-conditioned scales for one barrier case."""
    material = inputs.material
    experiment = inputs.experiment
    electron_dos = effective_density_of_states_m3(
        material.electron_effective_mass_kg, experiment.temperature_k
    )
    hole_dos = effective_density_of_states_m3(
        material.hole_effective_mass_kg, experiment.temperature_k
    )
    electron_contact = (
        inputs.background.electron_m3
        + contact_density_m3(electron_dos, case.electron_barrier_j, experiment.temperature_k)
        if electron_contact_m3 is None
        else electron_contact_m3
    )
    hole_contact = (
        inputs.background.hole_m3
        + contact_density_m3(hole_dos, case.hole_barrier_j, experiment.temperature_k)
        if hole_contact_m3 is None
        else hole_contact_m3
    )
    if electron_contact < 0 or hole_contact < 0:
        raise ValueError("contact densities must be nonnegative")

    vacancy_conductivity = vacancy_conductivity_s_per_m(
        material.oxygen_vacancy_m3, material.oxygen_vacancy_mobility_m2_per_v_s
    )
    # Keep each carrier's scale finite even when a limiting case has a zero contact
    # density, and avoid reconstructing a small carrier by catastrophic subtraction.
    electron_scale = max(
        electron_contact,
        inputs.background.electron_m3,
        vacancy_conductivity / (elementary_charge * material.electron_mobility_m2_per_v_s),
    )
    hole_scale = max(
        hole_contact,
        inputs.background.hole_m3,
        vacancy_conductivity / (elementary_charge * material.hole_mobility_m2_per_v_s),
    )
    electron_conductivity_scale = (
        elementary_charge * material.electron_mobility_m2_per_v_s * electron_scale
    )
    hole_conductivity_scale = elementary_charge * material.hole_mobility_m2_per_v_s * hole_scale
    conductivity_scale = (
        vacancy_conductivity + electron_conductivity_scale + hole_conductivity_scale
    )
    electric_field_scale = experiment.voltage_v / experiment.thickness_m
    permittivity = epsilon_0 * material.relative_permittivity

    return CasePhysics(
        name=case.name,
        length_m=experiment.thickness_m,
        voltage_v=experiment.voltage_v,
        electric_field_scale_v_per_m=electric_field_scale,
        permittivity_f_per_m=permittivity,
        electron_mobility_m2_per_v_s=material.electron_mobility_m2_per_v_s,
        hole_mobility_m2_per_v_s=material.hole_mobility_m2_per_v_s,
        recombination_m3_per_s=material.recombination_m3_per_s,
        vacancy_conductivity_s_per_m=vacancy_conductivity,
        electron_contact_m3=electron_contact,
        hole_contact_m3=hole_contact,
        equilibrium_electron_m3=inputs.background.electron_m3,
        equilibrium_hole_m3=inputs.background.hole_m3,
        electron_scale_m3=electron_scale,
        hole_scale_m3=hole_scale,
        conductivity_scale_s_per_m=conductivity_scale,
        current_scale_a_per_m2=conductivity_scale * electric_field_scale,
        gamma_vacancy=vacancy_conductivity / conductivity_scale,
        gamma_electron=electron_conductivity_scale / conductivity_scale,
        gamma_hole=hole_conductivity_scale / conductivity_scale,
        lambda_electron=(
            elementary_charge
            * electron_scale
            * experiment.thickness_m
            / (permittivity * electric_field_scale)
        ),
        lambda_hole=(
            elementary_charge
            * hole_scale
            * experiment.thickness_m
            / (permittivity * electric_field_scale)
        ),
        rho_hole=(
            material.recombination_m3_per_s
            * hole_scale
            * experiment.thickness_m
            / (material.electron_mobility_m2_per_v_s * electric_field_scale)
        ),
        rho_electron=(
            material.recombination_m3_per_s
            * electron_scale
            * experiment.thickness_m
            / (material.hole_mobility_m2_per_v_s * electric_field_scale)
        ),
        equilibrium_product_scaled=(
            inputs.background.electron_m3
            * inputs.background.hole_m3
            / (electron_scale * hole_scale)
        ),
        solve_for_electron=(
            inputs.background.electron_m3 * material.electron_mobility_m2_per_v_s
            <= inputs.background.hole_m3 * material.hole_mobility_m2_per_v_s
        ),
    )


def reconstruct_holes_scaled(
    electric_field_scaled: FloatArray,
    electron_scaled: FloatArray,
    log_current_scaled: float,
    physics: CasePhysics,
) -> FloatArray:
    """Reconstruct scaled holes from the spatially constant total current."""
    current_scaled = np.exp(log_current_scaled)
    return (
        current_scaled / electric_field_scaled
        - physics.gamma_vacancy
        - physics.gamma_electron * electron_scaled
    ) / physics.gamma_hole


def reconstruct_electrons_scaled(
    electric_field_scaled: FloatArray,
    hole_scaled: FloatArray,
    log_current_scaled: float,
    physics: CasePhysics,
) -> FloatArray:
    """Reconstruct scaled electrons from the spatially constant total current."""
    current_scaled = np.exp(log_current_scaled)
    return (
        current_scaled / electric_field_scaled
        - physics.gamma_vacancy
        - physics.gamma_hole * hole_scaled
    ) / physics.gamma_electron


def reconstruct_carriers_scaled(
    electric_field_scaled: FloatArray,
    carrier_scaled: FloatArray,
    log_current_scaled: float,
    physics: CasePhysics,
) -> tuple[FloatArray, FloatArray]:
    """Return scaled electron and hole densities for either state formulation."""
    if physics.solve_for_electron:
        electron = carrier_scaled
        hole = reconstruct_holes_scaled(
            electric_field_scaled, electron, log_current_scaled, physics
        )
    else:
        hole = carrier_scaled
        electron = reconstruct_electrons_scaled(
            electric_field_scaled, hole, log_current_scaled, physics
        )
    return electron, hole


def scaled_rhs(
    position_scaled: FloatArray,
    state_scaled: FloatArray,
    parameters: FloatArray,
    physics: CasePhysics,
    homotopy: float = 1.0,
) -> FloatArray:
    """Evaluate the dimensionless BVP equations."""
    del position_scaled
    electric_field, carrier, _voltage = state_scaled
    log_current = float(parameters[0])
    current = np.exp(log_current)
    electron, hole = reconstruct_carriers_scaled(electric_field, carrier, log_current, physics)
    field_derivative = homotopy * (
        physics.lambda_hole * (hole - physics.equilibrium_hole_scaled)
        - physics.lambda_electron * (electron - physics.equilibrium_electron_scaled)
    )
    reaction = electron * hole - physics.equilibrium_product_scaled
    if physics.solve_for_electron:
        carrier_derivative = (
            homotopy * physics.rho_hole * reaction / electric_field
            - (electron / electric_field)
            * (1.0 + physics.gamma_vacancy * electric_field / current)
            * field_derivative
        )
    else:
        electron_current_fraction = physics.gamma_electron * electron * electric_field / current
        carrier_derivative = (
            -homotopy * physics.rho_electron * reaction / electric_field
            - (
                hole / electric_field
                + physics.gamma_vacancy
                / (physics.gamma_hole * electric_field)
                * (1.0 - electron_current_fraction)
            )
            * field_derivative
        )
    return np.vstack((field_derivative, carrier_derivative, electric_field))


def scaled_boundary_residuals(
    state_at_anode: FloatArray,
    state_at_cathode: FloatArray,
    parameters: FloatArray,
    physics: CasePhysics,
) -> FloatArray:
    """Evaluate normalized contact and accumulated-voltage residuals."""
    electrons_at_cathode, _holes_at_cathode = reconstruct_carriers_scaled(
        state_at_cathode[0:1], state_at_cathode[1:2], float(parameters[0]), physics
    )
    _electrons_at_anode, holes_at_anode = reconstruct_carriers_scaled(
        state_at_anode[0:1], state_at_anode[1:2], float(parameters[0]), physics
    )
    hole_normalizer = physics.hole_contact_scaled if physics.hole_contact_scaled > 0.0 else 1.0
    electron_normalizer = (
        physics.electron_contact_scaled if physics.electron_contact_scaled > 0.0 else 1.0
    )
    return np.array(
        [
            (holes_at_anode[0] - physics.hole_contact_scaled) / hole_normalizer,
            (electrons_at_cathode[0] - physics.electron_contact_scaled) / electron_normalizer,
            state_at_anode[2],
            state_at_cathode[2] - 1.0,
        ],
        dtype=np.float64,
    )
