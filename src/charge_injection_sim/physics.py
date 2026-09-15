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

    @property
    def electron_contact_scaled(self) -> float:
        return self.electron_contact_m3 / self.electron_scale_m3

    @property
    def hole_contact_scaled(self) -> float:
        return self.hole_contact_m3 / self.hole_scale_m3


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
        contact_density_m3(electron_dos, case.electron_barrier_j, experiment.temperature_k)
        if electron_contact_m3 is None
        else electron_contact_m3
    )
    hole_contact = (
        contact_density_m3(hole_dos, case.hole_barrier_j, experiment.temperature_k)
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
        vacancy_conductivity / (elementary_charge * material.electron_mobility_m2_per_v_s),
    )
    hole_scale = max(
        hole_contact,
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


def scaled_rhs(
    position_scaled: FloatArray,
    state_scaled: FloatArray,
    parameters: FloatArray,
    physics: CasePhysics,
    homotopy: float = 1.0,
) -> FloatArray:
    """Evaluate the dimensionless BVP equations."""
    del position_scaled
    electric_field, electron, _voltage = state_scaled
    log_current = float(parameters[0])
    current = np.exp(log_current)
    hole = reconstruct_holes_scaled(electric_field, electron, log_current, physics)
    field_derivative = homotopy * (physics.lambda_hole * hole - physics.lambda_electron * electron)
    electron_derivative = (
        homotopy * physics.rho_hole * electron * hole / electric_field
        - (electron / electric_field)
        * (1.0 + physics.gamma_vacancy * electric_field / current)
        * field_derivative
    )
    return np.vstack((field_derivative, electron_derivative, electric_field))


def scaled_boundary_residuals(
    state_at_anode: FloatArray,
    state_at_cathode: FloatArray,
    parameters: FloatArray,
    physics: CasePhysics,
) -> FloatArray:
    """Evaluate normalized contact and accumulated-voltage residuals."""
    holes_at_anode = reconstruct_holes_scaled(
        state_at_anode[0:1], state_at_anode[1:2], float(parameters[0]), physics
    )[0]
    hole_normalizer = physics.hole_contact_scaled if physics.hole_contact_scaled > 0.0 else 1.0
    electron_normalizer = (
        physics.electron_contact_scaled if physics.electron_contact_scaled > 0.0 else 1.0
    )
    return np.array(
        [
            (holes_at_anode - physics.hole_contact_scaled) / hole_normalizer,
            (state_at_cathode[1] - physics.electron_contact_scaled) / electron_normalizer,
            state_at_anode[2],
            state_at_cathode[2] - 1.0,
        ],
        dtype=np.float64,
    )
