from pathlib import Path

import numpy as np
import pytest
from scipy.constants import electron_mass

from charge_injection_sim import load_config
from charge_injection_sim.physics import (
    conductivity_components_s_per_m,
    contact_density_m3,
    effective_density_of_states_m3,
    prepare_case,
    vacancy_conductivity_s_per_m,
)

BENCHMARK_PATH = Path(__file__).parents[1] / "configs" / "figure4b.toml"


def test_benchmark_density_of_states_and_contacts() -> None:
    inputs = load_config(BENCHMARK_PATH).to_si()
    temperature = inputs.experiment.temperature_k

    electron_dos = effective_density_of_states_m3(6.0 * electron_mass, temperature)
    hole_dos = effective_density_of_states_m3(12.0 * electron_mass, temperature)
    electron_contact = contact_density_m3(
        electron_dos, inputs.cases[0].electron_barrier_j, temperature
    )
    hole_contact = contact_density_m3(hole_dos, inputs.cases[0].hole_barrier_j, temperature)

    assert electron_dos == pytest.approx(7.538e26, rel=1e-4)
    assert hole_dos == pytest.approx(2.132e27, rel=1e-4)
    assert electron_contact == pytest.approx(2.794e16, rel=2e-4)
    assert hole_contact == pytest.approx(7.901e16, rel=2e-4)


def test_vacancy_baseline_and_conductivity_sum() -> None:
    inputs = load_config(BENCHMARK_PATH).to_si()
    material = inputs.material
    vacancy = vacancy_conductivity_s_per_m(
        material.oxygen_vacancy_m3, material.oxygen_vacancy_mobility_m2_per_v_s
    )
    electron = np.array([0.0, 1e20], dtype=np.float64)
    hole = np.array([0.0, 2e20], dtype=np.float64)

    electron_sigma, hole_sigma, vacancy_sigma, total = conductivity_components_s_per_m(
        electron,
        hole,
        material.electron_mobility_m2_per_v_s,
        material.hole_mobility_m2_per_v_s,
        vacancy,
    )

    assert vacancy == pytest.approx(1.713047e-6, rel=1e-6)
    assert vacancy / 100.0 == pytest.approx(1.713047e-8, rel=1e-6)
    np.testing.assert_allclose(total, electron_sigma + hole_sigma + vacancy_sigma)
    assert np.all(total > 0)


def test_case_scaling_has_unit_conductivity_fractions() -> None:
    inputs = load_config(BENCHMARK_PATH).to_si()
    physics = prepare_case(inputs, inputs.cases[0])

    assert physics.gamma_vacancy + physics.gamma_electron + physics.gamma_hole == pytest.approx(1.0)
    assert physics.current_scale_a_per_m2 == pytest.approx(
        physics.conductivity_scale_s_per_m * physics.electric_field_scale_v_per_m
    )
    assert physics.electron_contact_scaled > 0
    assert physics.hole_contact_scaled > 0
