from pathlib import Path

import numpy as np
import pytest
from scipy.constants import elementary_charge, epsilon_0
from scipy.integrate import quad
from scipy.optimize import brentq

from charge_injection_sim import load_config
from charge_injection_sim.config import InputConfig, ResolvedInputsSI
from charge_injection_sim.solver import (
    SolverError,
    SolverProgressEvent,
    solve_all_cases,
    solve_case,
)

BENCHMARK_PATH = Path(__file__).parents[1] / "configs" / "figure4b.toml"


def _synthetic_inputs(
    *,
    length_m: float,
    voltage_v: float,
    electron_mobility: float = 1e-8,
    hole_mobility: float = 3e-8,
    recombination: float = 0.0,
    tolerance: float = 1e-7,
    initial_nodes: int = 101,
) -> ResolvedInputsSI:
    inputs = load_config(BENCHMARK_PATH).to_si()
    material = inputs.material.model_copy(
        update={
            "relative_permittivity": 10.0,
            "oxygen_vacancy_m3": 3.120754537230381e22,
            "electron_mobility_m2_per_v_s": electron_mobility,
            "hole_mobility_m2_per_v_s": hole_mobility,
            "oxygen_vacancy_mobility_m2_per_v_s": 1e-10,
            "recombination_m3_per_s": recombination,
        }
    )
    experiment = inputs.experiment.model_copy(
        update={"voltage_v": voltage_v, "thickness_m": length_m}
    )
    settings = inputs.solver.model_copy(
        update={
            "residual_tolerance": tolerance,
            "initial_nodes": initial_nodes,
            "maximum_nodes": 10_000,
        }
    )
    return inputs.model_copy(
        update={"material": material, "experiment": experiment, "solver": settings}
    )


def test_charge_neutral_ohmic_limit() -> None:
    carrier_density = 1e20
    inputs = _synthetic_inputs(length_m=1e-6, voltage_v=0.4)

    result = solve_case(
        inputs,
        inputs.cases[0],
        electron_contact_m3=carrier_density,
        hole_contact_m3=carrier_density,
    )

    expected_field = inputs.experiment.voltage_v / inputs.experiment.thickness_m
    expected_sigma = 1e-6 + elementary_charge * (1e-8 + 3e-8) * carrier_density
    np.testing.assert_allclose(result.electric_field_v_per_m, expected_field, rtol=2e-7)
    np.testing.assert_allclose(result.electron_m3, carrier_density, rtol=2e-7)
    np.testing.assert_allclose(result.hole_m3, carrier_density, rtol=2e-7)
    assert result.diagnostics.current_density_a_per_m2 == pytest.approx(
        expected_sigma * expected_field, rel=2e-7
    )
    assert np.trapezoid(result.electric_field_v_per_m, result.position_m) == pytest.approx(
        inputs.experiment.voltage_v, rel=2e-7
    )


def test_quenched_equilibrium_ohmic_limit() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["background"]["model"] = "quenched_equilibrium"
    inputs = InputConfig.model_validate(data).to_si()
    background = inputs.background

    result = solve_case(
        inputs,
        inputs.cases[0],
        electron_contact_m3=background.electron_m3,
        hole_contact_m3=background.hole_m3,
    )

    expected_field = inputs.experiment.voltage_v / inputs.experiment.thickness_m
    expected_sigma = (
        2.0
        * elementary_charge
        * inputs.material.oxygen_vacancy_mobility_m2_per_v_s
        * inputs.material.oxygen_vacancy_m3
        + elementary_charge * inputs.material.electron_mobility_m2_per_v_s * background.electron_m3
        + elementary_charge * inputs.material.hole_mobility_m2_per_v_s * background.hole_m3
    )
    np.testing.assert_allclose(result.electric_field_v_per_m, expected_field, rtol=2e-7)
    np.testing.assert_allclose(result.electron_m3, background.electron_m3, rtol=2e-7)
    np.testing.assert_allclose(result.hole_m3, background.hole_m3, rtol=2e-7)
    assert result.diagnostics.current_density_a_per_m2 == pytest.approx(
        expected_sigma * expected_field, rel=2e-7
    )


def test_low_fe_n_type_background_converges_without_carrier_cancellation() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["background"]["model"] = "quenched_equilibrium"
    data["material"]["reported_total_fe_cm3"] = 1e17
    inputs = InputConfig.model_validate(data).to_si()

    result = solve_case(inputs, inputs.cases[0])

    assert inputs.background.electron_m3 / 1e6 == pytest.approx(4.76e18, rel=1e-12)
    assert result.diagnostics.minimum_electron_m3 > 0
    assert result.diagnostics.minimum_hole_m3 > 0
    assert result.diagnostics.poisson_scaled_residual <= 1e-4
    assert result.diagnostics.continuity_scaled_residual <= 1e-4


def test_recombination_free_bipolar_quadrature_reference() -> None:
    # Choose j, the Eq. (5) invariant, and endpoint fields, then derive L, V,
    # and contact densities independently by quadrature.
    current = 1.0
    invariant = 0.05
    vacancy_sigma = 1e-6
    electron_mobility = 1e-8
    hole_mobility = 3e-8
    permittivity = 10.0 * epsilon_0
    anode_field = 2e5
    cathode_field = 8e5

    def electron_current(field: float) -> float:
        return invariant * np.exp(-vacancy_sigma * field / current)

    def electron_density(field: float) -> float:
        return electron_current(field) / (elementary_charge * electron_mobility * field)

    def hole_density(field: float) -> float:
        return (current - electron_current(field) - vacancy_sigma * field) / (
            elementary_charge * hole_mobility * field
        )

    def field_gradient(field: float) -> float:
        return elementary_charge * (hole_density(field) - electron_density(field)) / permittivity

    length = quad(lambda field: 1.0 / field_gradient(field), anode_field, cathode_field)[0]
    voltage = quad(lambda field: field / field_gradient(field), anode_field, cathode_field)[0]
    inputs = _synthetic_inputs(
        length_m=length,
        voltage_v=voltage,
        electron_mobility=electron_mobility,
        hole_mobility=hole_mobility,
    )
    result = solve_case(
        inputs,
        inputs.cases[0],
        electron_contact_m3=electron_density(cathode_field),
        hole_contact_m3=hole_density(anode_field),
    )

    sample_positions = np.linspace(0.0, length, 21)

    def reference_field(position: float) -> float:
        if position == 0.0:
            return anode_field
        if position == length:
            return cathode_field
        return brentq(
            lambda field: (
                quad(lambda value: 1.0 / field_gradient(value), anode_field, field)[0] - position
            ),
            anode_field,
            cathode_field,
        )

    expected_field = np.array([reference_field(position) for position in sample_positions])
    actual_field = np.interp(sample_positions, result.position_m, result.electric_field_v_per_m)
    actual_electron = np.interp(sample_positions, result.position_m, result.electron_m3)
    actual_hole = np.interp(sample_positions, result.position_m, result.hole_m3)
    np.testing.assert_allclose(actual_field, expected_field, rtol=3e-5)
    np.testing.assert_allclose(
        actual_electron,
        [electron_density(field) for field in expected_field],
        rtol=4e-5,
    )
    np.testing.assert_allclose(
        actual_hole,
        [hole_density(field) for field in expected_field],
        rtol=4e-5,
    )
    assert result.diagnostics.current_density_a_per_m2 == pytest.approx(current, rel=2e-6)
    invariant_values = (
        result.electron_conductivity_s_per_m
        * result.electric_field_v_per_m
        * np.exp(
            result.vacancy_conductivity_s_per_m
            * result.electric_field_v_per_m
            / result.diagnostics.current_density_a_per_m2
        )
    )
    assert np.ptp(invariant_values) / invariant == pytest.approx(0.0, abs=2e-6)


@pytest.mark.regression
def test_benchmark_three_cases_converge_and_refine() -> None:
    inputs = load_config(BENCHMARK_PATH).to_si()
    results = solve_all_cases(inputs)

    assert [result.case_name for result in results] == [case.name for case in inputs.cases]
    assert [result.diagnostics.current_density_a_per_m2 for result in results] == pytest.approx(
        [0.19535093, 1.32160543, 2.39738471], rel=2e-5
    )
    assert all(result.diagnostics.status == 0 for result in results)
    assert all(result.diagnostics.maximum_boundary_residual <= 1e-5 for result in results)
    assert all(result.diagnostics.voltage_relative_error <= 1e-5 for result in results)
    assert all(result.diagnostics.poisson_scaled_residual <= 1e-4 for result in results)
    assert all(result.diagnostics.continuity_scaled_residual <= 1e-4 for result in results)
    assert all(np.all(np.diff(result.position_m) > 0) for result in results)

    refined_settings = inputs.solver.model_copy(
        update={"residual_tolerance": 3e-6, "initial_nodes": 301}
    )
    refined = solve_all_cases(inputs.model_copy(update={"solver": refined_settings}))
    for baseline, tighter in zip(results, refined, strict=True):
        common_position = np.linspace(0.0, inputs.experiment.thickness_m, 301)
        baseline_sigma = np.interp(
            common_position, baseline.position_m, baseline.total_conductivity_s_per_m
        )
        tighter_sigma = np.interp(
            common_position, tighter.position_m, tighter.total_conductivity_s_per_m
        )
        assert tighter.diagnostics.current_density_a_per_m2 == pytest.approx(
            baseline.diagnostics.current_density_a_per_m2, rel=5e-3
        )
        assert np.max(np.abs(tighter_sigma / baseline_sigma - 1.0)) <= 5e-3


def test_solver_reports_node_limit_failure() -> None:
    inputs = load_config(BENCHMARK_PATH).to_si()
    settings = inputs.solver.model_copy(
        update={"residual_tolerance": 1e-10, "initial_nodes": 2, "maximum_nodes": 2}
    )

    with pytest.raises(SolverError, match="maximum number of mesh nodes"):
        solve_case(inputs.model_copy(update={"solver": settings}), inputs.cases[-1])


@pytest.mark.regression
def test_solve_all_cases_emits_progress_events() -> None:
    inputs = load_config(BENCHMARK_PATH).to_si()
    events: list[SolverProgressEvent] = []

    results = solve_all_cases(inputs, progress=events.append)

    assert results
    total_cases = len(inputs.cases)
    continuation_steps = 5 if len(events) > total_cases else 1
    assert len(events) >= total_cases
    assert all(event.total_cases == total_cases for event in events)
    assert all(event.case_steps_total in (1, continuation_steps) for event in events)
    assert events[0].completed_cases == 0
    assert events[-1].completed_cases == total_cases - 1
    assert all(event.case_elapsed_seconds >= 0.0 for event in events)
    assert all(event.node_count >= inputs.solver.initial_nodes for event in events)
