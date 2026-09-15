"""Boundary-value solver and independent diagnostics for the injection model."""

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from numpy.typing import NDArray
from scipy.constants import elementary_charge
from scipy.integrate import solve_bvp

from charge_injection_sim.config import BarrierCaseSI, ResolvedInputsSI, SolverSettings
from charge_injection_sim.physics import (
    CasePhysics,
    conductivity_components_s_per_m,
    prepare_case,
    reconstruct_holes_scaled,
    scaled_boundary_residuals,
    scaled_rhs,
)

FloatArray = NDArray[np.float64]


class SolverError(RuntimeError):
    """Raised when a numerical result fails convergence or physical checks."""


@dataclass(frozen=True)
class SolverDiagnostics:
    """Compact convergence and physical-validity diagnostics."""

    status: int
    message: str
    current_density_a_per_m2: float
    node_count: int
    continuation_attempts: int
    elapsed_seconds: float
    maximum_boundary_residual: float
    voltage_relative_error: float
    current_relative_spread: float
    poisson_scaled_residual: float
    continuity_scaled_residual: float
    minimum_electric_field_v_per_m: float
    minimum_electron_m3: float
    minimum_hole_m3: float
    minimum_conductivity_s_per_m: float
    maximum_hole_cancellation_ratio: float


@dataclass(frozen=True)
class SimulationResult:
    """One accepted dimensional solution on the adaptive solver mesh."""

    case_name: str
    position_m: FloatArray
    electric_field_v_per_m: FloatArray
    electron_m3: FloatArray
    hole_m3: FloatArray
    electron_conductivity_s_per_m: FloatArray
    hole_conductivity_s_per_m: FloatArray
    vacancy_conductivity_s_per_m: FloatArray
    total_conductivity_s_per_m: FloatArray
    diagnostics: SolverDiagnostics


@dataclass
class _WarmStart:
    position_scaled: FloatArray
    state_scaled: FloatArray
    log_current_scaled: float
    physics: CasePhysics


def _initial_guess(mesh: FloatArray, physics: CasePhysics) -> tuple[FloatArray, FloatArray]:
    electron = np.full_like(mesh, physics.electron_contact_scaled)
    hole = physics.hole_contact_scaled
    current = (
        physics.gamma_vacancy
        + physics.gamma_electron * physics.electron_contact_scaled
        + physics.gamma_hole * hole
    )
    state = np.vstack((np.ones_like(mesh), electron, mesh))
    return state, np.array([np.log(current)], dtype=np.float64)


def _rescale_warm_start(
    mesh: FloatArray, physics: CasePhysics, warm_start: _WarmStart
) -> tuple[FloatArray, FloatArray]:
    old = warm_start.physics
    old_state = np.vstack(
        [np.interp(mesh, warm_start.position_scaled, row) for row in warm_start.state_scaled]
    )
    electric_field = old_state[0] * old.electric_field_scale_v_per_m
    electron = old_state[1] * old.electron_scale_m3
    voltage = old_state[2] * old.voltage_v
    current = np.exp(warm_start.log_current_scaled) * old.current_scale_a_per_m2
    state = np.vstack(
        (
            electric_field / physics.electric_field_scale_v_per_m,
            electron / physics.electron_scale_m3,
            voltage / physics.voltage_v,
        )
    )
    return state, np.array([np.log(current / physics.current_scale_a_per_m2)])


def _solve_scaled(
    physics: CasePhysics,
    settings: SolverSettings,
    warm_start: _WarmStart | None,
) -> tuple[object, int]:
    mesh = np.linspace(0.0, 1.0, settings.initial_nodes, dtype=np.float64)
    state, parameters = (
        _initial_guess(mesh, physics)
        if warm_start is None
        else _rescale_warm_start(mesh, physics, warm_start)
    )
    solution = None
    attempts = 0
    homotopy_steps = (0.0, 0.1, 0.3, 0.6, 1.0) if warm_start is None else (1.0,)
    for homotopy in homotopy_steps:
        attempts += 1
        solution = solve_bvp(
            lambda x, y, p, h=homotopy: scaled_rhs(x, y, p, physics, h),
            lambda ya, yb, p: scaled_boundary_residuals(ya, yb, p, physics),
            mesh,
            state,
            p=parameters,
            tol=settings.residual_tolerance,
            max_nodes=settings.maximum_nodes,
            verbose=0,
        )
        if not solution.success:
            raise SolverError(
                f"{physics.name}: continuation at {homotopy:g} failed: {solution.message}"
            )
        mesh = solution.x
        state = solution.y
        parameters = solution.p
    assert solution is not None
    return solution, attempts


def _diagnose(
    solution: object,
    physics: CasePhysics,
    elapsed_seconds: float,
    continuation_attempts: int,
) -> SolverDiagnostics:
    dense_x = np.linspace(0.0, 1.0, max(1001, 5 * solution.x.size), dtype=np.float64)
    scaled = solution.sol(dense_x)
    scaled_derivative = solution.sol(dense_x, 1)
    electric_field = scaled[0] * physics.electric_field_scale_v_per_m
    electron = scaled[1] * physics.electron_scale_m3
    hole_scaled = reconstruct_holes_scaled(scaled[0], scaled[1], float(solution.p[0]), physics)
    hole = hole_scaled * physics.hole_scale_m3
    current_density = np.exp(solution.p[0]) * physics.current_scale_a_per_m2
    electron_sigma, _hole_sigma, _vacancy_sigma, total_sigma = conductivity_components_s_per_m(
        electron,
        hole,
        physics.electron_mobility_m2_per_v_s,
        physics.hole_mobility_m2_per_v_s,
        physics.vacancy_conductivity_s_per_m,
    )
    local_current = total_sigma * electric_field
    voltage_integral = np.trapezoid(electric_field, dense_x * physics.length_m)

    field_derivative = (
        scaled_derivative[0] * physics.electric_field_scale_v_per_m / physics.length_m
    )
    electron_derivative = scaled_derivative[1] * physics.electron_scale_m3 / physics.length_m
    poisson_residual = field_derivative - elementary_charge * (hole - electron) / (
        physics.permittivity_f_per_m
    )
    electron_current = electron_sigma * electric_field
    electron_current_derivative = (
        elementary_charge
        * physics.electron_mobility_m2_per_v_s
        * (electron_derivative * electric_field + electron * field_derivative)
    )
    vacancy_current_derivative = physics.vacancy_conductivity_s_per_m * field_derivative
    continuity_residual = (
        electron_current_derivative
        + electron_current / current_density * vacancy_current_derivative
        - elementary_charge * physics.recombination_m3_per_s * electron * hole
    )
    boundary = scaled_boundary_residuals(solution.y[:, 0], solution.y[:, -1], solution.p, physics)
    cancellation = (
        np.abs(np.exp(solution.p[0]) / scaled[0])
        + physics.gamma_vacancy
        + np.abs(physics.gamma_electron * scaled[1])
    ) / np.maximum(np.abs(physics.gamma_hole * hole_scaled), np.finfo(np.float64).tiny)

    diagnostics = SolverDiagnostics(
        status=int(solution.status),
        message=str(solution.message),
        current_density_a_per_m2=float(current_density),
        node_count=int(solution.x.size),
        continuation_attempts=continuation_attempts,
        elapsed_seconds=elapsed_seconds,
        maximum_boundary_residual=float(np.max(np.abs(boundary))),
        voltage_relative_error=float(abs(voltage_integral / physics.voltage_v - 1.0)),
        current_relative_spread=float(np.ptp(local_current) / current_density),
        poisson_scaled_residual=float(
            np.max(np.abs(poisson_residual))
            * physics.length_m
            / physics.electric_field_scale_v_per_m
        ),
        continuity_scaled_residual=float(
            np.max(np.abs(continuity_residual)) * physics.length_m / current_density
        ),
        minimum_electric_field_v_per_m=float(np.min(electric_field)),
        minimum_electron_m3=float(np.min(electron)),
        minimum_hole_m3=float(np.min(hole)),
        minimum_conductivity_s_per_m=float(np.min(total_sigma)),
        maximum_hole_cancellation_ratio=float(np.max(cancellation)),
    )
    return diagnostics


def _dimensional_arrays(solution: object, physics: CasePhysics) -> dict[str, FloatArray]:
    position = solution.x * physics.length_m
    electric_field = solution.y[0] * physics.electric_field_scale_v_per_m
    electron = solution.y[1] * physics.electron_scale_m3
    hole = (
        reconstruct_holes_scaled(solution.y[0], solution.y[1], float(solution.p[0]), physics)
        * physics.hole_scale_m3
    )
    electron_sigma, hole_sigma, vacancy_sigma, total_sigma = conductivity_components_s_per_m(
        electron,
        hole,
        physics.electron_mobility_m2_per_v_s,
        physics.hole_mobility_m2_per_v_s,
        physics.vacancy_conductivity_s_per_m,
    )
    return {
        "position_m": position,
        "electric_field_v_per_m": electric_field,
        "electron_m3": electron,
        "hole_m3": hole,
        "electron_conductivity_s_per_m": electron_sigma,
        "hole_conductivity_s_per_m": hole_sigma,
        "vacancy_conductivity_s_per_m": vacancy_sigma,
        "total_conductivity_s_per_m": total_sigma,
    }


def _validate_diagnostics(diagnostics: SolverDiagnostics, tolerance: float) -> None:
    values = np.array(
        [
            diagnostics.current_density_a_per_m2,
            diagnostics.maximum_boundary_residual,
            diagnostics.voltage_relative_error,
            diagnostics.current_relative_spread,
            diagnostics.poisson_scaled_residual,
            diagnostics.continuity_scaled_residual,
            diagnostics.minimum_electric_field_v_per_m,
            diagnostics.minimum_electron_m3,
            diagnostics.minimum_hole_m3,
            diagnostics.minimum_conductivity_s_per_m,
        ]
    )
    if not np.all(np.isfinite(values)):
        raise SolverError("solution contains non-finite values or diagnostics")
    if diagnostics.minimum_electric_field_v_per_m <= 0:
        raise SolverError("solution has a non-positive electric field")
    if diagnostics.minimum_electron_m3 < 0 or diagnostics.minimum_hole_m3 < 0:
        raise SolverError("solution has a negative carrier concentration")
    if diagnostics.minimum_conductivity_s_per_m <= 0:
        raise SolverError("solution has non-positive conductivity")
    diagnostic_limit = max(1e-4, 10.0 * tolerance)
    errors = {
        "boundary": diagnostics.maximum_boundary_residual,
        "voltage": diagnostics.voltage_relative_error,
        "current": diagnostics.current_relative_spread,
        "Poisson": diagnostics.poisson_scaled_residual,
        "continuity": diagnostics.continuity_scaled_residual,
    }
    failed = {name: value for name, value in errors.items() if value > diagnostic_limit}
    if failed:
        details = ", ".join(f"{name}={value:.3g}" for name, value in failed.items())
        raise SolverError(f"solution failed independent diagnostics: {details}")


def _solve_case_with_warm_start(
    inputs: ResolvedInputsSI,
    case: BarrierCaseSI,
    *,
    electron_contact_m3: float | None = None,
    hole_contact_m3: float | None = None,
    warm_start: _WarmStart | None = None,
) -> tuple[SimulationResult, _WarmStart]:
    """Solve and validate one barrier case."""
    physics = prepare_case(
        inputs,
        case,
        electron_contact_m3=electron_contact_m3,
        hole_contact_m3=hole_contact_m3,
    )
    started = perf_counter()
    solution, attempts = _solve_scaled(physics, inputs.solver, warm_start)
    diagnostics = _diagnose(solution, physics, perf_counter() - started, attempts)
    _validate_diagnostics(diagnostics, inputs.solver.residual_tolerance)

    arrays = _dimensional_arrays(solution, physics)
    adaptive_scaled = solution.y.copy()
    for array in arrays.values():
        array.setflags(write=False)
    result = SimulationResult(case_name=case.name, diagnostics=diagnostics, **arrays)
    warm = _WarmStart(solution.x.copy(), adaptive_scaled, float(solution.p[0]), physics)
    return result, warm


def solve_case(
    inputs: ResolvedInputsSI,
    case: BarrierCaseSI,
    *,
    electron_contact_m3: float | None = None,
    hole_contact_m3: float | None = None,
) -> SimulationResult:
    """Solve and validate one barrier case without carrying continuation state."""
    result, _warm_start = _solve_case_with_warm_start(
        inputs,
        case,
        electron_contact_m3=electron_contact_m3,
        hole_contact_m3=hole_contact_m3,
    )
    return result


def solve_all_cases(inputs: ResolvedInputsSI) -> tuple[SimulationResult, ...]:
    """Solve cases from highest to lowest barrier, carrying each solution forward."""
    indexed_cases = sorted(
        enumerate(inputs.cases),
        key=lambda item: item[1].electron_barrier_j + item[1].hole_barrier_j,
        reverse=True,
    )
    results: dict[int, SimulationResult] = {}
    warm_start = None
    for index, case in indexed_cases:
        result, warm_start = _solve_case_with_warm_start(inputs, case, warm_start=warm_start)
        results[index] = result
    return tuple(results[index] for index in range(len(inputs.cases)))
