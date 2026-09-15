"""Charge-injection conductivity simulation."""

from charge_injection_sim.config import (
    MODEL_ASSUMPTIONS,
    InputConfig,
    ResolvedInputsSI,
    load_config,
    load_inputs,
)
from charge_injection_sim.output import save_run
from charge_injection_sim.solver import (
    SimulationResult,
    SolverDiagnostics,
    SolverError,
    solve_all_cases,
    solve_case,
)

__all__ = [
    "MODEL_ASSUMPTIONS",
    "InputConfig",
    "ResolvedInputsSI",
    "SimulationResult",
    "SolverDiagnostics",
    "SolverError",
    "load_config",
    "load_inputs",
    "save_run",
    "solve_all_cases",
    "solve_case",
]
__version__ = "0.1.0"
