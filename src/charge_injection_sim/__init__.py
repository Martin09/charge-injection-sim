"""Charge-injection conductivity simulation."""

from charge_injection_sim.config import (
    MODEL_ASSUMPTIONS,
    InputConfig,
    ResolvedInputsSI,
    load_config,
)
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
    "solve_all_cases",
    "solve_case",
]
__version__ = "0.1.0"
