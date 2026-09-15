import csv
import json
from pathlib import Path

import numpy as np
import pytest

from charge_injection_sim import load_config, load_inputs
from charge_injection_sim.output import CSV_COLUMNS, save_run
from charge_injection_sim.solver import SimulationResult, SolverDiagnostics

BENCHMARK_PATH = Path(__file__).parents[1] / "configs" / "figure4b.toml"


def _result(case_name: str) -> SimulationResult:
    position = np.array([0.0, 2e-4, 5e-4])
    field = np.array([7e4, 8e4, 9e4])
    electron = np.array([1e16, 2e16, 3e16])
    hole = np.array([3e16, 2e16, 1e16])
    electron_sigma = np.array([1e-7, 2e-7, 3e-7])
    hole_sigma = np.array([3e-7, 2e-7, 1e-7])
    vacancy_sigma = np.full(3, 1.7e-6)
    diagnostics = SolverDiagnostics(
        status=0,
        message="The algorithm converged to the desired accuracy.",
        current_density_a_per_m2=0.2,
        node_count=3,
        continuation_attempts=1,
        elapsed_seconds=0.01,
        maximum_boundary_residual=1e-8,
        voltage_relative_error=1e-8,
        current_relative_spread=1e-8,
        poisson_scaled_residual=1e-7,
        continuity_scaled_residual=1e-7,
        minimum_electric_field_v_per_m=7e4,
        minimum_electron_m3=1e16,
        minimum_hole_m3=1e16,
        minimum_conductivity_s_per_m=2.1e-6,
        maximum_hole_cancellation_ratio=1.0,
    )
    return SimulationResult(
        case_name=case_name,
        position_m=position,
        electric_field_v_per_m=field,
        electron_m3=electron,
        hole_m3=hole,
        electron_conductivity_s_per_m=electron_sigma,
        hole_conductivity_s_per_m=hole_sigma,
        vacancy_conductivity_s_per_m=vacancy_sigma,
        total_conductivity_s_per_m=electron_sigma + hole_sigma + vacancy_sigma,
        diagnostics=diagnostics,
    )


def test_saved_inputs_round_trip_and_csv_preserves_si_values(tmp_path: Path) -> None:
    inputs = load_config(BENCHMARK_PATH)
    resolved = inputs.to_si()
    results = tuple(_result(case.name) for case in inputs.cases)
    output = tmp_path / "run"

    saved = save_run(output, inputs, resolved, results)

    assert saved == output.resolve()
    assert load_inputs(saved / "inputs.json") == inputs
    assert json.loads((saved / "resolved_inputs_si.json").read_text())["experiment"][
        "thickness_m"
    ] == pytest.approx(5e-4)
    metadata = json.loads((saved / "metadata.json").read_text())
    assert metadata["status"] == "passed"
    assert metadata["model_assumptions"] == list(resolved.model_assumptions)
    assert metadata["solver_settings"] == resolved.solver.model_dump(mode="json")
    assert set(metadata["git"]) == {"revision", "dirty"}
    assert (saved / "conductivity.png").stat().st_size > 0

    with (saved / f"{inputs.cases[0].name}.csv").open(newline="") as data_file:
        rows = list(csv.DictReader(data_file))
    assert tuple(rows[0]) == CSV_COLUMNS
    assert [float(row["position_m"]) for row in rows] == [0.0, 2e-4, 5e-4]
    assert float(rows[1]["total_conductivity_s_per_m"]) == pytest.approx(2.1e-6)


def test_save_refuses_existing_directory(tmp_path: Path) -> None:
    inputs = load_config(BENCHMARK_PATH)
    output = tmp_path / "existing"
    output.mkdir()

    with pytest.raises(FileExistsError):
        save_run(
            output,
            inputs,
            inputs.to_si(),
            tuple(_result(case.name) for case in inputs.cases),
        )


def test_save_rejects_results_for_different_cases(tmp_path: Path) -> None:
    inputs = load_config(BENCHMARK_PATH)

    with pytest.raises(ValueError, match="result cases"):
        save_run(tmp_path / "run", inputs, inputs.to_si(), (_result("different"),))


def test_save_rejects_mismatched_resolved_snapshot(tmp_path: Path) -> None:
    inputs = load_config(BENCHMARK_PATH)
    changed = inputs.model_copy(
        update={
            "experiment": inputs.experiment.model_copy(update={"voltage_v": 20.0}),
        }
    )

    with pytest.raises(ValueError, match="resolved inputs"):
        save_run(
            tmp_path / "run",
            inputs,
            changed.to_si(),
            tuple(_result(case.name) for case in inputs.cases),
        )
