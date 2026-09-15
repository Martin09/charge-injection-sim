"""Serialization and static plotting for reproducible simulation runs."""

import csv
import json
import subprocess
from dataclasses import asdict
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from charge_injection_sim.config import InputConfig, ResolvedInputsSI
from charge_injection_sim.solver import SimulationResult

MODEL_VERSION = "figure4b-drift-bvp-v1"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CSV_COLUMNS = (
    "position_m",
    "electric_field_v_per_m",
    "electron_m3",
    "hole_m3",
    "electron_conductivity_s_per_m",
    "hole_conductivity_s_per_m",
    "vacancy_conductivity_s_per_m",
    "total_conductivity_s_per_m",
)


def _package_version() -> str:
    try:
        return version("charge-injection-sim")
    except PackageNotFoundError:
        return "unknown"


def _git_provenance() -> dict[str, str | bool | None]:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                check=True,
                capture_output=True,
                text=True,
                cwd=PROJECT_ROOT,
            ).stdout
        )
    except FileNotFoundError, subprocess.CalledProcessError:
        return {"revision": None, "dirty": None}
    return {"revision": revision, "dirty": dirty}


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_case_csv(path: Path, result: SimulationResult) -> None:
    arrays = [getattr(result, name) for name in CSV_COLUMNS]
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(CSV_COLUMNS)
        writer.writerows(zip(*arrays, strict=True))


def _write_figure(path: Path, results: tuple[SimulationResult, ...]) -> None:
    colors = ("black", "red", "blue")
    figure = Figure(figsize=(7.2, 4.8), constrained_layout=True)
    FigureCanvasAgg(figure)
    axes = figure.subplots()
    for index, result in enumerate(results):
        axes.plot(
            result.position_m * 100.0,
            result.total_conductivity_s_per_m * 0.01,
            color=colors[index] if index < len(colors) else None,
            label=result.case_name,
        )
    axes.set_yscale("log")
    axes.set_xlabel("Position, anode to cathode (cm)")
    axes.set_ylabel("Total conductivity (S/cm)")
    axes.legend()
    axes.grid(visible=True, which="both", alpha=0.2)
    figure.savefig(path, dpi=180)


def save_run(
    output_directory: str | Path,
    inputs: InputConfig,
    resolved_inputs: ResolvedInputsSI,
    results: tuple[SimulationResult, ...],
) -> Path:
    """Write one complete accepted run, refusing to replace an existing path."""
    if not results:
        raise ValueError("cannot save a run without results")
    if resolved_inputs != inputs.to_si():
        raise ValueError("resolved inputs do not match the validated input snapshot")
    expected_cases = tuple(case.name for case in inputs.cases)
    actual_cases = tuple(result.case_name for result in results)
    if actual_cases != expected_cases:
        raise ValueError("result cases do not match the input snapshot")

    output_path = Path(output_directory).resolve()
    output_path.mkdir(parents=True, exist_ok=False)
    _write_json(output_path / "inputs.json", inputs.model_dump(mode="json"))
    _write_json(output_path / "resolved_inputs_si.json", resolved_inputs.model_dump(mode="json"))

    diagnostics = {result.case_name: asdict(result.diagnostics) for result in results}
    metadata = {
        "created_at_utc": datetime.now(UTC).isoformat(),
        "status": "passed",
        "package_version": _package_version(),
        "model_version": MODEL_VERSION,
        "coordinate_convention": "position_m increases from anode (left) to cathode (right)",
        "raw_data": "CSV files contain exact adaptive-solver-mesh arrays in SI units",
        "model_assumptions": list(resolved_inputs.model_assumptions),
        "solver_settings": resolved_inputs.solver.model_dump(mode="json"),
        "git": _git_provenance(),
        "diagnostics": diagnostics,
    }
    _write_json(output_path / "metadata.json", metadata)
    for result in results:
        _write_case_csv(output_path / f"{result.case_name}.csv", result)
    _write_figure(output_path / "conductivity.png", results)
    return output_path
