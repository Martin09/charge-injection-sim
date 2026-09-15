"""Minimal local NiceGUI application for interactive benchmark runs."""

from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from queue import Empty as _EmptyQueue
from time import perf_counter
from typing import Any

import plotly.graph_objects as go
from nicegui import ui
from pydantic import ValidationError

from charge_injection_sim.config import (
    MODEL_ASSUMPTIONS,
    InputConfig,
    ResolvedInputsSI,
    load_config,
)
from charge_injection_sim.output import save_run
from charge_injection_sim.runner import start_solver_worker
from charge_injection_sim.solver import SimulationResult, SolverProgressEvent

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "figure4b.toml"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "outputs"

_PROGRESS_POLL_SECONDS = 0.1


@dataclass(frozen=True)
class _CompletedRun:
    inputs: InputConfig
    resolved: ResolvedInputsSI
    results: tuple[SimulationResult, ...]


def _figure(results: tuple[SimulationResult, ...]) -> go.Figure:
    colors = ("black", "red", "blue")
    figure = go.Figure()
    for index, result in enumerate(results):
        figure.add_scatter(
            x=result.position_m * 100.0,
            y=result.total_conductivity_s_per_m * 0.01,
            mode="lines",
            name=result.case_name,
            line={"color": colors[index] if index < len(colors) else None, "width": 2.5},
        )
    figure.update_layout(
        margin={"l": 70, "r": 25, "t": 30, "b": 60},
        xaxis_title="Position, anode to cathode (cm)",
        yaxis_title="Total conductivity (S/cm)",
        yaxis_type="log",
        hovermode="x unified",
        template="plotly_white",
        legend={"orientation": "h", "y": 1.04},
    )
    return figure


def _validation_message(error: ValidationError) -> str:
    lines = []
    for item in error.errors(include_url=False):
        location = ".".join(str(part) for part in item["loc"])
        lines.append(f"{location}: {item['msg']}")
    return "\n".join(lines)


def create_page(config_path: str | Path = DEFAULT_CONFIG) -> None:
    """Create the single-page UI from a validated startup configuration."""
    initial = load_config(config_path)
    completed: _CompletedRun | None = None
    in_flight = False

    ui.add_css("""
        body { background: #f4f1e8; color: #17211b; }
        .q-card { border: 1px solid #d7d0c2; box-shadow: none; }
        .scientific-note { border-left: 4px solid #a33b20; padding-left: 1rem; }
    """)
    with ui.column().classes("w-full max-w-7xl mx-auto p-4 md:p-8 gap-5"):
        ui.label("Charge Injection Conductivity").classes(
            "text-3xl md:text-4xl font-bold tracking-tight"
        )
        ui.label(
            "Approximate 1D Fe-doped SrTiO3 prototype. SI units are used internally; "
            "position runs from anode to cathode."
        ).classes("text-base text-gray-700")
        ui.label(
            "Background approximation: equilibrium electrons and holes are neglected "
            "(n0 = p0 = 0), with fixed compensating charge for the prescribed vacancy density."
        ).classes("scientific-note text-sm text-gray-800")

        status_label = ui.label("No completed run").classes("font-medium text-gray-700")

        with ui.row().classes("w-full items-start gap-5"):
            with ui.card().classes("w-full lg:w-80 p-5 gap-3"):
                ui.label("Run parameters").classes("text-xl font-bold")
                temperature = ui.number(
                    "Temperature (K)", value=initial.experiment.temperature_k, format="%.4g"
                ).classes("w-full")
                voltage = ui.number(
                    "Applied voltage (V)", value=initial.experiment.voltage_v, format="%.4g"
                ).classes("w-full")
                thickness = ui.number(
                    "Thickness (um)", value=initial.experiment.thickness_um, format="%.4g"
                ).classes("w-full")
                recombination = ui.number(
                    "Recombination (cm^3/s)",
                    value=initial.material.recombination_cm3_per_s,
                    format="%.4g",
                ).classes("w-full")
                barrier_inputs = [
                    ui.number(
                        f"{case.name} paired barrier (eV)",
                        value=case.electron_barrier_ev,
                        format="%.4g",
                    ).classes("w-full")
                    for case in initial.cases
                ]

                with ui.row().classes("items-center gap-3"):
                    run_button = ui.button("Run", icon="play_arrow")
                    save_button = ui.button("Save", icon="save").props("outline")
                    save_button.disable()
                    busy = ui.spinner(size="lg")
                    busy.set_visibility(False)
                error_label = ui.label().classes("whitespace-pre-wrap text-red-800")
                progress_row = ui.column().classes("w-full gap-1")
                with progress_row:
                    progress_bar = ui.linear_progress(value=0.0, show_value=False)
                    progress_detail = ui.label("").classes("text-xs text-gray-600")
                progress_row.set_visibility(False)

            with ui.column().classes("grow min-w-0 gap-5"):
                chart = ui.plotly(_figure(())).classes("w-full h-[32rem]")
                diagnostics_card = ui.card().classes("w-full p-5")
                with diagnostics_card:
                    ui.label("Diagnostics will appear after a successful run.")

        with ui.expansion("Fixed material and numerical parameters", icon="science").classes(
            "w-full bg-white"
        ):

            def _label(key: str) -> str:
                return key.replace("_", " ")

            def _number(label: str, value: float, integer: bool = False) -> Any:
                element = ui.number(label, value=value, format="%.6g").classes("w-full")
                if integer:
                    element.props("step=1")
                return element

            material_inputs = [
                _number(_label(key), value)
                for key, value in initial.material.model_dump().items()
                if key != "recombination_cm3_per_s"
            ]
            material_keys = [
                key for key in initial.material.model_dump() if key != "recombination_cm3_per_s"
            ]
            solver_inputs = []
            for key, value in initial.solver.model_dump().items():
                solver_inputs.append(
                    (
                        key,
                        _number(_label(key), value, integer=isinstance(value, int)),
                        isinstance(value, int),
                    )
                )

        with ui.expansion("Model assumptions", icon="info").classes(
            "w-full bg-white scientific-note"
        ):
            for assumption in MODEL_ASSUMPTIONS:
                ui.label(f"- {assumption}")

    def edited_snapshot() -> InputConfig:
        data: dict[str, Any] = initial.model_dump()
        data["experiment"].update(
            temperature_k=temperature.value,
            voltage_v=voltage.value,
            thickness_um=thickness.value,
        )
        data["material"]["recombination_cm3_per_s"] = recombination.value
        for key, material_input in zip(material_keys, material_inputs, strict=True):
            data["material"][key] = material_input.value
        for key, solver_input, _is_integer in solver_inputs:
            data["solver"][key] = solver_input.value
        for case, barrier_input in zip(data["cases"], barrier_inputs, strict=True):
            case["electron_barrier_ev"] = barrier_input.value
            case["hole_barrier_ev"] = barrier_input.value
        return InputConfig.model_validate(data)

    ui_inputs = [
        temperature,
        voltage,
        thickness,
        recombination,
        *barrier_inputs,
        *material_inputs,
        *(element for _key, element, _is_int in solver_inputs),
    ]

    def mark_edited() -> None:
        if completed is not None:
            status_label.set_text("Inputs changed; plot shows the last completed run")

    for field in ui_inputs:
        field.on_value_change(lambda _event: mark_edited())

    worker_process: Any = None
    progress_queue: Any = None
    pending_snapshot: InputConfig | None = None
    pending_resolved: ResolvedInputsSI | None = None
    run_started_at = 0.0

    def _progress_fraction(event: SolverProgressEvent) -> float:
        case_steps = max(event.case_steps_total, 1)
        within_case = min(event.case_steps_done, case_steps) / case_steps
        return (event.completed_cases + within_case) / max(event.total_cases, 1)

    def _eta_text(event: SolverProgressEvent) -> str:
        fraction = max(_progress_fraction(event), 0.02)
        elapsed = perf_counter() - run_started_at
        remaining = elapsed * (1.0 - fraction) / fraction
        return f"elapsed {elapsed:.0f} s, ~{remaining:.0f} s remaining"

    def _progress_detail(event: SolverProgressEvent) -> str:
        return (
            f"{event.current_case} (case {event.completed_cases + 1} of "
            f"{event.total_cases}, continuation step {event.case_steps_done} of "
            f"{event.case_steps_total}, {event.node_count} nodes): {_eta_text(event)}"
        )

    def accept_run(snapshot: InputConfig, resolved: ResolvedInputsSI, results) -> None:
        nonlocal completed
        completed = _CompletedRun(snapshot, resolved, results)
        chart.update_figure(_figure(results))
        diagnostics_card.clear()
        with diagnostics_card:
            ui.label("Accepted diagnostics").classes("text-xl font-bold")
            rows = []
            for result in results:
                diagnostics = result.diagnostics
                maximum_error = max(
                    diagnostics.maximum_boundary_residual,
                    diagnostics.voltage_relative_error,
                    diagnostics.poisson_scaled_residual,
                    diagnostics.continuity_scaled_residual,
                )
                rows.append(
                    {
                        "case": result.case_name,
                        "elapsed_s": f"{diagnostics.elapsed_seconds:.3f}",
                        "nodes": diagnostics.node_count,
                        "current": f"{diagnostics.current_density_a_per_m2:.6g}",
                        "max_residual": f"{maximum_error:.3g}",
                    }
                )
            ui.table(
                columns=[
                    {"name": "case", "label": "Case", "field": "case"},
                    {"name": "elapsed_s", "label": "Elapsed (s)", "field": "elapsed_s"},
                    {"name": "nodes", "label": "Nodes", "field": "nodes"},
                    {"name": "current", "label": "Current (A/m^2)", "field": "current"},
                    {
                        "name": "max_residual",
                        "label": "Max scaled error",
                        "field": "max_residual",
                    },
                ],
                rows=rows,
                row_key="case",
            ).classes("w-full")
        elapsed = sum(result.diagnostics.elapsed_seconds for result in results)
        status_label.set_text(f"Last completed run passed ({elapsed:.3f} s solver time)")
        save_button.enable()

    def reject_run(message: str) -> None:
        error_label.set_text(f"Run failed: {message}")
        status_label.set_text("Last run failed; no new figure was accepted")

    def finish_run() -> None:
        nonlocal in_flight, worker_process, progress_queue
        in_flight = False
        if worker_process is not None:
            worker_process.join()
        worker_process = None
        progress_queue = None
        progress_timer.deactivate()
        busy.set_visibility(False)
        progress_row.set_visibility(False)
        run_button.enable()
        for field in ui_inputs:
            field.enable()

    def poll_progress() -> None:
        nonlocal worker_process, progress_queue
        if not in_flight:
            return
        terminal: tuple[str, Any] | None = None
        if progress_queue is not None:
            while True:
                try:
                    kind, payload = progress_queue.get_nowait()
                except _EmptyQueue:
                    break
                if kind == "progress":
                    event: SolverProgressEvent = payload
                    fraction = max(_progress_fraction(event), 0.02)
                    progress_bar.set_value(min(fraction, 1.0))
                    progress_detail.set_text(_progress_detail(event))
                else:
                    terminal = (kind, payload)
                    break
        if terminal is not None:
            progress_bar.set_value(1.0)
            progress_detail.set_text("Finalizing...")
            snapshot, resolved = pending_snapshot, pending_resolved
            finish_run()
            if terminal[0] == "done":
                assert snapshot is not None and resolved is not None
                accept_run(snapshot, resolved, terminal[1])
            else:
                reject_run(terminal[1])
        elif (
            worker_process is not None
            and not worker_process.is_alive()
            and (progress_queue is None or progress_queue.empty())
        ):
            finish_run()
            reject_run("background solver exited without a result")

    progress_timer = ui.timer(_PROGRESS_POLL_SECONDS, poll_progress, active=False)

    def execute_run() -> None:
        nonlocal in_flight, worker_process, progress_queue
        nonlocal pending_snapshot, pending_resolved, run_started_at
        if in_flight:
            return
        try:
            snapshot = edited_snapshot()
        except ValidationError as error:
            error_label.set_text(_validation_message(error))
            return

        in_flight = True
        run_button.disable()
        for field in ui_inputs:
            field.disable()
        busy.set_visibility(True)
        error_label.set_text("")
        progress_row.set_visibility(True)
        progress_bar.set_value(0.02)
        progress_detail.set_text("Starting solver...")
        status_label.set_text("Solving barrier cases...")
        run_started_at = perf_counter()
        pending_snapshot = snapshot
        pending_resolved = snapshot.to_si()
        worker_process, progress_queue = start_solver_worker(pending_resolved)
        progress_timer.activate()

    def save_completed_run() -> None:
        if completed is None:
            return
        timestamp = datetime.now(UTC).strftime("figure4b-%Y%m%dT%H%M%S-%fZ")
        try:
            path = save_run(
                DEFAULT_OUTPUT_ROOT / timestamp,
                completed.inputs,
                completed.resolved,
                completed.results,
            )
        except Exception as error:
            ui.notify(f"Save failed: {error}", type="negative", timeout=0)
        else:
            ui.notify(f"Saved to {path}", type="positive", timeout=8000)

    run_button.on_click(execute_run)
    save_button.on_click(save_completed_run)


def run_app(config_path: str | Path = DEFAULT_CONFIG) -> None:
    """Build and serve the app on loopback only."""
    create_page(config_path)
    ui.run(
        host="127.0.0.1",
        port=8080,
        title="Charge Injection Simulation",
        show=False,
        reload=False,
    )
