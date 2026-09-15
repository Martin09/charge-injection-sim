from collections.abc import Callable

import pytest
from nicegui import ui

from charge_injection_sim.app import create_page


@pytest.mark.parametrize("invalid_nodes", [None, 201.5])
def test_invalid_node_count_is_reported_before_starting_worker(
    monkeypatch: pytest.MonkeyPatch, invalid_nodes: float | None
) -> None:
    callbacks: dict[str, Callable] = {}
    monkeypatch.setattr(
        ui.button, "on_click", lambda button, callback: callbacks.update({button.text: callback})
    )

    def unexpected_worker(_resolved):
        pytest.fail("invalid inputs must not start a solver worker")

    monkeypatch.setattr("charge_injection_sim.app.start_solver_worker", unexpected_worker)
    with ui.column() as page:
        create_page()
        elements = list(page.descendants())
        nodes = next(
            element
            for element in elements
            if isinstance(element, ui.number) and element.props["label"] == "initial nodes"
        )
        nodes.value = invalid_nodes
        callbacks["Run"]()
        assert any(
            isinstance(element, ui.label) and "solver.initial_nodes:" in element.text
            for element in elements
        )
    page.delete()
