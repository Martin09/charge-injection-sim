"""Solve every case, streaming progress events through the queue."""

import multiprocessing
from contextlib import suppress as contextlib_suppress
from typing import Any

from charge_injection_sim.config import ResolvedInputsSI
from charge_injection_sim.solver import solve_all_cases


def _configure_forkserver() -> None:
    """Preload heavy numeric modules into the fork server once.

    Workers then fork from an already-imported interpreter, so interactive
    runs do not repeatedly pay the module import cost. On platforms without
    fork support the default start method is used and this is a no-op.
    """
    if "forkserver" not in multiprocessing.get_all_start_methods():
        return
    # An already-running fork server keeps its established preload.
    with contextlib_suppress(RuntimeError):
        multiprocessing.set_forkserver_preload(
            ["numpy", "scipy.integrate", "charge_injection_sim.solver"]
        )


def _solver_worker(resolved: ResolvedInputsSI, queue: Any) -> None:
    """Solve every case, streaming progress events through the queue."""
    try:
        results = solve_all_cases(resolved, progress=lambda event: queue.put(("progress", event)))
    except Exception as error:  # sent back so the UI thread can present it
        queue.put(("error", f"{error}"))
    else:
        queue.put(("done", results))


def start_solver_worker(resolved: ResolvedInputsSI) -> tuple[Any, Any]:
    """Start one worker solve; return the process and its progress queue."""
    _configure_forkserver()
    if "forkserver" in multiprocessing.get_all_start_methods():
        context = multiprocessing.get_context("forkserver")
    else:
        context = multiprocessing.get_context()
    queue: Any = context.Queue()
    process: Any = context.Process(target=_solver_worker, args=(resolved, queue), daemon=True)
    process.start()
    return process, queue
