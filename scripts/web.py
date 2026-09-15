"""Launch the local Stage 3 web interface."""

from charge_injection_sim.app import run_app

if __name__ in {"__main__", "__mp_main__"}:
    run_app()
