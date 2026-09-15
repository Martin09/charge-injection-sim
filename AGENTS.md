# AGENTS.md

## Project

Python simulation of charge-injection-driven conductivity in dielectric materials, initially targeting the bipolar
Fe-doped SrTiO3 model and conductivity profiles in Figure 4(b) of Wang et al. (2017).

- Python >= 3.14, managed exclusively with `uv`.
- Simulation code: `src/charge_injection_sim/`.
- Tests: `tests/`.
- Reproducible plotting entry points: `scripts/`.
- Generated figures and data: `outputs/` (not source-controlled).
- Primary reference: `ref_paper.pdf`.

## Initial Setup

Run `uv sync`, then `uv run prek install` after initializing the Git repository.

## Checks

```bash
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run prek run --all-files
```

## Scientific Rules

- Use SI units internally. Convert literature values at explicit input/output boundaries and include units in parameter
  names, documentation, and plot labels.
- Keep model assumptions and approximations explicit, especially carrier injection, recombination, defect compensation,
  contact boundary conditions, and steady-state assumptions.
- Keep physical parameters separate from solver settings. Do not hide fitted values or numerical tolerances in model
  equations.
- Use `scipy.constants` for physical constants rather than duplicating approximate values.
- Check solver success, residuals, convergence, and physical invariants. Never silently clip negative concentrations or
  accept non-finite results.
- Verify that current density is spatially constant and that the integrated electric field matches the applied voltage,
  within documented numerical tolerances.
- Keep plotting separate from numerical computation. Plotting code consumes simulation results and does not alter them.
- Make simulations deterministic. Seed any stochastic method and record all parameters needed to reproduce an output.
- Preserve raw simulation values when post-processing; do not tune or smooth data solely to resemble the paper.

## Engineering Rules

- Prefer small, typed, testable functions and immutable parameter objects over module-level mutable state.
- Add tests for model behavior, unit conversions, boundary conditions, failure modes, and numerical regressions.
- Compare floating-point results using physically justified tolerances, not exact equality or image pixel matching.
- Validate the implementation first against an analytic or independently calculated limiting case, then against values
  reported in the paper.
- Do not call external services or download data in tests. Treat `ref_paper.pdf` as reference material, not a runtime
  dependency.
- Profile before optimizing; retain a clear reference implementation when introducing a faster numerical path.
- Keep `README.md` aligned with implemented capabilities and clearly distinguish reproduced, approximate, and planned
  results.
- Commit `uv.lock`; do not commit virtual environments, caches, or generated outputs.
- Use Conventional Commits and never skip hooks unless explicitly requested.
