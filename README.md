# Charge Injection Conductivity Simulation

This project will model how electron and hole injection at metal-dielectric interfaces changes the spatial conductivity
of a dielectric material. The first target is the bipolar charge-injection result shown in **Figure 4(b)** of Wang et al.,
where conductivity is plotted through a 500 um Fe-doped SrTiO3 crystal for several equal electron and hole Schottky
barrier heights.

The repository is currently a development scaffold. It does not yet contain a validated numerical implementation or a
reproduced figure.

## Reference

J.-J. Wang, T. J. M. Bayer, R. Wang, J. J. Carter, C. A. Randall, and L.-Q. Chen, "Unexpected significant increase in
bulk conductivity of a dielectric arising from charge injection," *Applied Physics Letters* **110**, 262902 (2017),
[doi:10.1063/1.4990677](https://doi.org/10.1063/1.4990677). A copy is included as [`ref_paper.pdf`](ref_paper.pdf).

## Intended Model

The paper combines four ingredients in a one-dimensional, steady-state model:

1. Drift conductivity from electrons, holes, and mobile oxygen vacancies.
2. Poisson's equation relating the electric-field gradient to injected space charge.
3. Schottky boundary concentrations at the cathode and anode.
4. A continuity equation coupling electron and hole currents through recombination.

The local total conductivity is expected to be evaluated as

```text
sigma(x) = e * [mu_n * n(x) + mu_p * p(x) + 2 * mu_VO * c_VO]
```

with charge signs handled in the transport equations. The first reproduction should use the paper's Fe-doped SrTiO3
case and compare profiles for equal electron and hole barriers of 1.00 eV, 0.85 eV, and 0.80 eV.

## Baseline Parameters

The initial benchmark described in the paper uses:

| Quantity | Value |
| --- | ---: |
| Temperature | 210 degC (483.15 K) |
| Applied voltage | 40 V |
| Sample thickness | 500 um |
| Relative permittivity | 220 |
| Acceptor concentration, primarily Fe | 5.58 x 10^18 cm^-3 |
| Oxygen-vacancy concentration | 2.43 x 10^18 cm^-3 |
| Electron effective mass | 6 m0 |
| Hole effective mass | 12 m0 |
| Electron mobility | 0.56 cm^2 V^-1 s^-1 |
| Hole mobility | 0.41 cm^2 V^-1 s^-1 |
| Oxygen-vacancy mobility | 2.20 x 10^-8 cm^2 V^-1 s^-1 |
| Recombination-rate constant | 1.0 x 10^-8 cm^3 s^-1 |

All values will be converted to SI units before calculation. Literature values, fitted values, and numerical settings
will remain distinguishable in code and generated metadata.

## Reproduction Criteria

The first implementation can be considered aligned with Figure 4(b) when it:

- produces finite, positive carrier concentrations and conductivity across the full sample;
- satisfies the applied-voltage and spatially constant-current constraints within stated tolerances;
- recovers the bulk oxygen-vacancy conductivity when injection is negligible;
- shows enhanced conductivity near both contacts for bipolar injection;
- reproduces the reported ordering and approximate scale of the 1.00, 0.85, and 0.80 eV profiles; and
- writes a figure and machine-readable data from one documented command.

Visual resemblance alone is not validation. Numerical comparisons should use values digitized from the paper or an
independent reference implementation where available.

## Open Modeling Questions

The article gives the governing equations and states that Equation (6) was solved with a Runge-Kutta method, but the
paper alone does not fully specify every closure relation, boundary condition, initial guess, or numerical continuation
strategy needed for an exact reproduction. Before treating results as quantitative, the implementation should resolve:

- how electron and hole concentrations are coupled when integrating the reduced equation;
- how the unknown current density is selected to satisfy both contact and voltage constraints;
- which defect-chemistry values define the equilibrium electron and hole concentrations; and
- whether the authors' supplemental Fortran implementation is available for cross-checking.

These uncertainties should be documented rather than absorbed into fitted constants.

## Planned Layout

```text
src/charge_injection_sim/  model, parameters, solver, and result types
scripts/                   thin reproducible simulation/plot entry points
tests/                     unit, limiting-case, and numerical regression tests
outputs/                   generated data and figures (ignored by Git)
ref_paper.pdf              primary reference paper
```

## Development

Install Python 3.14 and [`uv`](https://docs.astral.sh/uv/), then run:

```bash
uv sync
uv run prek install
```

Run the checks with:

```bash
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run prek run --all-files
```

The pre-commit hooks keep `uv.lock` synchronized, run Ruff, and validate common text/configuration errors. Tests disable
network access so that scientific results do not depend on external services.
