# Charge Injection Conductivity Simulation

This project models how electron and hole injection at metal-dielectric interfaces changes the spatial conductivity
of a dielectric material. The first target is the bipolar charge-injection result shown in **Figure 4(b)** of Wang et al.,
where conductivity is plotted through a 500 um Fe-doped SrTiO3 crystal for several equal electron and hole Schottky
barrier heights.

The repository contains the Stage 1 validated input contract, Stage 2 numerical solver, and Stage 3 local web and
reproducible-output workflow. All three benchmark cases converge with physical-domain, equation-residual, voltage,
current, and refinement checks. The [Stage 4 assessment](docs/reproduction.md) finds qualitative agreement with
Figure 4(b)'s ordering, shape, and conductivity scale under explicit approximations; quantitative agreement remains
unestablished. A real-browser smoke check remains to close the web MVP gate.

**First milestone:** an explicitly approximate, scientifically checked reproduction of Figure 4(b), driven by a
TOML configuration file validated with Pydantic v2 and explored through a local web interface on WSL/Linux.
The MVP uses one NiceGUI + Plotly page: edit parameters, click Run, inspect the three curves, and save the result.
A thin script supports repeatable runs. The author's code/data are unavailable; exact reproduction is not a
prerequisite for this milestone.

See the [critical review](docs/review.md) for the evidence and unresolved assumptions, and the
[staged implementation plan](docs/implementation-plan.md) for architecture, delivery gates, and validation.
The [prototype model contract](docs/model.md) records the equations, coordinate convention, scaling plan, assumptions,
and singular cases that constrain the solver implementation.

## Interactive Workflow

Start the loopback-only web page from the repository root:

```bash
uv run python scripts/web.py
```

Open `http://127.0.0.1:8080` in a browser. Under WSL 2, the Windows browser normally reaches this loopback address
directly. The page loads `configs/figure4b.toml`, exposes temperature, voltage, thickness, recombination, and the three
paired barrier heights, with the remaining material and solver values editable in the expandable parameter section.
Run performs full Pydantic validation and solves in a background process with polled progress and an estimated remaining
time. On WSL/Linux, a preloaded forkserver reduces repeated worker startup cost. A failed solve does not replace the last
accepted plot. Save writes the validated snapshot and exact adaptive-mesh result arrays under `outputs/`.

The first version prioritizes a working scientific loop. Auto-preview, caching, cancellation, pinned comparisons,
browser file editing, and publication polish are deferred until the MVP demonstrates what is useful.

## Reproducible Output

Run and export the benchmark without the web page:

```bash
uv run python scripts/reproduce_figure4b.py --config configs/figure4b.toml
```

The command prints the new output directory. Pass `--output outputs/my-run` to select it; an existing directory is
never overwritten. Relative config and output paths resolve from the current working directory. Every run contains:

- `inputs.json`: validated literature-facing inputs suitable for replay;
- `resolved_inputs_si.json`: the exact immutable SI solver snapshot;
- one CSV per case containing exact adaptive-mesh values with units in column names;
- `metadata.json`: diagnostics, assumptions, solver settings, package/model versions, and available Git provenance;
- `conductivity.png`: a basic static plot in Figure 4(b)'s cm and S/cm display units.

Replay any saved input snapshot with:

```bash
uv run python scripts/reproduce_figure4b.py --config outputs/<run>/inputs.json
```

## Reference

J.-J. Wang, T. J. M. Bayer, R. Wang, J. J. Carter, C. A. Randall, and L.-Q. Chen, "Unexpected significant increase in
bulk conductivity of a dielectric arising from charge injection," *Applied Physics Letters* **110**, 262902 (2017),
[doi:10.1063/1.4990677](https://doi.org/10.1063/1.4990677). A copy is included as [`ref_paper.pdf`](ref_paper.pdf).

## Intended Model

The paper combines four ingredients in a one-dimensional, electronic steady-state approximation with a prescribed
uniform vacancy background. This is not a simulation of long-time ionic redistribution or resistance degradation:

1. Drift conductivity from electrons, holes, and mobile oxygen vacancies.
2. Poisson's equation relating the electric-field gradient to injected space charge.
3. Schottky boundary concentrations at the cathode and anode.
4. A modified continuity equation coupling electron and hole currents through recombination, including the
   ionic-current correction in Equation (5).

The local total conductivity is expected to be evaluated as

```text
sigma(x) = e * [mu_n * n(x) + mu_p * p(x) + 2 * mu_VO * c_VO]
```

with charge signs handled in the transport equations. The first reproduction should use the paper's Fe-doped SrTiO3
case and compare profiles for equal electron and hole barriers of 1.00 eV, 0.85 eV, and 0.80 eV.
Figure 4(b) uses a linear position axis from 0 to 0.05 cm (anode to cathode) and a logarithmic conductivity axis
in S/cm. Equal barriers do not imply symmetric profiles because the effective masses and mobilities differ.

## Baseline Parameters

The initial benchmark described in the paper uses:

| Quantity | Value |
| --- | ---: |
| Temperature | 210 degC (483.15 K) |
| Applied voltage | 40 V |
| Sample thickness | 500 um |
| Relative permittivity | 220 |
| Reported total Fe concentration (not automatically ionized acceptors) | 5.58 x 10^18 cm^-3 |
| Oxygen-vacancy concentration | 2.43 x 10^18 cm^-3 |
| Electron effective mass | 6 m0 |
| Hole effective mass | 12 m0 |
| Electron mobility | 0.56 cm^2 V^-1 s^-1 |
| Hole mobility | 0.41 cm^2 V^-1 s^-1 |
| Oxygen-vacancy mobility | 2.20 x 10^-8 cm^2 V^-1 s^-1 |
| Recombination-rate constant | 1.0 x 10^-8 cm^3 s^-1 |

All values are converted to SI units before calculation. Literature values and numerical settings remain
distinguishable in code and generated metadata; the benchmark uses no fitted values.

Substituting total Fe for charged acceptors in Equation (2) implies a large hole background incompatible with the
stated vacancy-dominated conductivity. The simulator therefore offers two explicit background models:

- `simplified` (default) neglects equilibrium electrons and holes and uses fixed compensation equal to twice the
  specified vacancy density;
- `quenched_equilibrium` treats the specified vacancy density as frozen and solves Fe3+/Fe4+ ionization, intrinsic
  electron-hole equilibrium, Fe conservation, and charge neutrality at the simulation temperature using the Denk
  constants reported by Wang et al. (2016) in [`defectivity_paper.pdf`](defectivity_paper.pdf).

For the benchmark, the calculated mode gives approximately 87.1% charged Fe3+, `p0 = 1.48 x 10^10 cm^-3`, and
`n0 = 21 cm^-3`. Its equilibrium electronic conductivity is about 5.7% of the vacancy conductivity. The web app shows
these resolved values before a run. This mode calculates the equilibrium background only; it does not evolve vacancies
or continuously re-equilibrate Fe in response to injected charge.

As an independently calculated check, the tabulated vacancy concentration and mobility give
`sigma_VO = 1.713 x 10^-6 S/m = 1.713 x 10^-8 S/cm`. They imply a homogeneous dielectric relaxation frequency of
about 140 Hz, rather than the reported experimental 200 Hz. This discrepancy must be reported rather than
removed by silently adjusting mobility or concentration.

## Reproduction Criteria

The first implementation can be considered aligned with Figure 4(b) when it:

- produces finite, nonnegative carrier concentrations and strictly positive total conductivity across the full sample;
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
- whether the Denk quenched-equilibrium constants and the paper's total-density contact convention match the authors'
  unpublished implementation; and
- how sensitive the result is to explicit background and contact assumptions without access to the authors' code.

These uncertainties should be documented rather than absorbed into fitted constants.
The user has already sought the authors' code/data without success. The plan therefore derives a spatial boundary-value
problem from Equations (1), (3), and (5), rather than relying on the missing implementation or treating Equation (6)
as a complete initial-value problem.

## Project Layout

```text
src/charge_injection_sim/  validated inputs, solver, local app, plotting, and serialization
configs/                   validated TOML benchmark inputs
docs/                      review and staged implementation plan
scripts/                   thin reproducible simulation/plot entry points
tests/                     configuration, unit, limiting-case, solver, and refinement tests
tests/data/                digitized comparison curves with provenance (post-MVP)
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

The pre-commit hooks keep `uv.lock` synchronized, run Ruff, and validate common text/configuration errors.
The pytest configuration disables network access. Source packaging, Pydantic v2 validation, the benchmark TOML,
immutable SI conversion, numerical solver, exports, local UI, and scientific tests are implemented. The numerical API
remains available directly:

```python
from charge_injection_sim import load_config, solve_all_cases

inputs = load_config("configs/figure4b.toml").to_si()
results = solve_all_cases(inputs)
```

Each result contains immutable adaptive-mesh arrays in SI units and diagnostics including current density, solver
status, boundary and voltage errors, current spread, independently evaluated equation residuals, extrema, node count,
and elapsed time. Prefer `save_run` for a complete provenance-bearing output bundle.
Code will use type annotations; dedicated type-checker tooling can follow the MVP.
CI is out of scope at this point; checks will be run locally with uv and prek.
