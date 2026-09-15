# MVP implementation plan

Status: Stages 1–3 implemented; Stage 4 executed and assessed in [reproduction.md](reproduction.md).
Numerical/refinement checks and qualitative comparison pass. The remaining web MVP gate is a visual browser smoke
check; headless UI callbacks have been exercised. The scope remains: **prove the concept with a small working MVP**.
See [review.md](review.md) for the paper assessment and known scientific uncertainties.
The stage descriptions below retain the original design targets. Commands and configuration interfaces are implemented.
Subsequent UI changes made the remaining material/solver fields editable and added a background process with a progress
queue, ETA, and forkserver preloading, superseding the original read-only fields/framework-helper scope below.

## MVP goal and stopping point

**Edit parameters → click Run → inspect the three Figure 4(b) curves → save inputs and results.**

Build a local single-user web application on WSL/Linux, using NiceGUI + Plotly and a reusable Python solver.
TOML and Pydantic v2 provide reproducible inputs. Keep the paper's 1.00, 0.85, and 0.80 eV cases as the default.
The target is a numerically checked qualitative reproduction, with explicit approximations, rather than exact
agreement with unavailable author code/data.

The MVP is complete when this loop works, the scientific checks below pass, and the observed agreement and
disagreement with Figure 4(b) are documented. Stop there and evaluate usefulness before adding features.
CI is out of scope. Use local checks and the existing uv/prek workflow.

### Scope boundary

| Build for the MVP | Defer until a demonstrated need |
| --- | --- |
| One local page, numeric fields, Run button, one three-curve plot | Dashboards, polished layouts, slider systems, authentication, remote hosting |
| Pydantic errors shown at Run; unit labels and basic busy/error status | Continuous validation, auto-run, debounced previews, preview/validated modes |
| One in-flight calculation using the framework's background helper | Job manager, queues, custom process lifecycle, cancellation, progress streaming |
| Current result and its parameter snapshot | Caching, pinned comparisons, run history, UI-dependent warm starts |
| Load the benchmark TOML at startup; save the actual run inputs | Browser TOML upload/editor, full import/export workflow, config migrations |
| CSV results, JSON run metadata/config, a basic saved figure | NPZ bundles, artifact manifests/hashes, publication styling, multiple export formats |
| Thin reproducible script calling the same solver | Packaged multi-command CLI, service/API layer, plugin architecture |
| Necessary numerical tests and a local manual UI smoke check | Browser automation infrastructure, broad coverage targets, CI |
| NumPy/SciPy implementation and elapsed-time measurement | Numba/Cython, analytic Jacobians, multi-worker sweeps, performance budgets |

These deferrals reduce application complexity, not scientific checks. Non-finite or negative final densities,
failed solvers, and unmet voltage/current constraints must still produce failures rather than attractive plots.

## Stage 1 — Small package, inputs, and model contract

**Deliverables:** an importable `src/` package, one benchmark TOML file, a Pydantic input model, and a short model note.

1. Enable source packaging with a standard backend such as Hatchling; remove `tool.uv.package = false`.
   Keep Python >= 3.14 and use uv for dependencies. No console-script framework is needed.
2. Add Pydantic v2. Load TOML with standard-library `tomllib`; use Pydantic's JSON serialization for saved
   resolved inputs. A TOML writer is unnecessary for the MVP.
3. Use small frozen models for physical inputs and solver settings, with unknown keys forbidden and finite-value
   checks. Reject invalid temperature, dimensions, mobilities, densities, barriers, and solver limits.
   Allow zero recombination/background carriers for applicable limiting cases; reject unsupported voltage polarity.
4. Convert literature units once into immutable SI parameters. Use unit-suffixed field names, not free-form unit strings.
   Keep physical parameters, numerical tolerances, and plotting choices separate.
5. Record source/assumption provenance in `docs/model.md` and identify the model assumptions in every run's metadata.
   Use the baseline parameter table in the README; do not infer charged acceptors from total Fe.
6. Document the spatial formulation below, its boundary count, scaling, and singular cases before solving.

**Gate:** imports work without path hacks; config validation and unit-conversion tests pass. The model's equations,
units, coordinate direction, and assumptions are explicit. No web UI is needed to complete this stage.

### Configuration contract

Keep one `configs/figure4b.toml`, grouped into `[material]`, `[experiment]`, `[[cases]]`, and `[solver]`.
Use the README's complete baseline values; the following is only a sketch of the experiment/case sections:

```toml
[experiment]
temperature_k = 483.15
voltage_v = 40.0
thickness_um = 500.0

[[cases]]
name = "case_1"
electron_barrier_ev = 1.00
hole_barrier_ev = 1.00

[[cases]]
name = "case_2"
electron_barrier_ev = 0.85
hole_barrier_ev = 0.85

[[cases]]
name = "case_3"
electron_barrier_ev = 0.80
hole_barrier_ev = 0.80
```

Material fields should include `relative_permittivity`, `oxygen_vacancy_cm3`, `electron_effective_mass_m0`,
`hole_effective_mass_m0`, mobilities ending in `_cm2_per_v_s`, and `recombination_cm3_per_s`.
Keep reported total Fe as provenance rather than a misleading active input to the compensated-background model.
The single model initially fixes `n0=p0=0`; do not build a model-selection or defect-chemistry framework yet.

Expose only essential numerical settings: dimensionless residual tolerance, initial nodes, and maximum nodes.
Initial exploration can start at `1e-5`, 201, and 20000 respectively; these are provisional settings to validate,
not paper values or established defaults. A solver residual tolerance is not a guarantee of relative physical error.

Test conversions: `cm^-3 -> m^-3: x1e6`; `cm^2/(V s) -> m^2/(V s): x1e-4`;
`cm^3/s -> m^3/s: x1e-6`; `um -> m: x1e-6`; `eV -> J: xq`; `S/m -> S/cm: /100`.

### Scientific assumptions and proposed equations

Use homogeneous, isothermal, drift-only electronic transport with a prescribed uniform vacancy background and
its finite drift conductivity. Assume negligible equilibrium electrons/holes and compensating fixed charge
equivalent to `2 c_v`; this is an approximate background, not a solved Fe defect-chemistry model.
Retain the paper's modified continuity equation and constant bimolecular recombination coefficient.
Do not add diffusion, field-dependent barriers/mobility, trap kinetics, or vacancy evolution to the first model.

Set `x=0` at the anode and `x=L` at the cathode. For positive benchmark voltage, `E>0`, `j>0`, and
`integral(E dx)=V`. Let `q>0`, `epsilon=epsilon_0 epsilon_r`, `s_v=2 q mu_v c_v`, and `R=K(np-n0p0)`.
Keeping `n0,p0` in the derivation makes the source equations clear; the MVP evaluates them at zero.
Use `y(x)=(E,n,U)` and an unknown constant current parameter `j`:

```text
p   = (j/E - s_v - q mu_n n) / (q mu_p)
E'  = q [(p - p0) - (n - n0)] / epsilon
j_n = q mu_n n E
j_v = s_v E
n'  = {q R - [q mu_n n + (j_n/j) s_v] E'} / (q mu_n E)
U'  = E
```

The electron equation follows by expanding the printed Eq. (5), `j_n' + (j_n/j) j_v' = q R`.
It is not interchangeable with ordinary electron continuity. Four boundary residuals constrain three fields
plus the unknown current:

```text
p(0) = N_V exp[-q phi_p_eV / (k_B T)]
n(L) = N_C exp[-q phi_n_eV / (k_B T)]
U(0) = 0
U(L) = V
```

`U` is accumulated voltage drop, not electrostatic potential. Do not impose both species at both contacts.
Use the conventional parabolic-band DOS `N=2(2 pi m_eff k_B T/h^2)^(3/2)`, recording that convention.
Eq. (4) is interpreted as total contact carrier density; its distinction from injected density vanishes when
`n0=p0=0`. Equal barriers do not imply equal DOS or symmetric profiles.

Nondimensionalize before solving. This spatial formulation avoids using `n` as an integration coordinate at
carrier extrema, but divides by `E` and `j` and reconstructs `p` by subtraction. Check for singular values,
loss of significance, and invalid reconstructed densities. Zero voltage and arbitrary polarity are not supported
by this initial path. Change variables or formulation only if demonstrated numerical difficulty requires it,
with the new equations documented; never mask failures with density floors or silent clipping.

## Stage 2 — Get the physics and three-case solve working

**Deliverables:** pure numerical functions, one SciPy solver, essential diagnostics, and meaningful tests.

1. Implement DOS, contact densities, conductivity components, and the scaled equations using float64 NumPy arrays
   and `scipy.constants`. Precompute material constants once per case; no Pydantic or UI work in callbacks.
2. Independently verify the vacancy baseline (`1.713e-6 S/m`) and the DOS/contact scales recorded in the review.
3. Validate a compatible synthetic charge-neutral Ohmic limit: `E=V/L`, `j=sigma V/L`.
   Also implement one nontrivial independent transport reference before assessing paper agreement. With `R=0`,
   Eq. (5) yields `j_n exp(j_v/j)=constant`; combine it with Poisson and total current to construct an independent
   bipolar quadrature reference. Do not reuse the production residual as the reference. The originally proposed
   spatially varying unipolar specialization was removed because `p=0` also requires `j=j_n+j_v`; together with
   the invariant this permits only constant `j_v` and therefore no nontrivial field profile.
4. Start with dimensionless `scipy.integrate.solve_bvp`, including its unknown-current parameter. Use the
   high-barrier solution to initialize the next case. Add bounded intermediate continuation steps only if needed.
5. Return arrays plus a small diagnostic record: status/message, current, residuals, node count, and elapsed time.
   Limit nodes/continuation attempts; report failure clearly. Avoid implementing a second solver preemptively.
6. Check the numerical contract below, including one tighter-tolerance/different-mesh comparison of the benchmark.
   Document failures and assumptions before displaying profiles as results.

**Gate:** all three cases converge, physical-domain and residual checks pass, and results are stable under
refinement. Analytic/independent limits pass first. If this cannot be achieved, resolve the numerical or physical
issue before investing in the interface. Author code is not a dependency.

### Essential numerical contract

- Finite fields and nonnegative final carrier concentrations on the solver mesh and a denser evaluation grid;
  strictly positive total conductivity, correctly summed components, and the supported field direction.
- Contact-density residuals scaled to their boundary values.
- Voltage error from independent quadrature of the returned field.
- Spatial current consistency plus independently evaluated Poisson and modified-continuity residuals.
  Algebraic reconstruction enforces current by construction, so current consistency alone is insufficient.
- Refinement stability of current, extrema, and conductivity on a common grid.

Starting targets: contact/voltage relative error <= `1e-5`, scaled equation residuals <= `1e-4`, and current/profile
changes <= `0.5%` on refinement. Establish dimensionless scales and absolute tolerances near zero using the
limiting cases; document justified adjustments. These are proposed numerical targets, not paper-matching tolerances.
Run refinement tests when establishing/changing the solver and benchmark, not automatically twice for every UI edit.
Each interactive run still checks solver status, physical domain, and residuals; label it approximate and show failures.

## Stage 3 — Minimal local web page and reproducible output

**Deliverables:** one usable page, basic exports, and one thin reproduction script.

- Add NiceGUI and Plotly through uv, confirming Python >= 3.14 compatibility. Use a single `app.py` initially.
  Bind to loopback and document opening the WSL-hosted page from the Windows browser.
- Load `configs/figure4b.toml` at startup. Start with numeric inputs for temperature, voltage, thickness,
  recombination, and the three paired barrier heights. Other material values remain editable in TOML and are
  shown read-only in the UI; add more controls only when useful. Always show units and the background approximation.
- On Run, validate the whole edited parameter snapshot with the same Pydantic model used for file inputs.
  Show understandable validation errors and preserve the user's edits.
- Permit one in-flight run. Disable Run and show a simple busy indicator while awaiting NiceGUI's supported
  CPU-bound/background helper. Keep the numerical callable serializable and independent of UI state.
  Use framework-managed execution, not a custom worker manager. No cancellation or queue is required initially.
- Plot the three total-conductivity traces with Plotly: linear x in cm, log conductivity in S/cm, anode left,
  cathode right, and black/red/blue cases. Built-in zoom/hover is enough.
- Show elapsed time and pass/fail diagnostics. Preserve the parameter snapshot associated with the plot;
  label it as the last completed run when inputs change. Do not present failed/partial runs as a complete figure.
- Add Save for the current run: a fresh output directory containing the actual validated inputs and resolved SI
  parameters in JSON, a CSV per case, diagnostics/assumptions/version information, and one basic PNG.
  Use existing Matplotlib for PNG export to avoid a separate browser-based image-rendering dependency.
  Save the exact adaptive-mesh arrays; any plot sampling must be separately identified, not replace raw data.
- Keep serialization and plotting separate from calculation. Output units and coordinate conventions must be explicit.
  Record the package/model version, available Git revision and dirty status, and numerical settings needed to rerun.
- Add a thin script using the same loader/solver/export functions. Accept benchmark TOML or a saved validated
  JSON input snapshot so UI edits can be replayed without implementing a TOML writer or a multi-command CLI.

Proposed entry points:

```bash
uv run python scripts/web.py
uv run python scripts/reproduce_figure4b.py --config configs/figure4b.toml
uv run python scripts/reproduce_figure4b.py --config outputs/example-run/inputs.json
```

Document config/output path resolution and refuse overwriting an existing run directory. Use simple exceptions
and nonzero script exit status on failure; no artifact management framework is necessary.

**Gate:** manually exercise startup, edit, Run, invalid input, solver failure, Save, and saved-input replay locally.
The page remains responsive during solving. Automated tests cover the shared validation, numerical failures,
serialization round-trip, and exported units/values; no browser test infrastructure is required for this MVP.
No external services are required at runtime or in scientific tests.

## Stage 4 — Evaluate the MVP and decide the next step

**Deliverable:** a short `docs/reproduction.md` with run instructions, the generated result, observed timing,
qualitative assessment, and known discrepancies.

- Compare barrier ordering, near-flat 1.00 eV baseline, interior minima/contact enhancement for lower barriers,
  and approximate conductivity scale against Figure 4(b). Do not require artificial symmetry.
- Keep the 140 Hz versus experimental 200 Hz baseline discrepancy and unspecified defect chemistry visible.
  Do not fit or smooth values merely to improve visual agreement.
- Record solve time for the three-case benchmark and a representative parameter edit. Assess whether this
  edit/Run loop is useful before committing to latency targets or optimization infrastructure.
- If numerical checks pass but qualitative features fail, call it a solver prototype and investigate assumptions;
  do not declare Figure 4(b) reproduced. Exact quantitative agreement is not the MVP gate.

**MVP exit:** a useful end-to-end web workflow, stable numerical results, reproducible saved runs, and an honest
assessment of qualitative fidelity. Ask which observed limitation is worth addressing next.

## Small initial structure

```text
configs/figure4b.toml
src/charge_injection_sim/
    config.py       Pydantic models, file loading, conversion to SI parameters
    physics.py      constants, contacts, conductivity, scaled equations
    solver.py       SciPy solve, minimal continuation, diagnostics and result type
    output.py       CSV/JSON export and basic static plotting helpers
    app.py          NiceGUI parameter form and Plotly display
scripts/
    web.py
    reproduce_figure4b.py
tests/
    test_config.py
    test_physics.py
    test_solver.py
    test_output.py
```

Use typed functions and small immutable parameter/result objects. Avoid mutable global UI state. Physics and solver
code must not import NiceGUI/Plotly, perform file I/O, or depend on UI configuration objects. Pydantic validates at
the boundary; numerical callbacks use SI values. Result arrays need clear ownership; do not mutate them when plotting.
Split modules further only when size or genuinely separate responsibilities warrant it. No service layer, registry,
database, REST API, or generic solver abstraction is needed to validate the concept.

## Local verification

Run the existing required checks as implementation lands:

```bash
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run prek run --all-files
```

Meaningful tests cover units, invalid inputs, analytic/independent limits, boundaries, failure modes, and numerical
regressions. Use justified floating-point tolerances; no pixel comparisons or external data downloads in tests.
Keep annotations readable; adding a dedicated static checker can follow the MVP rather than block its first run.
Keep README aligned with actual commands. Do not add CI workflows or unneeded development tooling now.

## Later improvements — select based on MVP evidence

1. **Scientific fidelity:** digitize the paper with axis calibration and uncertainty, quantify log-conductivity
   errors, improve background defect chemistry, and assess sensitivity to K/barriers. Expand limiting-case coverage
   (including Mott–Gurney under its ideal-contact assumptions), then other panels/models as needed.
2. **Iteration usability:** first add whichever of more parameter fields, TOML upload/export, or pinned comparisons
   actually helps. Consider caching, auto-preview, cancellation, and richer progress only if measured delays justify them.
3. **Performance:** profile first. Improve scaling/continuation/vectorization before analytic Jacobians or optional
   compiled kernels; preserve a clear reference implementation and numerical equivalence. Add sweep parallelism only
   for independent workloads large enough to benefit. UI background execution is for responsiveness, not speedup.
4. **Engineering/polish:** type-checker integration, browser smoke automation, richer provenance/export formats,
   CLI packaging, and publication styling after the scientific loop is useful. CI remains deferred unless requested.

None of these later improvements is a prerequisite for MVP acceptance.
