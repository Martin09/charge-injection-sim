# Stage 4 — Figure 4(b) reproduction assessment

Evaluated 2026-09-15 using solver revision `dd6aff3ac955260294253547fdb77b2432f05048`,
model `figure4b-drift-bvp-v1`, and the locked dependencies. This assessment adds a UI validation fix;
the numerical model and benchmark inputs are unchanged.

**Outcome:** the prototype reproduces the qualitative Figure 4(b) ordering, shape, and conductivity
scale under its documented approximations. Numerical checks and saved-input replay pass. The measured
edit/Run cost supports interactive exploration without further optimization. Quantitative agreement
with the paper has not been established. A visual browser smoke check remains to close the web MVP gate.

## Run and inspect

From the repository root:

```bash
uv sync
uv run python scripts/reproduce_figure4b.py --config configs/figure4b.toml --output outputs/stage4-benchmark
uv run python scripts/web.py
```

Open `http://127.0.0.1:8080`. Run the defaults, then change **Applied voltage (V)** from **40** to
**44**, Run, and Save. Keep all other inputs unchanged. Replay the saved edit with:

```bash
uv run python scripts/reproduce_figure4b.py --config outputs/<saved-run>/inputs.json --output outputs/stage4-voltage44-replay
```

Paths passed to the reproduction script are relative to the working directory; the UI uses the
repository's `configs/figure4b.toml` and `outputs/`. Choose fresh output names when repeating these
commands: existing directories are intentionally not overwritten.

The generated benchmark figure is [conductivity.png](../outputs/stage4-benchmark/conductivity.png)
(available after running the command). Each directory also contains exact adaptive-mesh CSVs in SI,
validated inputs, resolved SI parameters, and metadata with diagnostics, assumptions, and Git state.
Generated artifacts remain ignored by Git; the table below records the assessed result in source control.

This evaluation saved runs under `outputs/stage4-benchmark`, `outputs/stage4-worker-benchmark`,
`outputs/stage4-voltage44`, `outputs/stage4-voltage44-replay`, and `outputs/stage4-refined`.
The CLI benchmark metadata records a clean tree at the revision above; later runs record a dirty tree
while the UI fix was being developed. The representative edit was validated through `InputConfig`
before solving, using the same worker entry point as the UI.

## Generated result and paper comparison

Compare against **PDF page 5, printed page 262902-4, Figure 4(b)** of `ref_paper.pdf`.
Both figures use anode left, cathode right, `x = 0–0.05 cm`, and logarithmic conductivity in S/cm.
The generated PNG autoscales its vertical range; the paper shows approximately `10^-8–10^-5 S/cm`.
The following values come from the raw adaptive mesh, with no smoothing or fitting. Minimum positions
are mesh estimates, not separately optimized extrema.

| Paired barrier (eV) | Anode (S/cm) | Minimum (S/cm) | Minimum x (cm) | Cathode (S/cm) | Current (A/m²) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1.00, black | 2.47248e-8 | 2.41439e-8 | 0.05000 | 2.41439e-8 | 0.195351 |
| 0.85, red | 2.35737e-7 | 1.50880e-7 | 0.04000 | 1.53412e-7 | 1.321605 |
| 0.80, blue | 6.76450e-7 | 2.60729e-7 | 0.02775 | 3.53935e-7 | 2.397385 |

- **Ordering passes:** blue > red > black throughout a 1001-point common-grid comparison.
- **Near-flat baseline passes:** the 1.00 eV curve varies by only 2.41% relative to its minimum,
  at a few `10^-8 S/cm`, consistent with the paper's black trace. It is not exactly vacancy-only:
  finite injection raises it above `sigma_v = 1.713047e-8 S/cm`.
- **Interior minima/contact enhancement pass for the lower barriers:** the 0.85 eV anode and cathode
  are respectively 1.56 and 1.017 times its minimum; the 0.80 eV ratios are 2.59 and 1.36.
  The red cathode enhancement is shallow. The black minimum lies at the cathode; requiring a U shape
  for the nearly flat high-barrier case would misstate the result.
- **Approximate scale passes:** red and blue are a few `10^-7 S/cm`, with blue's largest enhancement
  near the anode, as in the paper. The stronger injection raises the interior as well as the contacts.
- **Asymmetry is retained:** minima are displaced toward the cathode, especially for red. Unequal DOS
  and mobilities do not imply mirror symmetry even with paired barriers. The paper also shows asymmetry.

These are visual, order-of-magnitude comparisons to the published panel, not digitized measurements.
The exact differences in contact heights, minimum depths/locations, and curvature remain unquantified;
there is no defensible percent-error claim against Figure 4(b) yet.

## Numerical and reproducibility evidence

All three default cases passed solver status, finite/nonnegative carrier, positive field/conductivity,
component-sum, boundary, independent voltage quadrature, current, Poisson, and modified-continuity checks.
The largest observed default diagnostic across the three cases was:

| Diagnostic | Observed maximum | Benchmark assessment target |
| --- | ---: | ---: |
| Contact/boundary scaled residual | 1.18e-7 | 1e-5 |
| Voltage relative error | 5.17e-7 | 1e-5 |
| Spatial current relative spread | 9.27e-16 | 1e-5 |
| Independently evaluated Poisson scaled residual | 6.74e-5 | 1e-4 |
| Independently evaluated modified-continuity scaled residual | 1.21e-5 | 1e-4 |

These are the assessment targets; the production per-run diagnostic rejection threshold is
`max(1e-4, 10 * residual_tolerance)` for all five errors. The benchmark also meets the stricter
contact/voltage targets checked in the regression suite. See [model.md](model.md) for scaling and
dense-grid checks. Current reconstruction alone does not
validate the differential equations. Final adaptive node counts were 201, 498, and 202.

Refinement changed residual tolerance from `1e-5` to `3e-6` and initial nodes from 201 to 301,
with maximum nodes held at 20000. Comparing linearly interpolated raw conductivity on 1001 common
positions gave maximum relative changes of `2.06e-8`, `9.38e-6`, and `6.89e-4` respectively.
The worst profile change is **0.069%**, below the 0.5% target. Maximum relative changes in current,
raw-mesh minimum, and maximum across cases were `4.26e-8`, `2.92e-6`, and `3.75e-8`.
Interpolation here is a declared comparison operation; it does not replace saved raw values.

At 44 V the three currents were **0.215194, 1.502385, and 2.771742 A/m²**. All per-run checks passed.
Reloading the saved 44 V `inputs.json` and solving afresh reproduced the adaptive positions and total
conductivities within relative tolerance `1e-12` (zero absolute tolerance). Timing and creation/Git
metadata are not expected to match. The scientific suite additionally checks the analytic neutral
Ohmic limit and the independent recombination-free bipolar quadrature reference before paper fidelity
is inferred; these are stronger checks than visual resemblance.

## Timing and workflow assessment

Environment: WSL2 x86-64 Linux `5.15.153.1-microsoft-standard-WSL2`, Python 3.14.0,
NumPy 2.5.3, SciPy 1.18.1, NiceGUI 3.16.0; repository on the Windows-mounted `/mnt/c` filesystem.
These are single observations on this host, not latency guarantees or a statistical benchmark.

| Run | Sum of case solver times (s) | Worker start through received result and process join (s) |
| --- | ---: | ---: |
| Default 40 V, first worker in session | 0.182 | 3.574 |
| Representative 44 V edit, same session | 0.173 | 0.199 |

Times use `perf_counter`. The first worker measurement includes forkserver/module startup; the edit
reuses the preloaded forkserver but solves from fresh numerical initial conditions. Each emitted seven
progress events. Measurements exclude UI/server import, browser rendering, the UI's 0.1 s polling delay,
and CSV/PNG export. The separate CLI benchmark recorded 0.184 s total solver time. Solver time displayed
in the UI should therefore not be mistaken for total click-to-plot latency.

The subsecond warm edit path is useful for deliberate single-parameter exploration. Initial startup
cost is more noticeable than solving. There is no evidence here to justify another optimization layer,
caching, or parallel sweeps; keep explicit Run and the existing progress display.

Headless execution of the actual UI callbacks exercised default Run, 44 V edit/Run, Save, invalid input,
and a two-node/tight-tolerance solver failure. Successful runs produced plots/diagnostics and a saved
bundle; the failed run preserved Save for the last accepted snapshot. An empty node-count field exposed
premature `int()` conversion; fractional counts were also silently truncated. Passing the raw value to
Pydantic now reports both as validation errors before worker startup, covered by regression tests.

**Remaining manual check:** in a real browser, verify initial rendering, progress responsiveness during
Run, edit/Run/Save, validation messages, failed-run plot preservation, and opening the saved PNG.
Headless callback checks do not establish browser rendering or Windows-to-WSL responsiveness.

Local verification: `uv run pytest` passed all 32 tests; `uv run ruff check .`,
`uv run ruff format --check .`, and `uv run prek run --all-files` passed.

## Known limitations and next decision

- The tabulated vacancy baseline implies **139.96 Hz**, versus the stated experimental **200 Hz**:
  200 Hz is about 43% higher (the model value is about 30% lower). No mobility/density calibration was applied.
- The charged Fe fraction and equilibrium defect chemistry remain unspecified. The model sets `n0=p0=0`
  with fixed compensation equivalent to `2 c_v`; total Fe is not imposed as fully ionized acceptors.
- Uniform prescribed vacancies with finite drift current do not describe long-time ionic redistribution.
  This is a drift-only, isothermal, fixed-background electronic approximation, with conventional parabolic DOS.
- Constant bimolecular recombination approximates indirect-gap SrTiO3. The paper's Figure 4 caption has
  inconsistent K units; the implementation uses the body/equation-consistent `cm³/s` input convention.
- Author code/data and calibrated digitized curves are unavailable. No experimental impedance/time evolution,
  fitted defect chemistry, or exact quantitative reproduction is claimed.

**Recommended next scientific step:** digitize Figure 4(b) with axis calibration and uncertainty, then
quantify profile discrepancies before changing physics. Background/defect-chemistry provenance is the
priority if absolute agreement is the goal. Which observed limitation is most worth addressing next:
quantitative curve comparison, the background discrepancy, or an issue encountered in the browser workflow?
