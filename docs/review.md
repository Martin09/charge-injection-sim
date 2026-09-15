# Repository and reference review

Reviewed 2026-09-15 against commit `fb29f57` and the local `ref_paper.pdf`.
This is a planning review, not a claim of numerical reproduction.

## Agreed scope

- This pass delivers review and implementation documentation.
- First implementation milestone: qualitative Figure 4(b) prototype with explicit assumptions and numerical checks.
- Author code/data have already been sought and are unavailable. They are not a delivery dependency.
- Latest scope decision: a lightweight MVP with one local NiceGUI + Plotly page, explicit Run, and three curves.
  Load TOML at startup, validate with Pydantic v2, and save JSON inputs/CSV results for replay through a thin script.
  Primary platform WSL/Linux; Python >= 3.14, managed by uv.
- CI is out of scope. Retain local scientific tests, linting, and prek checks. Use type annotations;
  dedicated type-checker tooling and UI automation can follow the MVP.
- Quantitative comparison follows, without disguising model uncertainty as numerical error.

## Repository assessment

The working tree was clean at review start. Git is initialized and contains one scaffold commit.
Present: README, project instructions, reference PDF, locked NumPy/SciPy/Matplotlib dependencies,
pytest configuration, Ruff, coverage configuration, and prek hooks.
Absent at the initial review: source package, tests, configs, web interface, CLI, numerical implementation,
benchmark data, and CI. The lack of CI is accepted for the current scope.

Good foundations include a committed lockfile, explicit Python floor, strict warnings, network isolation for tests,
ignored generated outputs, and unusually clear scientific rules in `AGENTS.md`.

Concrete gaps:

| Finding | Impact | Recommendation |
| --- | --- | --- |
| `tool.uv.package = false`, no build backend | Future `src/` package will not be installed by this setup | Enable packaging with a standard backend and console entry point |
| Pydantic absent | No validated user-input boundary | Add Pydantic v2 with explicit finite-value/unit validation |
| No tests, configured `testpaths = ["tests"]` | `uv run pytest` currently fails with a warnings-as-errors configuration error | Add meaningful foundational tests in Stage 1 |
| No static type checker | Annotations would not be statically verified | Keep annotations; dedicated tooling can follow the MVP |
| Hook coverage excludes tests | Passing hooks is not scientific validation | Run the documented scientific checks locally |
| No machine-readable scientific provenance | Results cannot yet be independently regenerated | Versioned configs, resolved SI metadata, diagnostics, and comparison data |

No performance conclusion is possible before a solver exists. A three-case 1D calculation does not justify a GPU,
distributed execution, or native-extension maintenance in advance.

## Critical README assessment

The original README correctly called the repository a scaffold, identified the target barriers, recorded most
paper parameters accurately, separated SI conversion from literature inputs, and rejected visual-only validation.
Its main weaknesses were scientific specificity and actionable delivery detail:

1. **Total Fe versus charged acceptors:** the parameter table called the reported Fe concentration an acceptor
   concentration without distinguishing ionization/charge state. This can make the proposed baseline inconsistent.
2. **Steady state:** the model freezes the vacancy background while permitting its drift contribution. Calling it
   simply steady-state transport obscures the lack of long-time ionic mass conservation and evolving defects.
3. **Continuity:** Equation (5) is modified by an ionic-current term. Generic bipolar continuity is not equivalent.
4. **Figure contract:** electrode orientation, axis units, logarithmic scale, and expected asymmetry were unstated.
5. **Acceptance:** “approximate scale” had no operational comparison procedure, uncertainty handling, or convergence gate.
6. **Software contract:** the layout did not describe configuration validation, package installation, CLI failures,
   result metadata, solver diagnostics, or dependency direction.
7. **Development status:** the checks were listed without explaining that pytest fails before collection today.
8. **External code:** “supplemental Fortran” is imprecise; reference 42 links an external author-hosted package page.
   The plan must stand independently of that unavailable source.

The revised README is the short entry point; the implementation plan holds the engineering detail.
The subsequent scope refinements prioritize a minimal interactive proof of concept. Use framework-managed background
execution for one run and basic diagnostics. Cancellation, caching, preview modes, run history, and UI automation
are deferred; CI is out of scope. The earlier expanded application design is superseded by the MVP plan.

## What the paper actually specifies

Page references below use printed article pages 262902-1 through 262902-5. The PDF also has a publisher cover page.

| Source | Evidence relevant to implementation |
| --- | --- |
| p. 1, Eq. (1) | Drift-only current, local conductivity times field, spatially constant total current |
| p. 1, Eq. (2) | Charge neutrality of equilibrium carriers and charged defects |
| p. 2, Eq. (3) | Poisson source is injected hole density minus injected electron density |
| p. 2, Eq. (4) | Electron contact density at cathode, hole contact density at anode, thermally activated barriers |
| p. 2, Eq. (5) | Recombination `K (n p - n0 p0)` plus ionic-current correction on the left-hand side |
| p. 2, Eq. (6) | Reduced `dE/dn` equation; Runge–Kutta mentioned without a complete boundary-solving algorithm |
| p. 2, parameter paragraph | T, V, thickness, total Fe, vacancy density, permittivity, masses, mobilities, and assumed K |
| p. 4, Fig. 4(b) | Three equal-barrier cases, 1.00/0.85/0.80 eV; anode left, cathode right |
| p. 4, discussion | Direct band-to-band recombination approximates indirect-gap SrTiO3; total conductivity less sensitive to K |
| p. 5, ref. 42 | External code link, not sufficient documentation of algorithm or closure |

Figure 4(b) has `x = 0 ... 0.05 cm` and logarithmic `sigma` in `S/cm` (roughly `10^-8 ... 10^-5`).
The black 1.00 eV curve is nearly flat at a few `10^-8 S/cm`; the red 0.85 eV and blue 0.80 eV curves
are of order `10^-7 S/cm`, with enhanced contact regions and a broad interior minimum.
These are visual estimates, not digitized measurements or acceptance fixtures. Figure 4(d)'s 0.6 eV cases are different.
The 0/90/130 s symbols in panel (a) are experimental observations; they are not simulation time inputs for panel (b).

### Numerical sanity checks and discrepancies

Using `scipy.constants` and the tabulated values:

- `c_VO = 2.43e24 m^-3`, `mu_VO = 2.20e-12 m^2/(V s)`.
- `sigma_VO = 2 q mu_VO c_VO = 1.713047e-6 S/m = 1.713047e-8 S/cm`.
- `sigma_VO / (2 pi epsilon_0 epsilon_r) = 139.96 Hz` for a uniform medium, versus the stated experimental 200 Hz.
  These differ by about 43%; measurement/model provenance must be investigated before any calibration.
- If `[A-] = [Fe] = 5.58e24 m^-3`, Eq. (2) implies `p0 - n0 = 7.20e23 m^-3`.
  With the stated hole mobility, even the minimum implied hole conductivity is `4.73 S/m`, about
  `2.76e6` times the vacancy conductivity. This does not establish that the paper is wrong: the missing information
  may be the charged Fe fraction or other defect chemistry. It does establish that total Fe cannot safely be
  substituted as fully ionized acceptors in this benchmark.
- With the conventional parabolic-band DOS formula `N = 2 (2 pi m_eff k_B T / h^2)^(3/2)`,
  `N_C = 7.538e26 m^-3`, `N_V = 2.132e27 m^-3` at 483.15 K.
  The standard degeneracy convention is an explicit prototype assumption.

| Equal barriers (eV) | Cathode electron density (m^-3) | Anode hole density (m^-3) |
| --- | ---: | ---: |
| 1.00 | 2.794e16 | 7.901e16 |
| 0.85 | 1.025e18 | 2.900e18 |
| 0.80 | 3.407e18 | 9.636e18 |

These checks test scales and units; they are not solved profiles. The unequal DOS values explain why equal
barriers should not be enforced as mirror-symmetric boundary densities.

### Explicit prototype assumptions

- Homogeneous 1D, isothermal, drift-only transport; no diffusion, field-dependent mobility, tunneling, or
  field-induced Schottky barrier lowering beyond the supplied barrier values.
- Prescribed uniform vacancy density with finite ionic drift conductivity; no vacancy evolution equation.
- The default uses `n0 = p0 = 0` with fixed compensating charge equivalent to `2 c_VO`. The optional quenched mode
  calculates `n0`, `p0`, Fe3+, and Fe4+ from the Denk equilibria without treating total Fe as fully ionized.
- Equation (4)'s Schottky densities are treated as injected increments, following the prose. They are added to `n0` or
  `p0` in calculated-background mode; this differs from reading the typeset equation as a total-density boundary.
- Constant bimolecular `K = 1e-14 m^3/s`, not a fitted parameter or an SRH/Auger kinetics implementation.
  The Figure 4 caption prints mobility-like units for K; the body and dimensional balance of Eq. (5)
  require concentration-inverse time units (`cm^3/s` in the paper). Record this inconsistency.
- Retain the paper's modified continuity term for fidelity. Document its implications instead of asserting that
  this is a general, fully species-conserving mixed-conductor model.

Changing any of these defines a new model variant. Solver tuning must not silently change the physics.

## Review outcome

Proceed with an approximate prototype, not a claim of exact reproduction. The highest-priority work is a
dimensionally checked spatial formulation and independently verified numerical behavior. Once that exists,
compare all three curves on a common spatial grid and report remaining disagreement openly.
