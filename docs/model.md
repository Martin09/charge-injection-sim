# Prototype model contract

This note defines the first numerical model before its solver is implemented. It is an approximate interpretation of
Wang et al. (2017), not a complete Fe-doped SrTiO3 defect-chemistry model. The benchmark inputs are recorded in
[`configs/figure4b.toml`](../configs/figure4b.toml), validated in their literature-facing units, and converted once to
immutable SI values before numerical use.

## Provenance and assumptions

The temperature, voltage, thickness, relative permittivity, reported total Fe concentration, oxygen-vacancy
concentration, effective masses, mobilities, recombination coefficient, and contact barriers come from the parameter
paragraph and Figure 4(b) of the paper. The total Fe concentration is retained only as provenance. It is not treated as
the concentration of ionized acceptors because the charged fraction and complete defect chemistry are unspecified.

The prototype assumes:

- homogeneous, one-dimensional, isothermal, drift-only transport;
- a prescribed uniform oxygen-vacancy density with finite drift conductivity but no vacancy evolution;
- negligible equilibrium electrons and holes, `n0 = p0 = 0`;
- fixed compensating charge equivalent to twice the oxygen-vacancy density;
- a constant bimolecular recombination coefficient; and
- no diffusion, field-dependent barriers or mobilities, trap kinetics, or tunneling.

These assumptions are stored with every resolved input snapshot through `model_assumptions`. Changing one defines a
different model variant and must not be hidden as solver tuning.

## Coordinates, signs, and units

All equations are evaluated in SI units. Position `x` increases from the anode at `x = 0` to the cathode at `x = L`.
The initial implementation supports only positive applied voltage: `V > 0`, `E > 0`, and total current density
`j > 0`, with

```text
integral from 0 to L of E(x) dx = V.
```

The elementary charge `q` is positive. Electron and hole concentrations `n` and `p` are in `m^-3`, electric field `E`
is in `V/m`, current density `j` is in `A/m^2`, and conductivity is in `S/m`. The plotted output will convert position
to `cm` and conductivity to `S/cm` explicitly.

The conventional parabolic-band density of states is

```text
N = 2 [2 pi m_eff k_B T / h^2]^(3/2).
```

This factor-of-two convention is part of the model. Contact densities are interpreted as total carrier densities;
with `n0 = p0 = 0`, they are also the injected densities.

## Spatial boundary-value formulation

Let `epsilon = epsilon_0 epsilon_r`, vacancy conductivity `s_v = 2 q mu_v c_v`, and recombination
`R = K (n p - n0 p0)`. The state is `y(x) = (E, n, U)` and the spatially constant current density `j` is an unknown
scalar parameter. Hole concentration is reconstructed from total current:

```text
p   = [j/E - s_v - q mu_n n] / (q mu_p)
E'  = q [(p - p0) - (n - n0)] / epsilon
j_n = q mu_n n E
j_v = s_v E
n'  = {q R - [q mu_n n + (j_n/j) s_v] E'} / (q mu_n E)
U'  = E
```

The electron equation is the expansion of the paper's modified continuity equation,
`j_n' + (j_n/j) j_v' = q R`. Replacing it with ordinary electron continuity would change the model.

Three first-order fields plus the unknown `j` require four scalar boundary residuals:

```text
p(0) = N_V exp[-q phi_p,eV / (k_B T)]
n(L) = N_C exp[-q phi_n,eV / (k_B T)]
U(0) = 0
U(L) = V
```

Here `U` is accumulated voltage drop, not electrostatic potential. The formulation deliberately does not impose both
carrier species at both contacts. Equal barrier energies do not produce equal boundary densities because `N_C` and
`N_V` depend on different effective masses.

## Numerical scaling contract

The Stage 2 solver must nondimensionalize the equations before calling a boundary-value solver. Candidate fixed scales
are `x_s = L`, `E_s = V/L`, concentration scale `c_s = c_v`, conductivity scale `s_s = s_v`, voltage scale `U_s = V`,
and current scale `j_s = s_v V/L`. Contact densities and resulting dimensionless groups must be computed in float64
from the resolved SI values. If these provisional scales produce poorly conditioned boundary residuals, Stage 2 may
choose documented case-dependent scales, but physical inputs and equations must remain unchanged.

Boundary residuals must be normalized by their corresponding contact, field, or voltage scales. A solver's
dimensionless residual tolerance is a numerical control, not a guarantee of relative physical accuracy. Returned
solutions must also be checked independently in dimensional units.

## Singular and invalid cases

The spatial equations divide by `E` and `j`, and reconstruct `p` by subtracting current components. Consequently:

- zero voltage and negative polarity are unsupported by this formulation;
- a zero or sign-changing field, or zero current, is singular;
- cancellation in the reconstruction of `p` can lose precision or yield an invalid negative density;
- vanishing electron mobility makes the electron equation singular; and
- non-finite fields, currents, concentrations, or residuals invalidate a result.

Input validation rejects unsupported polarity and non-positive dimensions, temperature, densities, permittivity,
effective masses, and mobilities. Zero recombination and zero barriers remain valid for controlled limiting cases.
The solver must fail visibly rather than clip densities, insert hidden floors, or return partial curves as results.

## Required solution checks

Before a profile is accepted, Stage 2 must verify finite fields and nonnegative carrier densities on both the adaptive
mesh and a denser evaluation grid, positive total conductivity, the supported field direction, contact residuals,
independently integrated voltage, spatial current consistency, Poisson and modified-continuity residuals, and stability
under a tighter tolerance or changed mesh. The vacancy-only conductivity and at least one independent transport limit
must pass before comparison with Figure 4(b).
