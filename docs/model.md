# Prototype model contract

This note defines the first numerical model implemented in `physics.py` and `solver.py`. It is an approximate interpretation of
Wang et al. (2017), not a complete Fe-doped SrTiO3 defect-chemistry model. The benchmark inputs are recorded in
[`configs/figure4b.toml`](../configs/figure4b.toml), validated in their literature-facing units, and converted once to
immutable SI values before numerical use.

See [Differences from the paper](paper-discrepancies.md) for a short, non-specialist summary of where this prototype
differs from or cannot be verified against the published model.

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

## Numerical scaling

The solver nondimensionalizes before calling `scipy.integrate.solve_bvp`. It uses `x_s=L`, `E_s=V/L`, and `U_s=V`.
For each carrier, its concentration scale is the larger of its contact density and the density whose electronic
conductivity equals `s_v`. Thus `n_s=max(n_c,s_v/(q mu_n))` and `p_s=max(p_a,s_v/(q mu_p))`. This keeps the scales
finite for zero-contact limiting cases and avoids combining carrier variables near zero with coefficients of order
`10^7`, as the provisional `c_v` scale does for the benchmark.

Let `s_s=s_v+q mu_n n_s+q mu_p p_s`, `j_s=s_s E_s`, and define the conductivity fractions
`gamma_v=s_v/s_s`, `gamma_n=q mu_n n_s/s_s`, and `gamma_p=q mu_p p_s/s_s`. With `xi=x/L`, `e=E/E_s`,
`N=n/n_s`, `P=p/p_s`, `u=U/V`, and `J=j/j_s`, the implemented equations are:

```text
P    = (J/e - gamma_v - gamma_n N) / gamma_p
e'   = lambda_p P - lambda_n N
N'   = rho_p N P/e - (N/e)(1 + gamma_v e/J)e'
u'   = e
```

Here `lambda_n=q n_s L/(epsilon E_s)`, `lambda_p=q p_s L/(epsilon E_s)`, and
`rho_p=K p_s L/(mu_n E_s)`. The solver parameter is `log(J)`, guaranteeing a positive trial current without
clipping. Nonzero contact boundary residuals are normalized to their boundary values; a zero-contact limiting case uses
its finite carrier scale as the absolute normalization. The remaining residuals are `u(0)` and `u(1)-1`.

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

Before a profile is accepted, the solver verifies finite fields and nonnegative carrier densities on both the adaptive
mesh and a denser evaluation grid, positive total conductivity, the supported field direction, contact residuals,
independently integrated voltage, spatial current consistency, Poisson and modified-continuity residuals, and stability
under a tighter tolerance or changed mesh. The vacancy-only conductivity and at least one independent transport limit
must pass before comparison with Figure 4(b). Diagnostics are evaluated on at least 1001 points using the dimensional
spline and its derivative; accepted result arrays preserve the adaptive mesh.

The nontrivial independent test uses the compatible `R=0` bipolar invariant
`j_n exp(j_v/j)=constant`, total current, and Poisson's equation to integrate `dx/dE` and `E dx/dE` by quadrature.
A spatially varying `p=0` version is not a valid reference for this model: total-current closure would require
`j=j_n+j_v`, which is incompatible with the invariant unless `j_v`, and hence `E`, is constant.
