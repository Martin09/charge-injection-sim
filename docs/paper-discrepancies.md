# Differences from the paper

The current simulation is an approximate interpretation of Wang et al. (2017), not an exact copy of the authors'
model or software. The most important differences and uncertainties are summarized below.

## Main physical difference

The paper's general equations include equilibrium electron and hole concentrations, `n0` and `p0`, connected to
charged acceptors such as Fe through charge neutrality. The default simulation sets `n0 = p0 = 0` and assumes a fixed
negative background charge that compensates the oxygen vacancies. Optional calculated modes use the 2016 paper's Denk
constants to calculate Fe3+, Fe4+, `n0`, and `p0`, with the vacancy density either specified or derived from preparation
conditions.

In simplified mode, reported total Fe is saved as provenance but does not affect the graph. In calculated mode, it is
the conserved total `Fe3+ + Fe4+`; it is never substituted directly as charged Fe. For the benchmark the calculation
gives approximately `4.86e18 cm^-3` Fe3+, which nearly compensates `2 c_v`, while the remaining Fe is neutral Fe4+.

## Other differences and uncertainties

- **Vacancy concentration:** `simplified` and `quenched_equilibrium` accept the reported concentration directly.
  `preparation_equilibrium` calculates it from annealing temperature and oxygen partial pressure using the Denk R2
  equilibrium before applying the quench assumption.
- **Numerical method:** The paper reduces the equations to an equation for electric field versus electron concentration
  and mentions a Runge-Kutta solver. This project instead solves a spatial boundary-value problem with SciPy. It is
  intended to represent the same transport equations under the prototype assumptions, but the authors' exact algorithm
  and settings are unavailable.
- **Contact concentrations:** The paper's equation appears to prescribe total carrier density at each contact, while
  nearby text calls it injected density. These are identical when `n0 = p0 = 0`, but would differ in a model with
  nonzero equilibrium carriers.
- **Density of states:** This project uses the standard parabolic-band density-of-states formula, including a factor of
  two for spin. The paper gives effective masses but does not state every convention needed to verify an exact match.
- **Implemented scope:** The project currently calculates the three bipolar conductivity profiles corresponding to
  Figure 4(b). It does not reproduce the paper's impedance spectra, complete unipolar study, automatic barrier sweeps,
  or a time-dependent degradation model.
- **Validation:** The project uses independent limiting-case tests, residual checks, and mesh refinement. The paper
  mentions a trap-free analytic validation, but its original validation data and author code are unavailable.

## Approximations shared with the paper

Some simplifications are not unique to this project. Both approaches use drift-only transport, fixed material
properties, prescribed vacancies, thermally activated contact injection, and a constant bimolecular recombination rate.
The project also retains the paper's modified continuity equation containing the ionic-current correction.

## Practical interpretation

The simulation reproduces the qualitative ordering, shape, and approximate conductivity scale of Figure 4(b), but this
does not prove that its internal carrier and electric-field distributions match the authors' calculation. The optional
background closes the connection among total Fe, charged Fe, oxygen vacancies, and equilibrium carriers under the
stated Denk/quenching assumptions. Quantitative predictions still depend on those constants, the prescribed vacancy
inventory, and the unresolved total-versus-injected contact-density convention.
