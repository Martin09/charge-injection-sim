"""Quenched Fe defect-equilibrium relations for SrTiO3."""

from dataclasses import dataclass

import numpy as np
from scipy.constants import Boltzmann, elementary_charge
from scipy.optimize import brentq


@dataclass(frozen=True)
class EquilibriumBackground:
    """Equilibrium carrier and Fe charge-state concentrations in SI units."""

    electron_m3: float
    hole_m3: float
    charged_fe3_m3: float
    neutral_fe4_m3: float


@dataclass(frozen=True)
class AnnealedEquilibrium(EquilibriumBackground):
    """Full annealing equilibrium, including the oxygen-vacancy inventory."""

    oxygen_vacancy_m3: float


def denk_equilibrium_constants_si(temperature_k: float) -> tuple[float, float]:
    """Return Denk et al. constants for Fe ionization and electron-hole generation.

    The source expressions use eV and concentrations in cm^-3. Returned K_R3 and
    K_R4 have SI units m^-3 and m^-6, respectively.
    """
    boltzmann_ev_per_k = Boltzmann / elementary_charge
    kr3_cm3 = 2.77e21 * np.exp(
        -(1.18 - 3.7e-4 * temperature_k) / (boltzmann_ev_per_k * temperature_k)
    )
    kr4_cm6 = 7.67e42 * np.exp(
        -(3.3 - 6.0e-4 * temperature_k) / (boltzmann_ev_per_k * temperature_k)
    )
    return float(kr3_cm3 * 1e6), float(kr4_cm6 * 1e12)


def denk_oxygen_exchange_constant_si(temperature_k: float) -> float:
    """Return the Denk R2 constant in m^-9 Pa^1/2.

    Wang et al. (2016), Table 1, tabulates ``K_R2 = n^2 [V_O] sqrt(pO2)``
    in cm^-9 bar^1/2. Pressure and concentration are converted here so the
    equilibrium solver operates entirely in SI units.
    """
    if temperature_k <= 0:
        raise ValueError("temperature must be positive")
    boltzmann_ev_per_k = Boltzmann / elementary_charge
    kr2_cm9_bar_half = 1.82e60 * np.exp(
        -(4.97 - 1.2e-3 * temperature_k) / (boltzmann_ev_per_k * temperature_k)
    )
    return float(kr2_cm9_bar_half * 1e18 * np.sqrt(1e5))


def annealed_fe_equilibrium(
    temperature_k: float,
    oxygen_partial_pressure_pa: float,
    total_fe_m3: float,
) -> AnnealedEquilibrium:
    """Solve the Denk annealing equilibrium including oxygen exchange.

    Doubly ionized oxygen vacancies, Fe3+/Fe4+, electrons, and holes satisfy
    reactions R2-R4, Fe conservation, and charge neutrality.
    """
    if temperature_k <= 0 or oxygen_partial_pressure_pa <= 0 or total_fe_m3 <= 0:
        raise ValueError("temperature, oxygen partial pressure, and total Fe must be positive")

    kr2 = denk_oxygen_exchange_constant_si(temperature_k)
    kr3, kr4 = denk_equilibrium_constants_si(temperature_k)
    constants = np.array([kr2, kr3, kr4])
    if not np.all(np.isfinite(constants)) or np.any(constants <= 0):
        raise ValueError("Denk equilibrium constants are not finite and positive")

    log_kr2 = np.log(kr2)
    log_kr3 = np.log(kr3)
    log_kr4 = np.log(kr4)
    log_total_fe = np.log(total_fe_m3)
    log_pressure = np.log(oxygen_partial_pressure_pa)

    def log_charge_balance(log_hole: float) -> float:
        log_electron = log_kr4 - log_hole
        log_vacancy = log_kr2 - 2.0 * log_electron - 0.5 * log_pressure
        log_charged_fe3 = log_total_fe + log_kr3 - np.logaddexp(log_hole, log_kr3)
        positive_charge = np.logaddexp(np.log(2.0) + log_vacancy, log_hole)
        negative_charge = np.logaddexp(log_charged_fe3, log_electron)
        return float(positive_charge - negative_charge)

    log_hole = brentq(log_charge_balance, -700.0, 700.0, xtol=1e-12, rtol=1e-14)
    hole = float(np.exp(log_hole))
    electron = float(kr4 / hole)
    oxygen_vacancy = float(kr2 / (electron**2 * np.sqrt(oxygen_partial_pressure_pa)))
    charged_fe3 = float(total_fe_m3 * kr3 / (hole + kr3))
    neutral_fe4 = float(total_fe_m3 * hole / (hole + kr3))
    values = np.array([electron, hole, charged_fe3, neutral_fe4, oxygen_vacancy])
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("annealed defect equilibrium produced invalid concentrations")
    return AnnealedEquilibrium(
        electron,
        hole,
        charged_fe3,
        neutral_fe4,
        oxygen_vacancy,
    )


def quenched_fe_equilibrium(
    temperature_k: float,
    total_fe_m3: float,
    oxygen_vacancy_m3: float,
) -> EquilibriumBackground:
    """Solve local Fe redox, carrier generation, conservation, and neutrality.

    The specified doubly ionized oxygen-vacancy concentration is treated as
    frozen during quenching. Fe3+ is the negatively charged acceptor state and
    Fe4+ is neutral on a Ti site.
    """
    if temperature_k <= 0 or total_fe_m3 <= 0 or oxygen_vacancy_m3 <= 0:
        raise ValueError("temperature, total Fe, and oxygen-vacancy density must be positive")

    kr3, kr4 = denk_equilibrium_constants_si(temperature_k)
    concentration_scale = max(total_fe_m3, 2.0 * oxygen_vacancy_m3, kr3, np.sqrt(kr4))

    def neutrality_at_log_hole(log_hole: float) -> float:
        hole = np.exp(log_hole)
        electron = kr4 / hole
        charged_fe3 = total_fe_m3 * kr3 / (hole + kr3)
        return (2.0 * oxygen_vacancy_m3 + hole - charged_fe3 - electron) / concentration_scale

    center = np.log(concentration_scale)
    lower = center - 100.0
    upper = center + 10.0
    if neutrality_at_log_hole(lower) >= 0 or neutrality_at_log_hole(upper) <= 0:
        raise ValueError("could not bracket the quenched defect-equilibrium solution")
    log_hole = brentq(neutrality_at_log_hole, lower, upper, xtol=1e-12, rtol=1e-14)
    hole = float(np.exp(log_hole))
    electron = float(kr4 / hole)
    charged_fe3 = float(total_fe_m3 * kr3 / (hole + kr3))
    # Evaluate both fractions directly: subtracting Fe3+ from total Fe loses
    # the very small but physically meaningful Fe4+ population in n-type cases.
    neutral_fe4 = float(total_fe_m3 * hole / (hole + kr3))
    values = np.array([electron, hole, charged_fe3, neutral_fe4])
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("quenched defect equilibrium produced invalid concentrations")
    return EquilibriumBackground(electron, hole, charged_fe3, neutral_fe4)
