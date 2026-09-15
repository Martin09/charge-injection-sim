"""Validated configuration loading and conversion to SI units."""

import tomllib
from enum import StrEnum
from pathlib import Path
from typing import Annotated, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator
from scipy.constants import electron_mass, elementary_charge

from charge_injection_sim.defect_chemistry import annealed_fe_equilibrium, quenched_fe_equilibrium

COMMON_MODEL_ASSUMPTIONS = (
    "homogeneous one-dimensional isothermal drift-only transport",
    "uniform frozen oxygen-vacancy density with finite drift conductivity",
    "constant bimolecular recombination coefficient",
    "no diffusion, field-dependent transport, trap kinetics, or vacancy evolution",
)
MODEL_ASSUMPTIONS = (
    *COMMON_MODEL_ASSUMPTIONS,
    "negligible equilibrium electron and hole densities (n0 = p0 = 0)",
    "fixed compensating charge equivalent to twice the oxygen-vacancy density",
)

PositiveFloat = Annotated[float, Field(gt=0, allow_inf_nan=False)]
NonNegativeFloat = Annotated[float, Field(ge=0, allow_inf_nan=False)]
PositiveInt = Annotated[int, Field(gt=0)]


class FrozenModel(BaseModel):
    """Base for immutable models with a closed configuration schema."""

    model_config = ConfigDict(frozen=True, extra="forbid")


class BackgroundModel(StrEnum):
    """Available treatments of the neutral equilibrium background."""

    SIMPLIFIED = "simplified"
    QUENCHED_EQUILIBRIUM = "quenched_equilibrium"
    PREPARATION_EQUILIBRIUM = "preparation_equilibrium"


class BackgroundInputs(FrozenModel):
    """Choice of equilibrium-background treatment."""

    model: BackgroundModel = BackgroundModel.SIMPLIFIED


class PreparationInputs(FrozenModel):
    """Sample preparation conditions used to establish the vacancy inventory."""

    annealing_temperature_k: PositiveFloat
    oxygen_partial_pressure_bar: PositiveFloat


class MaterialInputs(FrozenModel):
    """Material parameters expressed in the literature-facing units."""

    relative_permittivity: PositiveFloat
    reported_total_fe_cm3: PositiveFloat
    oxygen_vacancy_cm3: PositiveFloat
    electron_effective_mass_m0: PositiveFloat
    hole_effective_mass_m0: PositiveFloat
    electron_mobility_cm2_per_v_s: PositiveFloat
    hole_mobility_cm2_per_v_s: PositiveFloat
    oxygen_vacancy_mobility_cm2_per_v_s: PositiveFloat
    recombination_cm3_per_s: NonNegativeFloat


class ExperimentInputs(FrozenModel):
    """Experimental conditions expressed in user-facing units."""

    temperature_k: PositiveFloat
    voltage_v: PositiveFloat
    thickness_um: PositiveFloat


class BarrierCaseInputs(FrozenModel):
    """Contact barriers for one conductivity profile."""

    name: str = Field(min_length=1, pattern=r"^[a-z0-9][a-z0-9_-]*$")
    electron_barrier_ev: NonNegativeFloat
    hole_barrier_ev: NonNegativeFloat


class SolverSettings(FrozenModel):
    """Numerical controls, kept separate from physical inputs."""

    residual_tolerance: Annotated[float, Field(gt=0, lt=1, allow_inf_nan=False)]
    initial_nodes: Annotated[int, Field(ge=2)]
    maximum_nodes: PositiveInt

    @model_validator(mode="after")
    def validate_node_limits(self) -> Self:
        if self.maximum_nodes < self.initial_nodes:
            raise ValueError("maximum_nodes must be greater than or equal to initial_nodes")
        return self


class MaterialSI(FrozenModel):
    """Material parameters resolved to SI units."""

    relative_permittivity: PositiveFloat
    reported_total_fe_m3: PositiveFloat
    oxygen_vacancy_m3: PositiveFloat
    electron_effective_mass_kg: PositiveFloat
    hole_effective_mass_kg: PositiveFloat
    electron_mobility_m2_per_v_s: PositiveFloat
    hole_mobility_m2_per_v_s: PositiveFloat
    oxygen_vacancy_mobility_m2_per_v_s: PositiveFloat
    recombination_m3_per_s: NonNegativeFloat


class ExperimentSI(FrozenModel):
    """Experimental conditions resolved to SI units."""

    temperature_k: PositiveFloat
    voltage_v: PositiveFloat
    thickness_m: PositiveFloat


class PreparationSI(FrozenModel):
    """Sample preparation conditions resolved to SI units."""

    annealing_temperature_k: PositiveFloat
    oxygen_partial_pressure_pa: PositiveFloat


class BarrierCaseSI(FrozenModel):
    """One contact-barrier case resolved to joules."""

    name: str
    electron_barrier_j: NonNegativeFloat
    hole_barrier_j: NonNegativeFloat


class BackgroundSI(FrozenModel):
    """Resolved equilibrium background concentrations in SI units."""

    model: BackgroundModel
    electron_m3: NonNegativeFloat
    hole_m3: NonNegativeFloat
    charged_fe3_m3: NonNegativeFloat | None
    neutral_fe4_m3: NonNegativeFloat | None


class ResolvedInputsSI(FrozenModel):
    """Complete immutable numerical input snapshot in SI units."""

    material: MaterialSI
    experiment: ExperimentSI
    preparation: PreparationSI
    background: BackgroundSI
    cases: tuple[BarrierCaseSI, ...]
    solver: SolverSettings
    model_assumptions: tuple[str, ...]


class InputConfig(FrozenModel):
    """Complete validated configuration in literature-facing units."""

    material: MaterialInputs
    experiment: ExperimentInputs
    preparation: PreparationInputs
    background: BackgroundInputs = BackgroundInputs()
    cases: tuple[BarrierCaseInputs, ...]
    solver: SolverSettings

    @model_validator(mode="after")
    def validate_cases(self) -> Self:
        if not self.cases:
            raise ValueError("at least one barrier case is required")
        names = [case.name for case in self.cases]
        if len(names) != len(set(names)):
            raise ValueError("barrier case names must be unique")
        return self

    def to_si(self) -> ResolvedInputsSI:
        """Convert literature-facing values once at the numerical boundary."""
        material = self.material
        experiment = self.experiment
        total_fe_m3 = material.reported_total_fe_cm3 * 1e6
        preparation = PreparationSI(
            annealing_temperature_k=self.preparation.annealing_temperature_k,
            oxygen_partial_pressure_pa=self.preparation.oxygen_partial_pressure_bar * 1e5,
        )
        if self.background.model is BackgroundModel.PREPARATION_EQUILIBRIUM:
            annealed = annealed_fe_equilibrium(
                preparation.annealing_temperature_k,
                preparation.oxygen_partial_pressure_pa,
                total_fe_m3,
            )
            oxygen_vacancy_m3 = annealed.oxygen_vacancy_m3
        else:
            oxygen_vacancy_m3 = material.oxygen_vacancy_cm3 * 1e6
        if self.background.model in (
            BackgroundModel.QUENCHED_EQUILIBRIUM,
            BackgroundModel.PREPARATION_EQUILIBRIUM,
        ):
            equilibrium = quenched_fe_equilibrium(
                experiment.temperature_k,
                total_fe_m3,
                oxygen_vacancy_m3,
            )
            background = BackgroundSI(
                model=self.background.model,
                electron_m3=equilibrium.electron_m3,
                hole_m3=equilibrium.hole_m3,
                charged_fe3_m3=equilibrium.charged_fe3_m3,
                neutral_fe4_m3=equilibrium.neutral_fe4_m3,
            )
            if self.background.model is BackgroundModel.PREPARATION_EQUILIBRIUM:
                background_assumptions = (
                    "oxygen-vacancy density follows Denk equilibrium at the specified "
                    "annealing temperature and oxygen partial pressure",
                    "the annealed oxygen-vacancy density is frozen while electronic and "
                    "Fe equilibria re-establish at the simulation temperature",
                )
            else:
                background_assumptions = (
                    "equilibrium carriers and Fe charge states follow the quenched Denk "
                    "defect chemistry",
                    "specified oxygen-vacancy density is frozen while electronic and Fe "
                    "equilibria re-establish",
                )
        else:
            background = BackgroundSI(
                model=self.background.model,
                electron_m3=0.0,
                hole_m3=0.0,
                charged_fe3_m3=None,
                neutral_fe4_m3=None,
            )
            background_assumptions = (
                "negligible equilibrium electron and hole densities (n0 = p0 = 0)",
                "fixed compensating charge equivalent to twice the oxygen-vacancy density",
            )
        return ResolvedInputsSI(
            material=MaterialSI(
                relative_permittivity=material.relative_permittivity,
                reported_total_fe_m3=total_fe_m3,
                oxygen_vacancy_m3=oxygen_vacancy_m3,
                electron_effective_mass_kg=(material.electron_effective_mass_m0 * electron_mass),
                hole_effective_mass_kg=material.hole_effective_mass_m0 * electron_mass,
                electron_mobility_m2_per_v_s=(material.electron_mobility_cm2_per_v_s * 1e-4),
                hole_mobility_m2_per_v_s=material.hole_mobility_cm2_per_v_s * 1e-4,
                oxygen_vacancy_mobility_m2_per_v_s=(
                    material.oxygen_vacancy_mobility_cm2_per_v_s * 1e-4
                ),
                recombination_m3_per_s=material.recombination_cm3_per_s * 1e-6,
            ),
            experiment=ExperimentSI(
                temperature_k=experiment.temperature_k,
                voltage_v=experiment.voltage_v,
                thickness_m=experiment.thickness_um * 1e-6,
            ),
            preparation=preparation,
            background=background,
            cases=tuple(
                BarrierCaseSI(
                    name=case.name,
                    electron_barrier_j=case.electron_barrier_ev * elementary_charge,
                    hole_barrier_j=case.hole_barrier_ev * elementary_charge,
                )
                for case in self.cases
            ),
            solver=self.solver,
            model_assumptions=COMMON_MODEL_ASSUMPTIONS + background_assumptions,
        )


def load_config(path: str | Path) -> InputConfig:
    """Load and validate a TOML input file."""
    config_path = Path(path)
    with config_path.open("rb") as config_file:
        data = tomllib.load(config_file)
    return InputConfig.model_validate(data)


def load_inputs(path: str | Path) -> InputConfig:
    """Load validated literature-facing inputs from TOML or a saved JSON snapshot."""
    input_path = Path(path)
    if input_path.suffix.lower() == ".toml":
        return load_config(input_path)
    if input_path.suffix.lower() == ".json":
        return InputConfig.model_validate_json(input_path.read_text(encoding="utf-8"))
    raise ValueError(f"unsupported input format {input_path.suffix!r}; expected .toml or .json")
