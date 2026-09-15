"""Validated configuration loading and conversion to SI units."""

import tomllib
from pathlib import Path
from typing import Annotated, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator
from scipy.constants import electron_mass, elementary_charge

MODEL_ASSUMPTIONS = (
    "homogeneous one-dimensional isothermal drift-only transport",
    "uniform prescribed oxygen-vacancy density with finite drift conductivity",
    "negligible equilibrium electron and hole densities (n0 = p0 = 0)",
    "fixed compensating charge equivalent to twice the oxygen-vacancy density",
    "constant bimolecular recombination coefficient",
    "no diffusion, field-dependent transport, trap kinetics, or vacancy evolution",
)

PositiveFloat = Annotated[float, Field(gt=0, allow_inf_nan=False)]
NonNegativeFloat = Annotated[float, Field(ge=0, allow_inf_nan=False)]
PositiveInt = Annotated[int, Field(gt=0)]


class FrozenModel(BaseModel):
    """Base for immutable models with a closed configuration schema."""

    model_config = ConfigDict(frozen=True, extra="forbid")


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


class BarrierCaseSI(FrozenModel):
    """One contact-barrier case resolved to joules."""

    name: str
    electron_barrier_j: NonNegativeFloat
    hole_barrier_j: NonNegativeFloat


class ResolvedInputsSI(FrozenModel):
    """Complete immutable numerical input snapshot in SI units."""

    material: MaterialSI
    experiment: ExperimentSI
    cases: tuple[BarrierCaseSI, ...]
    solver: SolverSettings
    model_assumptions: tuple[str, ...]


class InputConfig(FrozenModel):
    """Complete validated configuration in literature-facing units."""

    material: MaterialInputs
    experiment: ExperimentInputs
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
        return ResolvedInputsSI(
            material=MaterialSI(
                relative_permittivity=material.relative_permittivity,
                reported_total_fe_m3=material.reported_total_fe_cm3 * 1e6,
                oxygen_vacancy_m3=material.oxygen_vacancy_cm3 * 1e6,
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
            cases=tuple(
                BarrierCaseSI(
                    name=case.name,
                    electron_barrier_j=case.electron_barrier_ev * elementary_charge,
                    hole_barrier_j=case.hole_barrier_ev * elementary_charge,
                )
                for case in self.cases
            ),
            solver=self.solver,
            model_assumptions=MODEL_ASSUMPTIONS,
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
