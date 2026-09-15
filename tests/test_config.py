from pathlib import Path

import pytest
from pydantic import ValidationError
from scipy.constants import electron_mass, elementary_charge

from charge_injection_sim.config import BackgroundModel, InputConfig, load_config, load_inputs

BENCHMARK_PATH = Path(__file__).parents[1] / "configs" / "figure4b.toml"


def test_benchmark_loads_with_expected_cases() -> None:
    config = load_config(BENCHMARK_PATH)

    assert config.experiment.temperature_k == pytest.approx(483.15)
    assert [case.electron_barrier_ev for case in config.cases] == [1.0, 0.85, 0.8]
    assert config.solver.maximum_nodes == 20_000
    assert config.background.model is BackgroundModel.SIMPLIFIED


def test_literature_units_convert_to_si() -> None:
    resolved = load_config(BENCHMARK_PATH).to_si()

    assert resolved.material.reported_total_fe_m3 == pytest.approx(5.58e24)
    assert resolved.material.oxygen_vacancy_m3 == pytest.approx(2.43e24)
    assert resolved.material.electron_effective_mass_kg == pytest.approx(6 * electron_mass)
    assert resolved.material.hole_effective_mass_kg == pytest.approx(12 * electron_mass)
    assert resolved.material.electron_mobility_m2_per_v_s == pytest.approx(0.56e-4)
    assert resolved.material.hole_mobility_m2_per_v_s == pytest.approx(0.41e-4)
    assert resolved.material.oxygen_vacancy_mobility_m2_per_v_s == pytest.approx(2.20e-12)
    assert resolved.material.recombination_m3_per_s == pytest.approx(1.0e-14)
    assert resolved.experiment.thickness_m == pytest.approx(500e-6)
    assert resolved.cases[0].electron_barrier_j == pytest.approx(elementary_charge)


def test_resolved_inputs_serialize_to_json() -> None:
    resolved = load_config(BENCHMARK_PATH).to_si()

    restored = type(resolved).model_validate_json(resolved.model_dump_json())

    assert restored == resolved
    assert restored.model_assumptions


def test_quenched_equilibrium_background_matches_benchmark_defect_chemistry() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["background"]["model"] = "quenched_equilibrium"

    resolved = InputConfig.model_validate(data).to_si()
    background = resolved.background

    assert background.electron_m3 == pytest.approx(2.074e7, rel=5e-4)
    assert background.hole_m3 == pytest.approx(1.4764e16, rel=5e-4)
    assert background.charged_fe3_m3 == pytest.approx(4.86e24, rel=5e-8)
    assert background.neutral_fe4_m3 == pytest.approx(7.20e23, rel=5e-8)
    assert background.charged_fe3_m3 is not None
    charge_residual = (
        2.0 * resolved.material.oxygen_vacancy_m3
        + background.hole_m3
        - background.charged_fe3_m3
        - background.electron_m3
    )
    assert abs(charge_residual) / resolved.material.reported_total_fe_m3 < 1e-12


def test_quenched_equilibrium_preserves_small_fe4_population() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["background"]["model"] = "quenched_equilibrium"
    data["material"]["reported_total_fe_cm3"] = 1e18

    resolved = InputConfig.model_validate(data).to_si()
    background = resolved.background

    assert background.charged_fe3_m3 is not None
    assert background.neutral_fe4_m3 is not None
    assert background.neutral_fe4_m3 > 0
    assert background.charged_fe3_m3 + background.neutral_fe4_m3 == pytest.approx(
        resolved.material.reported_total_fe_m3, rel=1e-14
    )
    assert background.electron_m3 > background.hole_m3


def test_unknown_background_model_is_rejected() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["background"]["model"] = "fully_ionized_fe"

    with pytest.raises(ValidationError):
        InputConfig.model_validate(data)


def test_missing_background_uses_simplified_model() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    del data["background"]

    config = InputConfig.model_validate(data)

    assert config.background.model is BackgroundModel.SIMPLIFIED


def test_load_inputs_rejects_unknown_file_type(tmp_path: Path) -> None:
    path = tmp_path / "inputs.yaml"

    with pytest.raises(ValueError, match=r"expected \.toml or \.json"):
        load_inputs(path)


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("material", "relative_permittivity", 0),
        ("material", "oxygen_vacancy_cm3", -1),
        ("material", "electron_mobility_cm2_per_v_s", float("nan")),
        ("material", "recombination_cm3_per_s", -1),
        ("experiment", "temperature_k", 0),
        ("experiment", "voltage_v", -40),
        ("experiment", "thickness_um", float("inf")),
        ("solver", "residual_tolerance", 1),
        ("solver", "initial_nodes", 1),
    ],
)
def test_invalid_values_are_rejected(section: str, field: str, value: float) -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data[section][field] = value

    with pytest.raises(ValidationError):
        InputConfig.model_validate(data)


def test_zero_recombination_and_barriers_are_supported() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["material"]["recombination_cm3_per_s"] = 0
    data["cases"][0]["electron_barrier_ev"] = 0

    config = InputConfig.model_validate(data)

    assert config.material.recombination_cm3_per_s == 0
    assert config.cases[0].electron_barrier_ev == 0


def test_unknown_keys_are_rejected() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["material"]["untracked_parameter"] = 1

    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        InputConfig.model_validate(data)


def test_duplicate_case_names_are_rejected() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["cases"][1]["name"] = data["cases"][0]["name"]

    with pytest.raises(ValidationError, match="barrier case names must be unique"):
        InputConfig.model_validate(data)


def test_maximum_nodes_cannot_be_less_than_initial_nodes() -> None:
    data = load_config(BENCHMARK_PATH).model_dump()
    data["solver"]["maximum_nodes"] = data["solver"]["initial_nodes"] - 1

    with pytest.raises(ValidationError, match="maximum_nodes must be"):
        InputConfig.model_validate(data)


def test_models_are_frozen() -> None:
    config = load_config(BENCHMARK_PATH)

    with pytest.raises(ValidationError, match="Instance is frozen"):
        config.experiment.voltage_v = 20.0
