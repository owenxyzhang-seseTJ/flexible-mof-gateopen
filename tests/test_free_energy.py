from __future__ import annotations

import math

from flexible_mof_gateopen.free_energy import (
    co2_fugacity_coefficient,
    cumulative_integral,
    integrate_required_host_penalty,
    solve_crossing_pressure,
)


def test_henry_limit_is_exact_below_first_point() -> None:
    value = cumulative_integral([1.0, 10.0], [2.0, 8.0], target_pressure_bar=0.25)
    assert math.isclose(value, 0.5, rel_tol=0.0, abs_tol=1e-12)


def test_formula_unit_normalization_and_error() -> None:
    value, error = integrate_required_host_penalty(
        [1.0, 2.0, 4.0],
        [2.0, 4.0, 8.0],
        [0.0, 0.0, 0.0],
        temperature_K=300.0,
        open_formula_units=2.0,
        closed_formula_units=1.0,
        open_error=[0.2, 0.2, 0.2],
        closed_error=[0.1, 0.1, 0.1],
        target_pressure_bar=1.0,
    )
    assert value > 0
    assert error is not None and error > 0


def test_crossing_is_interpolated_in_log_pressure() -> None:
    crossing = solve_crossing_pressure([1.0, 10.0], [0.0, 10.0], 5.0)
    assert crossing is not None
    assert math.isclose(crossing, math.sqrt(10.0), rel_tol=1e-12)


def test_co2_eos_sensitivity_is_finite() -> None:
    assert co2_fugacity_coefficient(19.0, 268.0, "ideal") == 1.0
    assert 0.7 < co2_fugacity_coefficient(19.0, 268.0, "pr") < 1.0
    assert 0.7 < co2_fugacity_coefficient(19.0, 268.0, "srk") < 1.0
