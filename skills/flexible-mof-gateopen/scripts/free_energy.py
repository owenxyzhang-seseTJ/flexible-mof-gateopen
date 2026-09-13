"""Free-energy bookkeeping for fixed-host two-state screening.

The functions here intentionally implement the smallest defensible model used
by the example workflow.  They integrate the *adsorption* grand-potential
difference from absolute loadings and leave the empty-host free-energy term as
an explicit input.  They do not claim to sample a flexible framework or a
kinetic gate-crossing barrier.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Iterable, Sequence

R_KJ_MOL_K = 8.314462618e-3
FugacityFunction = Callable[[float], float]


def _as_float_list(values: Iterable[float], name: str) -> list[float]:
    result = [float(value) for value in values]
    if not result:
        raise ValueError(f"{name} cannot be empty")
    if any(not math.isfinite(value) for value in result):
        raise ValueError(f"{name} contains a non-finite value")
    return result


def _validate_grid(pressure_bar: Sequence[float], delta_loading: Sequence[float]) -> None:
    if len(pressure_bar) != len(delta_loading):
        raise ValueError("pressure and loading arrays must have the same length")
    if len(pressure_bar) < 2:
        raise ValueError("at least two pressure points are required")
    if any(p <= 0 for p in pressure_bar):
        raise ValueError("pressure values must be positive")
    if any(p2 <= p1 for p1, p2 in zip(pressure_bar, pressure_bar[1:])):
        raise ValueError("pressure values must be strictly increasing")


def _fugacity_grid(pressure_bar: Sequence[float], fugacity: FugacityFunction | None) -> list[float]:
    values = [float(p if fugacity is None else fugacity(float(p))) for p in pressure_bar]
    if any(f <= 0 or not math.isfinite(f) for f in values):
        raise ValueError("fugacity values must be finite and positive")
    if any(f2 <= f1 for f1, f2 in zip(values, values[1:])):
        raise ValueError("fugacity must increase with pressure")
    return values


def _interp_logx(x: Sequence[float], y: Sequence[float], target: float) -> float:
    """Linear interpolation in log(x), with endpoint extrapolation disabled."""
    if target < x[0] or target > x[-1]:
        raise ValueError("target is outside the supplied grid")
    if target == x[0]:
        return float(y[0])
    for i in range(1, len(x)):
        if target <= x[i]:
            x0, x1 = math.log(x[i - 1]), math.log(x[i])
            fraction = (math.log(target) - x0) / (x1 - x0)
            return float(y[i - 1] + fraction * (y[i] - y[i - 1]))
    return float(y[-1])


def cumulative_integral(
    pressure_bar: Sequence[float],
    delta_loading: Sequence[float],
    *,
    fugacity: FugacityFunction | None = None,
    target_pressure_bar: float | None = None,
) -> float:
    """Return ``integral_0^f ΔN(f) dln(f)`` in molecules per formula unit.

    The first measured point is treated with the Henry limit ``N proportional
    to f``.  Above that point, ``ΔN`` is linearly interpolated in log fugacity
    and integrated segment by segment.  This avoids the common, biased
    practice of adding an artificial zero-pressure point and applying a
    trapezoid in log pressure.
    """
    p = _as_float_list(pressure_bar, "pressure_bar")
    n = _as_float_list(delta_loading, "delta_loading")
    _validate_grid(p, n)
    f = _fugacity_grid(p, fugacity)
    target = p[-1] if target_pressure_bar is None else float(target_pressure_bar)
    if target <= 0:
        raise ValueError("target pressure must be positive")
    ft = float(target if fugacity is None else fugacity(target))
    if ft <= 0 or not math.isfinite(ft):
        raise ValueError("target fugacity must be finite and positive")

    if ft <= f[0]:
        return n[0] * ft / f[0]

    integral = n[0]  # exact Henry contribution from zero to f[0]
    last_f, last_n = f[0], n[0]
    for next_f, next_n in zip(f[1:], n[1:]):
        if ft <= next_f:
            current_n = _interp_logx([last_f, next_f], [last_n, next_n], ft)
            return integral + 0.5 * (last_n + current_n) * math.log(ft / last_f)
        integral += 0.5 * (last_n + next_n) * math.log(next_f / last_f)
        last_f, last_n = next_f, next_n
    return integral


def cumulative_error(
    pressure_bar: Sequence[float],
    delta_error: Sequence[float],
    *,
    fugacity: FugacityFunction | None = None,
    target_pressure_bar: float | None = None,
) -> float:
    """Conservative independent-error propagation for :func:`cumulative_integral`."""
    p = _as_float_list(pressure_bar, "pressure_bar")
    e = _as_float_list(delta_error, "delta_error")
    _validate_grid(p, e)
    f = _fugacity_grid(p, fugacity)
    target = p[-1] if target_pressure_bar is None else float(target_pressure_bar)
    ft = float(target if fugacity is None else fugacity(target))
    if ft <= f[0]:
        return e[0] * ft / f[0]
    variance = e[0] ** 2  # exact Henry contribution
    last_f, last_e = f[0], e[0]
    for next_f, next_e in zip(f[1:], e[1:]):
        if ft <= next_f:
            dx = math.log(ft / last_f)
            fraction = math.log(ft / last_f) / math.log(next_f / last_f)
            current_e = last_e + fraction * (next_e - last_e)
            variance += (0.5 * dx) ** 2 * (last_e**2 + current_e**2)
            return math.sqrt(variance)
        dx = math.log(next_f / last_f)
        variance += (0.5 * dx) ** 2 * (last_e**2 + next_e**2)
        last_f, last_e = next_f, next_e
    return math.sqrt(variance)


def integrate_required_host_penalty(
    pressure_bar: Sequence[float],
    open_loading: Sequence[float],
    closed_loading: Sequence[float],
    *,
    temperature_K: float,
    open_formula_units: float = 1.0,
    closed_formula_units: float = 1.0,
    open_error: Sequence[float] | None = None,
    closed_error: Sequence[float] | None = None,
    fugacity: FugacityFunction | None = None,
    target_pressure_bar: float | None = None,
) -> tuple[float, float | None]:
    """Return the host penalty required for a two-state crossing.

    Loadings are supplied per crystallographic cell and normalized to one
    framework formula unit before integration.  A positive result means that
    the open state needs that much larger empty-host free energy to be at
    equilibrium at the target pressure.
    """
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if open_formula_units <= 0 or closed_formula_units <= 0:
        raise ValueError("formula-unit counts must be positive")
    p = _as_float_list(pressure_bar, "pressure_bar")
    op = [float(v) / open_formula_units for v in open_loading]
    cl = [float(v) / closed_formula_units for v in closed_loading]
    if len(op) != len(cl):
        raise ValueError("open and closed loading arrays must have the same length")
    delta = [a - b for a, b in zip(op, cl)]
    error = None
    if open_error is not None or closed_error is not None:
        if open_error is None or closed_error is None:
            raise ValueError("open_error and closed_error must be supplied together")
        oe = [float(v) / open_formula_units for v in open_error]
        ce = [float(v) / closed_formula_units for v in closed_error]
        error = [math.hypot(a, b) for a, b in zip(oe, ce)]
    factor = R_KJ_MOL_K * float(temperature_K)
    value = factor * cumulative_integral(p, delta, fugacity=fugacity, target_pressure_bar=target_pressure_bar)
    uncertainty = None
    if error is not None:
        uncertainty = factor * cumulative_error(p, error, fugacity=fugacity, target_pressure_bar=target_pressure_bar)
    return value, uncertainty


def required_penalty_curve(
    pressure_bar: Sequence[float],
    open_loading: Sequence[float],
    closed_loading: Sequence[float],
    *,
    temperature_K: float,
    open_formula_units: float = 1.0,
    closed_formula_units: float = 1.0,
    fugacity: FugacityFunction | None = None,
) -> list[float]:
    """Return the cumulative required host penalty at each pressure point."""
    return [
        integrate_required_host_penalty(
            pressure_bar,
            open_loading,
            closed_loading,
            temperature_K=temperature_K,
            open_formula_units=open_formula_units,
            closed_formula_units=closed_formula_units,
            fugacity=fugacity,
            target_pressure_bar=p,
        )[0]
        for p in pressure_bar
    ]


def solve_crossing_pressure(
    pressure_bar: Sequence[float],
    required_penalty_kj_per_formula_unit: Sequence[float],
    host_penalty_kj_per_formula_unit: float,
) -> float | None:
    """Find a log-pressure crossing of a host penalty and a cumulative curve."""
    p = _as_float_list(pressure_bar, "pressure_bar")
    y = _as_float_list(required_penalty_kj_per_formula_unit, "required_penalty")
    _validate_grid(p, y)
    target = float(host_penalty_kj_per_formula_unit)
    for p0, p1, y0, y1 in zip(p, p[1:], y, y[1:]):
        if y0 == target:
            return p0
        if (y0 - target) * (y1 - target) <= 0:
            if y1 == y0:
                return p0
            fraction = (target - y0) / (y1 - y0)
            return math.exp(math.log(p0) + fraction * (math.log(p1) - math.log(p0)))
    if y[-1] == target:
        return p[-1]
    return None


def _cubic_roots_real(coefficients: Sequence[float]) -> list[float]:
    """Small dependency-free cubic solver for EOS fugacity calculations."""
    # Cardano implementation for a monic cubic x^3 + a x^2 + b x + c.
    a, b, c = [float(v) for v in coefficients[1:]]
    p = b - a * a / 3.0
    q = 2.0 * a**3 / 27.0 - a * b / 3.0 + c
    discriminant = (q / 2.0) ** 2 + (p / 3.0) ** 3
    if discriminant >= 0:
        root = math.copysign(abs(-q / 2.0 + math.sqrt(discriminant)) ** (1 / 3), -q / 2.0 + math.sqrt(discriminant))
        other = math.copysign(abs(-q / 2.0 - math.sqrt(discriminant)) ** (1 / 3), -q / 2.0 - math.sqrt(discriminant))
        return [root + other - a / 3.0]
    radius = 2.0 * math.sqrt(-p / 3.0)
    angle = math.acos(max(-1.0, min(1.0, (3.0 * q / (2.0 * p)) * math.sqrt(-3.0 / p)))) / 3.0
    return [radius * math.cos(angle - 2.0 * math.pi * k / 3.0) - a / 3.0 for k in range(3)]


def co2_fugacity_coefficient(pressure_bar: float, temperature_K: float, model: str = "ideal") -> float:
    """Return a pure-CO2 fugacity coefficient from ideal/PR/SRK EOS.

    This is a post-processing sensitivity helper.  It is not a replacement for
    running RASPA with a validated reservoir fugacity.
    """
    model = model.lower()
    if model == "ideal":
        return 1.0
    if pressure_bar <= 0 or temperature_K <= 0:
        raise ValueError("pressure and temperature must be positive")
    tc, pc, omega = 304.1282, 7377300.0, 0.22394
    r = 8.314462618
    p = pressure_bar * 1.0e5
    if model == "pr":
        kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
        alpha = (1.0 + kappa * (1.0 - math.sqrt(temperature_K / tc))) ** 2
        a = 0.45724 * r**2 * tc**2 * alpha / pc
        b = 0.07780 * r * tc / pc
        aa = a * p / (r**2 * temperature_K**2)
        bb = b * p / (r * temperature_K)
        roots = _cubic_roots_real([1.0, -(1.0 - bb), aa - 3 * bb**2 - 2 * bb, -(aa * bb - bb**2 - bb**3)])
        z = max(root for root in roots if root > bb)
        return math.exp(z - 1.0 - math.log(z - bb) - aa / (2.0 * math.sqrt(2.0) * bb) * math.log((z + (1 + math.sqrt(2.0)) * bb) / (z + (1 - math.sqrt(2.0)) * bb)))
    if model == "srk":
        m = 0.480 + 1.574 * omega - 0.176 * omega**2
        alpha = (1.0 + m * (1.0 - math.sqrt(temperature_K / tc))) ** 2
        a = 0.42747 * r**2 * tc**2 * alpha / pc
        b = 0.08664 * r * tc / pc
        aa = a * p / (r**2 * temperature_K**2)
        bb = b * p / (r * temperature_K)
        roots = _cubic_roots_real([1.0, -1.0, aa - bb - bb**2, -aa * bb])
        z = max(root for root in roots if root > bb)
        return math.exp(z - 1.0 - math.log(z - bb) - aa / bb * math.log(1.0 + bb / z))
    raise ValueError(f"unknown fugacity model: {model}")


@dataclass(frozen=True)
class ScreeningResult:
    """Machine-readable summary for a quick fixed-host screen."""

    temperature_K: float
    target_pressure_bar: float
    fugacity_model: str
    required_host_penalty_kj_per_formula_unit: float
    required_host_penalty_error_kj_per_formula_unit: float | None
    assumed_host_penalty_kj_per_formula_unit: float | None
    crossing_pressure_bar: float | None
    classification: str = "conditional_fixed_host_screen"
