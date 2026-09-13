#!/usr/bin/env python3
"""Run a formula-unit-normalized two-state fixed-host free-energy screen.

Example:
  python scripts/quick_free_energy_screen.py \
    --csv examples/xpcu5zn_co2_268K/data/fixed_host_isotherms_268K.csv \
    --temperature 268 --open-phase open_alpha --closed-phase closed_gamma \
    --open-formula-units 2 --closed-formula-units 4 \
    --target-pressure 19 --host-penalty 67.9 --eos ideal
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from flexible_mof_gateopen.free_energy import (
        ScreeningResult,
        co2_fugacity_coefficient,
        integrate_required_host_penalty,
        required_penalty_curve,
        solve_crossing_pressure,
    )
except ModuleNotFoundError:  # standalone copy inside a skill directory
    from free_energy import (  # type: ignore[no-redef]
        ScreeningResult,
        co2_fugacity_coefficient,
        integrate_required_host_penalty,
        required_penalty_curve,
        solve_crossing_pressure,
    )


def read_rows(path: Path, open_phase: str, closed_phase: str, refine_path: Path | None = None, replace_below: float = 1.0) -> tuple[list[float], list[float], list[float], list[float], list[float]]:
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    if refine_path is not None:
        refined = list(csv.DictReader(refine_path.open(newline="", encoding="utf-8")))
        rows = [row for row in rows if not (row.get("phase") in {open_phase, closed_phase} and float(row["pressure_bar"]) <= replace_below)]
        rows.extend(row for row in refined if row.get("phase") in {open_phase, closed_phase} and float(row["pressure_bar"]) <= replace_below)
    selected = [row for row in rows if row.get("phase") in {open_phase, closed_phase}]
    if not selected:
        raise ValueError(f"no rows for {open_phase!r}/{closed_phase!r} in {path}")
    grouped: dict[str, dict[float, dict[str, str]]] = {open_phase: {}, closed_phase: {}}
    for row in selected:
        pressure = float(row["pressure_bar"])
        if pressure in grouped[row["phase"]]:
            raise ValueError(f"duplicate pressure {pressure:g} for phase {row['phase']}")
        grouped[row["phase"]][pressure] = row
    if set(grouped[open_phase]) != set(grouped[closed_phase]):
        raise ValueError("open and closed phases do not share a pressure grid")
    pressures = sorted(grouped[open_phase])

    def value(row: dict[str, str], preferred: tuple[str, ...]) -> float:
        for key in preferred:
            if key in row:
                return float(row[key])
        raise ValueError(f"missing one of columns {preferred}")

    open_loading = [value(grouped[open_phase][p], ("absolute_loading_molecules_per_cell", "absolute_loading")) for p in pressures]
    closed_loading = [value(grouped[closed_phase][p], ("absolute_loading_molecules_per_cell", "absolute_loading")) for p in pressures]
    open_error = [value(grouped[open_phase][p], ("absolute_loading_error_molecules_per_cell", "absolute_loading_error")) for p in pressures]
    closed_error = [value(grouped[closed_phase][p], ("absolute_loading_error_molecules_per_cell", "absolute_loading_error")) for p in pressures]
    return pressures, open_loading, closed_loading, open_error, closed_error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--temperature", type=float, required=True, dest="temperature_K")
    parser.add_argument("--open-phase", default="open_alpha")
    parser.add_argument("--closed-phase", default="closed_gamma")
    parser.add_argument("--open-formula-units", type=float, default=1.0)
    parser.add_argument("--closed-formula-units", type=float, default=1.0)
    parser.add_argument("--target-pressure", type=float, required=True, dest="target_pressure_bar")
    parser.add_argument("--host-penalty", type=float, default=None, dest="host_penalty")
    parser.add_argument("--eos", choices=("ideal", "pr", "srk"), default="ideal")
    parser.add_argument("--refine-csv", type=Path, default=None, help="optional independent low-pressure CSV replacing points below --replace-below")
    parser.add_argument("--replace-below", type=float, default=1.0)
    parser.add_argument("--out", type=Path, default=None)
    return parser


def trim_nonmonotonic_eos_grid(
    pressures: list[float],
    open_loading: list[float],
    closed_loading: list[float],
    open_error: list[float],
    closed_error: list[float],
    *,
    temperature_K: float,
    eos: str,
    target_pressure_bar: float,
) -> tuple[list[float], list[float], list[float], list[float], list[float], bool]:
    """Trim a gas-EOS grid after a non-monotonic pure-fluid branch.

    Below the CO2 saturation region the vapor-root fugacity is monotonic.  A
    supplied grid can extend into a two-phase region (notably for SRK at 268 K
   ); those points are not silently integrated as if they were one gas phase.
    """
    if eos == "ideal":
        return pressures, open_loading, closed_loading, open_error, closed_error, False
    keep = len(pressures)
    while keep >= 2:
        fugacities = [p * co2_fugacity_coefficient(p, temperature_K, eos) for p in pressures[:keep]]
        if all(b > a for a, b in zip(fugacities, fugacities[1:])):
            break
        keep -= 1
    if keep < 2 or target_pressure_bar > pressures[keep - 1]:
        raise ValueError(
            f"{eos.upper()} fugacity is not monotonic through the target pressure; "
            "restrict the grid to one fluid branch or use ideal fugacity"
        )
    return (
        pressures[:keep],
        open_loading[:keep],
        closed_loading[:keep],
        open_error[:keep],
        closed_error[:keep],
        keep != len(pressures),
    )


def main() -> None:
    args = build_parser().parse_args()
    pressures, open_loading, closed_loading, open_error, closed_error = read_rows(args.csv, args.open_phase, args.closed_phase, args.refine_csv, args.replace_below)
    pressures, open_loading, closed_loading, open_error, closed_error, grid_trimmed = trim_nonmonotonic_eos_grid(
        pressures,
        open_loading,
        closed_loading,
        open_error,
        closed_error,
        temperature_K=args.temperature_K,
        eos=args.eos,
        target_pressure_bar=args.target_pressure_bar,
    )
    fugacity = lambda p: p * co2_fugacity_coefficient(p, args.temperature_K, args.eos)
    value, error = integrate_required_host_penalty(
        pressures,
        open_loading,
        closed_loading,
        temperature_K=args.temperature_K,
        open_formula_units=args.open_formula_units,
        closed_formula_units=args.closed_formula_units,
        open_error=open_error,
        closed_error=closed_error,
        fugacity=fugacity,
        target_pressure_bar=args.target_pressure_bar,
    )
    curve = required_penalty_curve(
        pressures,
        open_loading,
        closed_loading,
        temperature_K=args.temperature_K,
        open_formula_units=args.open_formula_units,
        closed_formula_units=args.closed_formula_units,
        fugacity=fugacity,
    )
    crossing = None
    if args.host_penalty is not None:
        crossing = solve_crossing_pressure(pressures, curve, args.host_penalty)
    result = ScreeningResult(
        temperature_K=args.temperature_K,
        target_pressure_bar=args.target_pressure_bar,
        fugacity_model=args.eos,
        required_host_penalty_kj_per_formula_unit=value,
        required_host_penalty_error_kj_per_formula_unit=error,
        assumed_host_penalty_kj_per_formula_unit=args.host_penalty,
        crossing_pressure_bar=crossing,
    )
    payload = {
        "result": result.__dict__,
        "pressure_bar": pressures,
        "required_host_penalty_curve_kJ_per_formula_unit": curve,
        "source_csv": str(args.csv),
        "refine_csv": str(args.refine_csv) if args.refine_csv else None,
        "pressure_grid_trimmed_for_eos_branch": grid_trimmed,
        "limitations": [
            "fixed-host adsorption grand-potential integration only",
            "host penalty is an explicit input or calibration, not supplied by PACMOF charges",
            "a crossing is an equilibrium scenario and not P_open/P_close kinetics",
            "spinodals and barriers require intermediate framework states or enhanced sampling",
        ],
    }
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        args.out.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
