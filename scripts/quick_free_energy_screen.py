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
        integrate_required_host_penalty,
        required_penalty_curve,
        solve_crossing_pressure,
    )
    from flexible_mof_gateopen.preflight import load_gas_profile, pressure_gate
except ModuleNotFoundError:  # standalone copy inside a skill directory
    from free_energy import (  # type: ignore[no-redef]
        ScreeningResult,
        integrate_required_host_penalty,
        required_penalty_curve,
        solve_crossing_pressure,
    )
    from preflight import load_gas_profile, pressure_gate


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
    parser.add_argument("--gas-profile", type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    pressures, open_loading, closed_loading, open_error, closed_error = read_rows(args.csv, args.open_phase, args.closed_phase, args.refine_csv, args.replace_below)
    profile = load_gas_profile(args.gas_profile)
    reservoir = pressure_gate(profile, args.temperature_K, pressures)
    pressure_gate(profile, args.temperature_K, [args.target_pressure_bar])
    if args.eos != "ideal":
        raise ValueError("PR/SRK reweighting of loadings sampled at ideal fugacity is not a consistent reservoir correction. Use actual sampled fugacities; legacy EOS scenarios are withdrawn.")
    fugacity = None  # These input isotherms were sampled with f=P.
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
        "reservoir": reservoir,
        "pressure_bar": pressures,
        "required_host_penalty_curve_kJ_per_formula_unit": curve,
        "source_csv": str(args.csv),
        "refine_csv": str(args.refine_csv) if args.refine_csv else None,
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
