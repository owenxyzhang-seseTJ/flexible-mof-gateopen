#!/usr/bin/env python3
"""Parse final absolute/excess loading and move diagnostics from RASPA text."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_loading(text: str, label: str) -> tuple[float, float]:
    pattern = re.compile(rf"{re.escape(label)} loading average\s+([+\-0-9.eE]+)\s+\+/-\s+([+\-0-9.eE]+)\s+\[molecules/cell\]")
    matches = pattern.findall(text)
    if not matches:
        raise ValueError(f"missing {label} loading")
    return tuple(float(value) for value in matches[-1])


def parse_move_stat(text: str, move: str, statistic: str) -> list[int] | None:
    """Parse a numbered statistic from a RASPA move-diagnostics line."""
    pattern = re.compile(
        rf"^\s*{re.escape(move)}\s+{re.escape(statistic)}:\s+(.+)$",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if not match:
        return None
    values = []
    for token in match.group(1).split():
        try:
            values.append(int(token))
        except ValueError:
            break
    return values or None


def parse_output(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    absolute, absolute_error = parse_loading(text, "Abs.")
    excess, excess_error = parse_loading(text, "Excess")
    return {
        "path": str(path),
        "absolute_loading_molecules_per_cell": absolute,
        "absolute_loading_error_molecules_per_cell": absolute_error,
        "excess_loading_molecules_per_cell": excess,
        "excess_loading_error_molecules_per_cell": excess_error,
        "move_diagnostic": {
            "translation_total": parse_move_stat(text, "Translation", "total"),
            "translation_constructed": parse_move_stat(text, "Translation", "constructed"),
            "swap_constructed": parse_move_stat(text, "Swap (CBMC)", "constructed"),
            "swap_accepted": parse_move_stat(text, "Swap (CBMC)", "accepted"),
            "widom_constructed": parse_move_stat(text, "Widom", "constructed"),
        },
        "warning": "zero absolute loading requires an accessibility/setup audit using constructed and accepted guest-move diagnostics; it is not a precise zero-uncertainty measurement",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(parse_output(args.output), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
