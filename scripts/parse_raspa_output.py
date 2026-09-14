#!/usr/bin/env python3
"""Parse final absolute/excess loading and move diagnostics from RASPA text."""

from __future__ import annotations

import argparse
import json
import re
import math
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


def parse_point(directory, gas, formula_units, unit_cells, expected_seed=None):
    """Parse paired native JSON/text; cell means simulation box, uc means CIF cell.

    Native averageEnergies are total box energies in kelvin (U/kB), not an
    adsorption enthalpy per molecule. Divide by Z * prod(unit_cells) and
    multiply by R to obtain kJ/mol of framework formula units.
    """
    directory = Path(directory)
    texts = list((directory / "output").glob("output*.s0.txt"))
    if len(texts) != 1: raise ValueError("exactly one single-system final output required")
    path=texts[0]
    data=json.loads(path.with_suffix(".json").read_text())
    parsed=parse_output(path)
    text=path.read_text()
    uc=re.findall(r"Abs\. loading average\s+([+\-0-9.eE]+)\s+\+/-\s+([+\-0-9.eE]+)\s+\[molecules/uc\]",text)
    if not uc or not math.isclose(parsed["absolute_loading_molecules_per_cell"]/math.prod(unit_cells),float(uc[-1][0]),rel_tol=2e-6,abs_tol=1e-8):
        raise ValueError("simulation-cell/unit-cell loading normalization mismatch")
    henry=re.findall(r"Average Henry coefficient:\s+([+\-0-9.eE]+)\s+\+/-\s+([+\-0-9.eE]+)\s+\[molec\./uc/Pa\]",text)
    if expected_seed is not None and data.get("seed") != expected_seed:
        raise ValueError("RASPA seed differs from requested independent seed")
    z=float(formula_units)*math.prod(unit_cells)
    if not math.isfinite(z) or z<=0: raise ValueError("positive formula units required")
    energy=data["properties"]["averageEnergies"]["Framework-Molecule"][gas]
    total=energy["total"]
    parts=[energy[k]["mean"] for k in ("vanDerWaals","coulombReal","coulombFourier","tailCorrection")]
    if not math.isclose(sum(parts),total["mean"],rel_tol=1e-6,abs_tol=1e-6):
        raise ValueError("host-guest energy decomposition does not sum to total")
    loading=parsed["absolute_loading_molecules_per_cell"]/z
    error=parsed["absolute_loading_error_molecules_per_cell"]/z
    u=total["mean"]*.008314462618/z
    ue=total["confidence"]*.008314462618/z
    if not all(math.isfinite(v) for v in (loading,error,u,ue)) or min(loading,error,ue)<0:
        raise ValueError("invalid final statistics")
    accepted=parsed["move_diagnostic"]["swap_accepted"]
    accessibility=bool(accepted and sum(accepted)>0 and loading>0)
    return dict(loading_molecules_fu=loading, loading_confidence_molecules_fu=error,
                host_guest_energy_kJ_mol_fu=u, host_guest_energy_confidence_kJ_mol_fu=ue,
                mean_energy_per_adsorbed_molecule_kJ_mol=u/loading if loading>0 else None,
                accessibility_passed=accessibility, seed=data.get("seed"),raspa_version=data.get("version"),
                sampling_observation="positive_loading_with_accepted_swaps" if accessibility else "insertion_or_loading_unresolved_not_proof_of_exclusion",
                energy_components_K=energy,raw=parsed,
                widom_henry_molecules_fu_bar=[float(x)*1e5/float(formula_units) for x in henry[-1]] if henry else None,
                normalization="total simulation box / (formula_units_per_CIF * supercell_product)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(parse_output(args.output), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
