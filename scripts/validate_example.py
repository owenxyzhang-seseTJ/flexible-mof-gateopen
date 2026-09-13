#!/usr/bin/env python3
"""Validate the bundled X-pcu-5-Zn example without running RASPA."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from flexible_mof_gateopen.cif_audit import audit_cif
    from flexible_mof_gateopen.free_energy import integrate_required_host_penalty
except ModuleNotFoundError:  # standalone copy inside a skill directory
    from cif_audit import audit_cif  # type: ignore[no-redef]
    from free_energy import integrate_required_host_penalty  # type: ignore[no-redef]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example", type=Path, default=Path("examples/xpcu5zn_co2_268K"))
    args = parser.parse_args()
    root = args.example
    provenance = json.loads((root / "provenance.json").read_text(encoding="utf-8"))
    hashes = {name: sha256(root / "structures" / name) for name in provenance["derived_sha256"]}
    if hashes != provenance["derived_sha256"]:
        raise SystemExit(f"derived CIF hash mismatch: {hashes}")
    audits = {
        name: audit_cif(root / "structures" / name)
        for name in ("alpha_pacmof_charged.cif", "gamma_pacmof_charged.cif")
    }
    if not all(item["charge_neutral"] for item in audits.values()):
        raise SystemExit(f"charge neutrality audit failed: {audits}")
    rows = list(csv.DictReader((root / "data" / "fixed_host_isotherms_268K.csv").open(encoding="utf-8")))
    phases = {phase: sorted((row for row in rows if row["phase"] == phase), key=lambda row: float(row["pressure_bar"])) for phase in ("open_alpha", "closed_gamma")}
    pressures = [float(row["pressure_bar"]) for row in phases["open_alpha"]]
    value, error = integrate_required_host_penalty(
        pressures,
        [float(row["absolute_loading_molecules_per_cell"]) for row in phases["open_alpha"]],
        [float(row["absolute_loading_molecules_per_cell"]) for row in phases["closed_gamma"]],
        temperature_K=268.0,
        open_formula_units=2.0,
        closed_formula_units=4.0,
        open_error=[float(row["absolute_loading_error_molecules_per_cell"]) for row in phases["open_alpha"]],
        closed_error=[float(row["absolute_loading_error_molecules_per_cell"]) for row in phases["closed_gamma"]],
        target_pressure_bar=19.0,
    )
    report = json.loads((root / "reports" / "validation_report_268K.json").read_text(encoding="utf-8"))
    result = {
        "ok": report["status_counts"]["fail"] == 0 and abs(value - 67.94) < 0.05 and all(item["charge_neutral"] for item in audits.values()),
        "charge_audits": audits,
        "derived_hashes": hashes,
        "computed_required_host_penalty_at_19_bar_kJ_per_formula_unit": value,
        "computed_error_kJ_per_formula_unit": error,
        "archived_status_counts": report["status_counts"],
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if not result["ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
