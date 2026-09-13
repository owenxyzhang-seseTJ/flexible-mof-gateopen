from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from flexible_mof_gateopen.free_energy import integrate_required_host_penalty


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "xpcu5zn_co2_268K"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_public_example_has_no_machine_paths() -> None:
    for path in EXAMPLE.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".csv", ".json", ".md", ".cif"}:
            text = path.read_text(encoding="utf-8", errors="replace")
            assert "/Users/" not in text
            assert "Mobile Documents" not in text


def test_provenance_hashes_and_charges() -> None:
    provenance = json.loads((EXAMPLE / "provenance.json").read_text(encoding="utf-8"))
    for name, expected in provenance["derived_sha256"].items():
        assert _sha256(EXAMPLE / "structures" / name) == expected
    pytest.importorskip("gemmi")
    from flexible_mof_gateopen.cif_audit import audit_cif

    for name in ("alpha_pacmof_charged.cif", "gamma_pacmof_charged.cif"):
        result = audit_cif(EXAMPLE / "structures" / name)
        assert result["charge_tag_present"] is True
        assert result["charge_neutral"] is True
        assert result["max_occupancy_deviation"] == 0.0


def test_example_reproduces_19_bar_screen() -> None:
    rows = list(csv.DictReader((EXAMPLE / "data" / "fixed_host_isotherms_268K.csv").open(encoding="utf-8")))
    phases = {phase: sorted((row for row in rows if row["phase"] == phase), key=lambda row: float(row["pressure_bar"])) for phase in ("open_alpha", "closed_gamma")}
    pressures = [float(row["pressure_bar"]) for row in phases["open_alpha"]]
    assert pressures == [float(row["pressure_bar"]) for row in phases["closed_gamma"]]
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
    assert abs(value - 67.94) < 0.05
    assert error is not None and abs(error - 0.94) < 0.05


def test_archived_validation_has_no_failures() -> None:
    report = json.loads((EXAMPLE / "reports" / "validation_report_268K.json").read_text(encoding="utf-8"))
    assert report["status_counts"]["fail"] == 0
    assert report["analysis"]["required_host_at_19_bar_kJ_per_formula_unit"] > 60
