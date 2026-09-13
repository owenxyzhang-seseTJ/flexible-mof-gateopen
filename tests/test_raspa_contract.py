from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_raspa_template_requires_cif_charges(tmp_path) -> None:
    module = _load("prepare_raspa_inputs", ROOT / "scripts" / "prepare_raspa_inputs.py")
    cif = tmp_path / "host.cif"
    force_field = tmp_path / "force_field.json"
    cif.write_text("data_host\n", encoding="utf-8")
    force_field.write_text("{}\n", encoding="utf-8")
    output = tmp_path / "raspa"
    # Exercise the same writer used by the CLI without launching RASPA.
    output.mkdir()
    point = output / "P_1bar"
    point.mkdir()
    (point / "host.cif").write_text(cif.read_text(encoding="utf-8"), encoding="utf-8")
    (point / "force_field.json").write_text(force_field.read_text(encoding="utf-8"), encoding="utf-8")
    (point / "simulation.json").write_text(json.dumps(module.simulation(268.0, 1.0, "host")), encoding="utf-8")
    sim = json.loads((point / "simulation.json").read_text(encoding="utf-8"))
    system = sim["Systems"][0]
    assert system["ChargeMethod"] == "Ewald"
    assert system["UseChargesFrom"] == "CIF_File"
    assert sim["Components"][0]["FugacityCoefficient"] == 1.0


def test_raspa_parser_exposes_zero_move_flag(tmp_path) -> None:
    module = _load("parse_raspa_output", ROOT / "scripts" / "parse_raspa_output.py")
    output = tmp_path / "output.s0.txt"
    output.write_text(
        """Abs. loading average   0.000000e+00 +/-  0.000000e+00 [molecules/cell]
Excess loading average  -6.2e-01 +/-  0.000000e+00 [molecules/cell]
    Translation          total:                 0          0          0
    Swap (CBMC)          constructed:      131020          0          0
    Swap (CBMC)          accepted:              0          0          0
""",
        encoding="utf-8",
    )
    result = module.parse_output(output)
    assert result["absolute_loading_molecules_per_cell"] == 0.0
    assert result["move_diagnostic"]["translation_total"] == [0, 0, 0]
    assert result["move_diagnostic"]["translation_constructed"] is None
    assert result["move_diagnostic"]["swap_constructed"] == [131020, 0, 0]
    assert result["move_diagnostic"]["swap_accepted"] == [0, 0, 0]
