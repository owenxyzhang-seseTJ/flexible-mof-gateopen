from __future__ import annotations

import pytest


def test_cif_audit_charge_neutral(tmp_path) -> None:
    pytest.importorskip("gemmi")
    from flexible_mof_gateopen.cif_audit import audit_cif

    path = tmp_path / "charged.cif"
    path.write_text(
        """data_test
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_charge
A1 C 0 0 0 1 0.25
A2 O 0.5 0.5 0.5 1 -0.25
""",
        encoding="utf-8",
    )
    result = audit_cif(path)
    assert result["atom_count"] == 2
    assert result["charge_neutral"] is True
    assert abs(result["charge_sum_e"]) < 1e-12
