import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import os

import pytest

from flexible_mof_gateopen.preflight import pressure_gate, structure_gate, sha256
from flexible_mof_gateopen.free_energy import cumulative_error, cumulative_integral

ROOT = Path(__file__).resolve().parents[1]


def voc():
    return {"name": "test_vapor", "kind": "VOC", "source": "synthetic test", "saturation": {
        "method": "value", "temperature_K": 298.15, "pressure_bar": .1, "source": "synthetic test"}}


def test_saturation_bound_and_temperature():
    assert pressure_gate(voc(), 298.15, [.01, .1])["relative_pressure"] == pytest.approx([.1, 1])
    for temperature, pressures in [(298.15, [.10001]), (299, [.01]), (298.15, [float("nan")]), (298.15, [0])]:
        with pytest.raises(ValueError):
            pressure_gate(voc(), temperature, pressures)
    with pytest.raises(ValueError):
        pressure_gate({"name": "ethanol", "kind": "VOC"}, 298, [.01])


def test_antoine_units_and_validity():
    profile = voc()
    profile["saturation"] = {"method": "antoine", "A": 5., "B": 1000., "C": 0.,
        "temperature_unit": "K", "pressure_unit": "kPa", "valid_temperature_K": [290, 310], "source": "synthetic"}
    expected = 10**(5-1000/300) * .01
    assert pressure_gate(profile, 300, [.01])["psat_bar"] == pytest.approx(expected)
    with pytest.raises(ValueError):
        pressure_gate(profile, 320, [.01])


def test_shared_node_uncertainty_and_bounds():
    # Integral = (1 + ln(2)/2)*N0 + ln(2)*N1 + ln(2)/2*N2.
    w = [1 + math.log(2)/2, math.log(2), math.log(2)/2]
    assert cumulative_error([1, 2, 4], [1, 2, 3]) == pytest.approx(math.sqrt(sum((a*b)**2 for a,b in zip(w,[1,2,3]))))
    with pytest.raises(ValueError):
        cumulative_integral([1, 2], [1, 2], target_pressure_bar=3)


@pytest.fixture
def reviewed(tmp_path):
    pytest.importorskip("pymatgen.io.cif")
    cif = tmp_path / "water.cif"
    cif.write_text("""data_water
_cell_length_a 30
_cell_length_b 30
_cell_length_c 30
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_charge
O1 O .5 .5 .5 1 -0.8
H1 H .532 .5 .5 1 0.4
H2 H .492 .531 .5 1 0.4
""")
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps({"fixture": "synthetic water, not a real MOF cleaning certificate"}))
    data = {"schema_version": 1, "status": "passed", "source_sha256": sha256(cif), "candidate_sha256": sha256(cif),
        "transformation_record": "Synthetic ordered water test fixture, no transformation",
        "source_audit": {"path": "audit.json", "sha256": sha256(audit)},
        "candidate_audit": {"path": "audit.json", "sha256": sha256(audit)},
        "disorder": {"status": "absent", "reason": "synthetic ordered fixture"},
        "solvent": {"status": "retained_with_reason", "reason": "water fixture tests the gate, not MOF activation"},
        "chemistry": {k: {"status": "passed", "evidence": "synthetic test geometry"} for k in
            ("connectivity", "hydrogen_hosts", "metal_coordination", "nonbonded_contacts", "composition_charge_balance")}}
    report = tmp_path / "cleaning.json"
    report.write_text(json.dumps(data))
    return cif, report, data


def test_structure_gate_rejects_unresolved_or_changed_model(reviewed):
    cif, report, data = reviewed
    assert structure_gate(cif, report, charges=True)["atom_count"] == 3
    data["solvent"]["status"] = "unknown"
    report.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="solvent"):
        structure_gate(cif, report)
    data["solvent"]["status"] = "removed"
    report.write_text(json.dumps(data))
    cif.write_text(cif.read_text().replace("H1 H .532", "H1 H .533"))
    with pytest.raises(ValueError, match="hash"):
        structure_gate(cif, report)
    data["candidate_sha256"] = sha256(cif)
    report.write_text(json.dumps(data))
    cif.write_text(cif.read_text().replace(".5 1 -0.8", ".5 0.5 -0.8"))
    data["candidate_sha256"] = sha256(cif)
    report.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="occupancies"):
        structure_gate(cif, report)


def test_host_cli_blocks_unreviewed_cif_before_ml_load(reviewed,tmp_path):
    pytest.importorskip("scipy")
    cif,report,data=reviewed
    data["solvent"]["status"]="unknown"
    report.write_text(json.dumps(data))
    spec=tmp_path/"host.json"
    spec.write_text(json.dumps({"states":[{"name":"test","cif":str(cif),"cleaning_report":str(report),"formula_units":1}]}))
    out=tmp_path/"blocked.json"
    run=subprocess.run([sys.executable,str(ROOT/"scripts/advanced_free_energy.py"),"host-energy",
                        "--spec",str(spec),"--out",str(out)],capture_output=True,text=True)
    assert run.returncode!=0 and "unresolved solvent" in run.stderr
    assert not out.exists()


def test_cli_voc_limit_and_valid_custom_guest(reviewed, tmp_path):
    cif, report, _ = reviewed
    profile = tmp_path / "gas.json"
    profile.write_text(json.dumps(voc()))
    molecule = tmp_path / "molecule.json"
    molecule.write_text(json.dumps({"Type": "rigid", "pseudoAtoms": [["O_co2", [0,0,1.149]], ["C_co2", [0,0,0]], ["O_co2", [0,0,-1.149]]]}))
    out = tmp_path / "sim"
    cmd = [sys.executable, str(ROOT / "scripts/prepare_raspa_inputs.py"), "--cif", str(cif),
        "--cleaning-report", str(report), "--gas-profile", str(profile), "--molecule-json", str(molecule),
        "--force-field", str(ROOT / "examples/xpcu5zn_co2_268K/data/force_field.json"),
        "--temperature", "298.15", "--output-dir", str(out), "--unit-cells", "1", "1", "1", "--cutoff", "12", "--pressures"]
    failed = subprocess.run(cmd+["1"], capture_output=True, text=True)
    assert failed.returncode != 0 and "Psat" in failed.stderr and not out.exists()
    success = subprocess.run(cmd+["0.1"], capture_output=True, text=True)
    assert success.returncode == 0, success.stderr
    sim = json.loads((out / "P_0p1bar/simulation.json").read_text())
    assert sim["Components"][0]["Name"] == "test_vapor"
    assert (out / "P_0p1bar/test_vapor.json").exists()
    assert subprocess.run(cmd+["0.1"], capture_output=True).returncode != 0


def test_rdkit_ethanol(tmp_path):
    pytest.importorskip("rdkit")
    spec = importlib.util.spec_from_file_location("guest", ROOT / "scripts/prepare_guest_rdkit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.prepare("CCO", tmp_path / "ethanol", conformers=3)
    from rdkit import Chem
    mol = next(iter(Chem.SDMolSupplier(str(tmp_path / "ethanol/guest.sdf"), removeHs=False)))
    assert mol.GetNumAtoms() == 9
    assert result["force_field"] == "MMFF94s"
    assert any(c["id"] == result["selected_conformer"] and c["status"] == 0 for c in result["conformers"])


@pytest.mark.skipif(not os.environ.get("RASPA_SMOKE_BIN"), reason="optional native RASPA integration")
@pytest.mark.parametrize("reservoir", ["ideal", "coolprop"])
def test_native_raspa_preflight(reviewed, tmp_path, reservoir):
    if reservoir == "coolprop": pytest.importorskip("CoolProp")
    cif, report, _ = reviewed
    out = tmp_path / "native"
    command = [sys.executable, str(ROOT / "scripts/prepare_raspa_inputs.py"), "--cif", str(cif),
        "--cleaning-report", str(report), "--force-field", str(ROOT / "examples/xpcu5zn_co2_268K/data/force_field.json"),
        "--temperature", "298.15", "--output-dir", str(out), "--unit-cells", "1", "1", "1", "--cutoff", "12",
        "--pressures", "0.01", "--raspa-bin", os.environ["RASPA_SMOKE_BIN"], "--run", "--reservoir", reservoir]
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    outputs = list((out / "P_0p01bar/output").glob("output_*.txt"))
    assert outputs
    text = outputs[0].read_text()
    assert "Abs. loading average" in text and "Ewald alpha" in text
