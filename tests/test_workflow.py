import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from flexible_mof_gateopen.adsorption_fit import fit_isotherm
from parse_raspa_output import parse_point
from run_workflow import Workflow
from prepare_guest_rdkit import validate_geometry


def test_isotherm_and_extensive_energy_fit():
    f=np.geomspace(.001,100,15);n=3*2*f/(1+2*f);u=-25*n+2*n*n
    r=fit_isotherm(f,n,u)
    assert r["n_sat_molecules_fu"]==pytest.approx(3,rel=1e-5)
    assert r["b_per_bar"]==pytest.approx(2,rel=1e-5)
    assert r["interaction_energy_fit"]["a_kJ_mol_adsorbate"]==pytest.approx(-25)
    assert r["interaction_energy_fit"]["b"]==pytest.approx(2)
    assert fit_isotherm([1,2,3],[0,0,0])["status"]=="unidentifiable_zero_loading"


def test_native_energy_and_supercell_normalization(tmp_path):
    out=tmp_path/"output";out.mkdir()
    path=out/"output_268_1e5.s0.txt"
    path.write_text("Abs. loading average 24 +/- 2 [molecules/cell]\n"
                    "Abs. loading average 3 +/- .25 [molecules/uc]\n"
                    "Excess loading average 23 +/- 2 [molecules/cell]\n"
                    "Swap (CBMC) accepted: 12 5 0\n"
                    "Average Henry coefficient: 0.00004 +/- 0.000002 [molec./uc/Pa]\n")
    energy={k:{"mean":v,"confidence":1} for k,v in
            dict(vanDerWaals=-900,coulombReal=-30,coulombFourier=-70,tailCorrection=0,total=-1000).items()}
    path.with_suffix(".json").write_text(json.dumps(dict(seed=417,version="3.1.0",properties=dict(averageEnergies={"Framework-Molecule":{"CO2":energy}}))))
    r=parse_point(tmp_path,"CO2",4,[2,2,2],417)
    assert r["loading_molecules_fu"]==.75
    assert r["host_guest_energy_kJ_mol_fu"]==pytest.approx(-1000*.008314462618/32)
    assert r["widom_henry_molecules_fu_bar"][0]==pytest.approx(1)
    assert r["accessibility_passed"]
    with pytest.raises(ValueError,match="normalization"):
        parse_point(tmp_path,"CO2",4,[1,1,1],417)
    with pytest.raises(ValueError,match="seed"):
        parse_point(tmp_path,"CO2",4,[2,2,2],99)


def test_checkpoint_cache_and_tamper_detection(tmp_path):
    spec=tmp_path/"spec.json";spec.write_text("{}")
    w=Workflow(spec,tmp_path/"run")
    def step(d):
        (d/"file.txt").write_text("actual result");return {"answer":42}
    assert w.stage("test",step)=={"answer":42}
    assert w.stage("test",lambda _:pytest.fail("cache should avoid re-execution"))=={"answer":42}
    (tmp_path/"run/test/attempt_001/file.txt").write_text("changed")
    with pytest.raises(ValueError,match="checkpoint"):
        w.stage("test",step)
    spec.write_text('{"temperature_K":268}')
    with pytest.raises(ValueError,match="inputs/code changed"):
        Workflow(spec,tmp_path/"run")


def test_co2_mmff_converges_but_fails_geometry_gate():
    Chem=pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    mol=Chem.MolFromSmiles("O=C=O")
    AllChem.EmbedMolecule(mol,randomSeed=417)
    assert AllChem.MMFFOptimizeMolecule(mol)==0
    with pytest.raises(ValueError,match="bond|CO2"):
        validate_geometry(mol)


def test_analysis_will_not_invent_host_free_energy(tmp_path):
    w=Workflow.__new__(Workflow)
    w.spec={"pressure_bar":[.001,.01,.1],"temperature_K":268,"henry_review":{"closed":"test","open":"test"}}
    row=lambda p:dict(pressure_bar=p,fugacity_bar=p,loading_molecules_fu=p,
                     host_guest_energy_kJ_mol_fu=-20*p,loading_confidence_molecules_fu=p*.1,accessibility_passed=True)
    ads={name:dict(rows=[row(p) for p in w.spec["pressure_bar"]],volume_A3_fu=200) for name in ("closed","open")}
    result=w.analyze(None,ads,tmp_path)
    assert result["status"]=="blocked"
    assert "host free energy unavailable" in result["blocking_reasons"]
    assert "equilibrium" not in result
