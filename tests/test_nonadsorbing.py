"""Analytical checks: exclusion is explicit, not inferred from failed swaps."""
import copy
import sys
from pathlib import Path
import numpy as np
import pytest

from flexible_mof_gateopen.advanced_thermo import state_adsorption_omega,closed_adsorption_sensitivity,R
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_workflow import Workflow


def model(kind="nonadsorbing"):
    return dict(kind=kind,interpretation="conditional_scenario",evidence="analytic test only",
                max_fugacity_bar=100,K_molecules_fu_bar=.02,
                finite_pressure_assumption="linear uptake over 0 to 100 bar for a test")


def test_exact_zero_and_linear_integrals_and_guards():
    state={"adsorption_model":model()}
    assert state_adsorption_omega(state,[.001,2,99],300)==pytest.approx([0,0,0])
    state["adsorption_model"]["kind"]="henry_scenario"
    assert state_adsorption_omega(state,[.001,2,99],300)==pytest.approx(-R*300*.02*np.array([.001,2,99]))
    with pytest.raises(ValueError,match="extrapolation"):
        state_adsorption_omega(state,[101],300)
    state["adsorption_model"].pop("finite_pressure_assumption")
    with pytest.raises(ValueError,match="finite_pressure_assumption"):
        state_adsorption_omega(state,[1],300)
    state["adsorption_model"]=model();state["adsorption_model"]["evidence"]=""
    with pytest.raises(ValueError,match="evidence"):
        state_adsorption_omega(state,[1],300)


def test_residual_adsorption_shifts_opening_higher_and_preserves_spec():
    p=np.geomspace(.001,100,501)
    closed=dict(name="closed",host_F_kJ_mol_fu=0,host_source="test",volume_A3_fu=100,adsorption_model=model())
    opened=dict(name="open",host_F_kJ_mol_fu=3,host_source="test",volume_A3_fu=100,
                adsorption_model=dict(model("henry_scenario"),K_molecules_fu_bar=.1))
    spec=dict(temperature_K=300,pressure_bar=p.tolist(),states=[closed,opened]);original=copy.deepcopy(spec)
    r=closed_adsorption_sensitivity(spec,[dict(fugacity_bar=x) for x in p],[.01,.02])
    roots=[s["equilibrium_opening_pressures_bar"][0] for s in r["scenarios"]]
    assert roots==pytest.approx([3/(R*300*(.1-k)) for k in [0,.01,.02]],rel=1e-5)
    assert roots[0]<roots[1]<roots[2]
    assert spec==original
    assert r["scenarios"][-1]["maximum_closed_adsorption_stabilization_kJ_mol_fu"]==pytest.approx(R*300*2)


def test_workflow_explicit_zero_not_automatic_pass(tmp_path):
    w=Workflow.__new__(Workflow);w.base=tmp_path
    w.spec=dict(pressure_bar=[.001,.01,.1],temperature_K=268,henry_review={"open":"synthetic Henry data"})
    ads={}
    for name in ("closed","open"):
        rows=[]
        for p in w.spec["pressure_bar"]:
            n=0 if name=="closed" else p
            rows.append(dict(pressure_bar=p,fugacity_bar=p,loading_molecules_fu=n,
                host_guest_energy_kJ_mol_fu=-20*n,loading_confidence_molecules_fu=0,accessibility_passed=n>0))
        ads[name]=dict(rows=rows,volume_A3_fu=100)
    host=dict(harmonic_state_differences=[dict(status="harmonic_estimate",thermal=[dict(delta_F_host_harmonic_kJ_mol_fu=.1)])])
    assert w.analyze(host,ads,tmp_path)["status"]=="blocked"
    w.spec["adsorption_models"]={"closed":dict(model(),sensitivity_K_molecules_fu_bar=[.001,.01])}
    r=w.analyze(host,ads,tmp_path)
    assert not r["blocking_reasons"]
    assert r["equilibrium_opening_pressures_bar"]
    assert len(r["closed_adsorption_sensitivity"]["scenarios"])==3
    # An explicit closed model must not suppress the OPEN-phase sampling gate.
    ads["open"]["rows"][0]["accessibility_passed"]=False
    assert w.analyze(host,ads,tmp_path)["status"]=="blocked"
    ads["open"]["rows"][0]["accessibility_passed"]=True
    ads["closed"]["rows"][0]["loading_molecules_fu"]=.01
    assert any("contradicts" in x for x in w.analyze(host,ads,tmp_path)["blocking_reasons"])
