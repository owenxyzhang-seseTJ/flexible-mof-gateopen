import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import pytest
np=pytest.importorskip("numpy")
pytest.importorskip("scipy")
from flexible_mof_gateopen.advanced_thermo import (
    adsorption_omega,all_crossings,phase_diagram,qha_minimize,landscape_extrema,landscape_scan,PV,R)
from flexible_mof_gateopen.reservoir import reservoir_point
from flexible_mof_gateopen.framework_ti import mc_ti,KB
from flexible_mof_gateopen.host_models import energy_difference
from flexible_mof_gateopen.host_models import harmonic_state_differences

ROOT=Path(__file__).resolve().parents[1]


def test_exact_piecewise_integral_and_fugacity_remapping():
    f=[1,math.e,math.e**2]; n=[1,2,3]
    # first Henry area=1, then integral of 1+ln f
    result=adsorption_omega(f,n,[.1,math.exp(.5),math.e**2],300,henry_confirmed=True)
    assert result == pytest.approx(-R*300*np.array([.1,1+.5+.5*.5**2,5]))
    with pytest.raises(ValueError,match="Henry"):
        adsorption_omega(f,n,[1],300)
    with pytest.raises(ValueError,match="extrapolation"):
        adsorption_omega(f,n,[100],300,henry_confirmed=True)


def test_all_roots_tangent_endpoint_and_plateau():
    assert all_crossings([1,2,4,8],[1,-1,1,-1])["pressure_bar"] == pytest.approx([2**.2928932188,2**1.5,2**2.7071067812])
    assert all_crossings([1,2,3],[1,0,1])["pressure_bar"] == [2]
    assert all_crossings([1,2,3],[0,0,1])["coexistence_intervals_bar"] == [[1,2]]


def test_multiphase_mechanics_and_metastable_crossings():
    p=np.geomspace(.1,100,100)
    states=[]
    for name,f0,v in [("A",0,100),("B",.1,10),("C",100,50)]:
        states.append(dict(name=name,host_F_kJ_mol_fu=f0,host_source="analytic test",volume_A3_fu=v,
                           fugacity_bar=[.01,101],loading_molecules_fu=[0,0],henry_confirmed=True))
    result=phase_diagram(dict(temperature_K=300,pressure_bar=p,states=states),
                         [reservoir_point(300,x) for x in p])
    root=result["crossings"][0]
    assert root["pressure_bar"] == pytest.approx(.1/(PV*90),rel=1e-4)
    assert root["globally_stable"] and root["stable_before"]=="A" and root["stable_after"]=="B"
    assert "C" not in result["stable_state"]


def test_qha_no_cross_basin_extrapolation():
    v=np.linspace(90,110,21)
    r=qha_minimize(v,(v-100)**2)
    assert r["volume_A3_fu"]==pytest.approx(100,abs=.001)
    assert not r["boundary_minimum"]
    assert qha_minimize(v,np.arange(len(v)))["boundary_minimum"]
    with pytest.raises(ValueError): qha_minimize([1,2,3],[1,2,3])


def test_double_well_barrier_and_spinodal():
    q=np.linspace(-2,2,201)
    r=landscape_extrema(q,(q*q-1)**2)
    assert len(r["barriers"])==2
    assert [b["barrier_kJ_mol_fu"] for b in r["barriers"]]==pytest.approx([1,1],abs=1e-6)
    h=np.linspace(0,2,41)
    p=np.exp(h)
    r=landscape_scan(q,p,[(q*q-1)**2-x*q for x in h])
    critical=8/(3*math.sqrt(3))
    assert len(r["spline_spinodals"])==1
    sp=r["spline_spinodals"][0]
    assert math.log(sp["pressure_bar"])==pytest.approx(critical,rel=.002)
    assert sp["q"]==pytest.approx(-1/math.sqrt(3),rel=.002)


def test_mc_ti_identical_reference_is_exact():
    k=2.;temp=300.;bound=.3;d=3
    r=mc_ti(lambda x:k*(x@x)/2,d,temperature_K=temp,spring_eV_A2=k,bound_A=bound,
            lambdas=[0,.5,1],warmup=10,samples=100,stride=1)
    assert r["delta_F_TI_eV"]==0
    integral=math.sqrt(2*math.pi*KB*temp/k)*math.erf(bound*math.sqrt(k/(2*KB*temp)))
    assert r["F_conf_eV"]==pytest.approx(-d*KB*temp*math.log(integral))


def test_mc_ti_nontrivial_harmonic_target():
    # Broad bounds approximate full Gaussian; deterministic seeded statistical validation.
    r=mc_ti(lambda x:2*(x@x)/2,3,spring_eV_A2=1,bound_A=2,step_A=.3,
            lambdas=np.linspace(0,1,21),warmup=500,samples=2500,stride=3,seed=32)
    exact=1.5*KB*300*math.log(2)
    assert r["delta_F_TI_eV"]==pytest.approx(exact,abs=.0025)


def test_energy_composition_gate():
    with pytest.raises(ValueError,match="composition"):
        energy_difference([dict(composition_per_fu={"C":1},energy_kJ_mol_fu=1),
                           dict(composition_per_fu={"C":2},energy_kJ_mol_fu=2)])


def test_static_plus_vibrational_free_energy_not_entropy_only():
    def state(name,e,f,s):
        return dict(name=name,base=dict(composition_per_fu={"C":1},energy_kJ_mol_fu=e),
                    min_frequency_THz=0,thermal=[dict(temperature_K=300,F_vib_kJ_mol_fu=f,
                    S_J_mol_fu_K=s,F_total_kJ_mol_fu=e+f)])
    phases=[state("closed",0,2,10),state("open",10,1,20)]
    result=harmonic_state_differences(phases)[0]["thermal"][0]
    assert result["delta_E_static_kJ_mol_fu"]==10
    assert result["minus_T_delta_S_kJ_mol_fu"]==-3
    assert result["delta_U_vib_kJ_mol_fu"]==2
    assert result["delta_F_host_harmonic_kJ_mol_fu"]==9  # not 10-3=7
    phases[1]["min_frequency_THz"]=-.01
    assert harmonic_state_differences(phases)[0]["status"]=="blocked_unstable_phonons"


def test_coolprop_vapor_and_saturation():
    pytest.importorskip("CoolProp")
    r=reservoir_point(268,19,backend="coolprop")
    assert .8<r["fugacity_coefficient"]<.95
    assert r["fugacity_bar"]==pytest.approx(19*r["fugacity_coefficient"])
    with pytest.raises(ValueError,match="Psat"):
        reservoir_point(268,100,backend="coolprop")
    saturated=reservoir_point(268,r["psat_bar"],backend="coolprop")
    assert saturated["fugacity_bar"]>0


def test_advanced_cli_and_preserves_output(tmp_path):
    spec=tmp_path/"r.json";out=tmp_path/"out.json"
    spec.write_text(json.dumps(dict(temperature_K=300,pressure_bar=[1,2])))
    cmd=[sys.executable,str(ROOT/"scripts/advanced_free_energy.py"),"reservoir","--spec",str(spec),"--out",str(out)]
    r=subprocess.run(cmd,capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    assert json.loads(out.read_text())["points"][1]["fugacity_bar"]==2
    assert subprocess.run(cmd,capture_output=True).returncode!=0


@pytest.mark.parametrize("command,filename",[("phase-diagram","phase_diagram.json"),("qha","qha.json")])
def test_bundled_advanced_cli(tmp_path,command,filename):
    out=tmp_path/"out.json"
    cmd=[sys.executable,str(ROOT/"scripts/advanced_free_energy.py"),command,"--spec",
         str(ROOT/"examples/advanced_validation"/filename),"--out",str(out)]
    r=subprocess.run(cmd,capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    result=json.loads(out.read_text())
    assert result["slices"]


def test_landscape_cli_analytic_and_host_gate(tmp_path):
    q=np.linspace(-2,2,101);tilts=np.linspace(0,2,11)
    spec=tmp_path/"landscape.json";out=tmp_path/"out.json"
    spec.write_text(json.dumps(dict(analytic_test=True,landscape_source="analytic tilted double well",
        coordinate_measure="dimensionless q",q=q.tolist(),pressure_bar=np.exp(tilts).tolist(),
        omega_kJ_mol_fu=[((q*q-1)**2-h*q).tolist() for h in tilts])))
    cmd=[sys.executable,str(ROOT/"scripts/advanced_free_energy.py"),"landscape","--spec",str(spec),"--out",str(out)]
    r=subprocess.run(cmd,capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    assert len(json.loads(out.read_text())["spline_spinodals"])==1
    spec.write_text(json.dumps(dict(states=[])))
    failed=subprocess.run([sys.executable,str(ROOT/"scripts/advanced_free_energy.py"),"host-energy",
                           "--spec",str(spec),"--out",str(tmp_path/"blocked.json")],capture_output=True,text=True)
    assert failed.returncode!=0 and "states required" in failed.stderr


@pytest.mark.skipif(not os.environ.get("RUN_ML_SMOKE"),reason="optional pretrained checkpoint inference")
@pytest.mark.parametrize("backend",["chgnet","mace"])
def test_actual_pretrained_energy(backend):
    from ase.build import bulk
    from flexible_mof_gateopen.host_models import calculator,energy_record
    calc,meta=calculator(backend)
    a=bulk("Si","diamond",a=5.43);a.calc=calc
    r=energy_record(a,2)
    assert r["energy_eV_cell"]<0 and r["max_force_eV_A"]<.05
    assert len(meta["weights_sha256"])==64


@pytest.mark.skipif(not os.environ.get("RUN_ML_SMOKE"),reason="optional finite-displacement ML phonons")
def test_actual_chgnet_phonopy():
    from ase.build import bulk
    from flexible_mof_gateopen.host_models import calculator,phonon_free_energy
    a=bulk("Si","diamond",a=5.43);a.calc=calculator()[0]
    r=phonon_free_energy(a,2,[0,300],supercell=[2,2,2],mesh=[4,4,4])
    assert r["status"]=="harmonic_screen"
    assert r["thermal"][0]["F_vib_kJ_mol_fu"]>0
    assert r["thermal"][1]["F_vib_kJ_mol_fu"]<r["thermal"][0]["F_vib_kJ_mol_fu"]
