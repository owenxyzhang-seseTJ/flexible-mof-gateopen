#!/usr/bin/env python3
"""Portable analytic regressions; no material data, model loading or simulation."""
import json
import math
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen.advanced_thermo import adsorption_omega,state_adsorption_omega,closed_adsorption_sensitivity,R
except ModuleNotFoundError:
    from advanced_thermo import adsorption_omega,state_adsorption_omega,closed_adsorption_sensitivity,R


def main():
    result=adsorption_omega([1,math.e,math.e**2],[1,2,3],[.1,math.e**2],300,henry_confirmed=True)
    np.testing.assert_allclose(result,-R*300*np.array([.1,5]))
    model=dict(kind="nonadsorbing",interpretation="conditional_scenario",evidence="analytic regression",
        max_fugacity_bar=100,finite_pressure_assumption="linear synthetic adsorption")
    np.testing.assert_equal(state_adsorption_omega(dict(adsorption_model=model),[1,100],300),[0,0])
    closed=dict(name="closed",host_F_kJ_mol_fu=0,volume_A3_fu=100,host_source="analytic regression",adsorption_model=model)
    opened=dict(name="open",host_F_kJ_mol_fu=3,volume_A3_fu=100,host_source="analytic regression",
        adsorption_model=dict(model,kind="henry_scenario",K_molecules_fu_bar=.1))
    p=np.geomspace(.001,100,501)
    r=closed_adsorption_sensitivity(dict(temperature_K=300,pressure_bar=p,states=[closed,opened]),
        [dict(fugacity_bar=x) for x in p],[.01,.02])
    roots=[s["equilibrium_opening_pressures_bar"][0] for s in r["scenarios"]]
    np.testing.assert_allclose(roots,[3/(R*300*(.1-k)) for k in [0,.01,.02]],rtol=1e-5)
    try:state_adsorption_omega(dict(adsorption_model=model),[101],300)
    except ValueError:pass
    else:raise AssertionError("fugacity extrapolation was not rejected")
    print(json.dumps(dict(passed=True,checks=["Henry tail and piecewise integral","zero branch","residual shifts opening higher","extrapolation rejection"],
        scope="software only; not material accuracy"),indent=2))


if __name__=="__main__":main()
