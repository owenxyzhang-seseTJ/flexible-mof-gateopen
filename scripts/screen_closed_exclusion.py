#!/usr/bin/env python3
"""Diagnostic zero/residual closed-adsorption sensitivity from verified run data.

This command deliberately does not certify short open-isotherm data or proxy
host free energies. Its outputs are conditional scenarios, not gate pressures.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen.advanced_thermo import state_adsorption_omega,R,PV
    from flexible_mof_gateopen.reservoir import reservoir_point
    from flexible_mof_gateopen.adsorption_fit import fit_isotherm
    from flexible_mof_gateopen.preflight import sha256
except ModuleNotFoundError:
    from advanced_thermo import state_adsorption_omega,R,PV
    from reservoir import reservoir_point
    from adsorption_fit import fit_isotherm
    from preflight import sha256


def run(spec,base):
    paths={k:(base/spec[k]).resolve() for k in ("closed_checkpoint","open_checkpoint","host_json","widom_json")}
    raw={k:json.loads(v.read_text()) for k,v in paths.items()}
    ads={}
    for name in ("closed","open"):
        checkpoint=raw[name+"_checkpoint"]
        if any(not Path(p).is_file() or sha256(p)!=h for p,h in checkpoint["artifacts"].items()):
            raise ValueError(f"{name} checkpoint artifact modified/missing")
        ads[name]=checkpoint["result"]
    host=raw["host_json"];widom=raw["widom_json"]
    t=float(spec["temperature_K"])
    if host["temperature_K"]!=t or widom["temperature_K"]!=t or widom["gas"]!=spec["fluid"]:
        raise ValueError("temperature/guest mismatch")
    if not spec.get("open_isotherm_assumption") or host["method"]!="static_energy_plus_elastic_network_entropy_proxy":
        raise ValueError("explicit open isotherm assumption and supported host proxy output required")
    if any(row["loading_molecules_fu"]!=0 for row in ads["closed"]["rows"]):
        raise ValueError("closed zero reference contradicts positive observed loading")
    # Bind Widom to the closed model actually used for the adsorption checkpoint.
    if any(raw["closed_checkpoint"]["artifacts"].get(p)!=h or sha256(p)!=h for p,h in widom["source_inputs_sha256"].items()):
        raise ValueError("Widom input does not match closed adsorption model")
    states=[]
    for name in ("closed","open"):
        data=ads[name];pressures=sorted(set(r["pressure_bar"] for r in data["rows"]))
        f=[];n=[];errors=[]
        for p in pressures:
            rows=[r for r in data["rows"] if r["pressure_bar"]==p]
            eos=reservoir_point(t,p,backend=spec["reservoir"],fluid=spec["fluid"])
            if any(not np.isclose(r["fugacity_bar"],eos["fugacity_bar"],rtol=1e-8) for r in rows):
                raise ValueError("simulation fugacity differs from requested reservoir")
            vals=[r["loading_molecules_fu"] for r in rows]
            f.append(rows[0]["fugacity_bar"]);n.append(float(np.mean(vals)))
            errors.append(max([r["loading_confidence_molecules_fu"] for r in rows]+[float(np.ptp(vals))/2]))
        record=next(s for s in host["states"] if s["name"]==name)
        if not np.isclose(record["volume_A3_fu"],data["volume_A3_fu"],rtol=1e-8):
            raise ValueError("host/adsorption volume normalization mismatch")
        states.append(dict(name=name,pressure_bar=pressures,fugacity_bar=f,loading_molecules_fu=n,
                           loading_native_or_seed_sensitivity=errors,volume_A3_fu=data["volume_A3_fu"],
                           henry_confirmed=True,
                           henry_assumption=spec["open_isotherm_assumption"] if name=="open" else "analytical zero scenario"))
    lo=max(s["pressure_bar"][0] for s in states);hi=min(s["pressure_bar"][-1] for s in states)
    maxf=reservoir_point(t,hi,backend=spec["reservoir"],fluid=spec["fluid"])["fugacity_bar"]
    states[0]["adsorption_model"]=dict(kind="nonadsorbing",interpretation="conditional_scenario",
        evidence=f"{len(widom['rows'])} independent empty-host Widom records in {paths['widom_json']}; zero branch explicitly assumed, not proven",
        max_fugacity_bar=maxf)
    scale=float(widom["sensitivity_K_scale_molecules_fu_bar"])
    baseline=float(host["parameter_scenarios"][0]["delta_F_host_proxy_kJ_mol_fu"])
    delta_v=states[1]["volume_A3_fu"]-states[0]["volume_A3_fu"]
    def evaluate(k,df=baseline,opened=None,langmuir=False):
        opened=opened or states[1]
        closed=deepcopy(states[0])
        if k:
            closed["adsorption_model"].update(kind="henry_scenario",K_molecules_fu_bar=k,
                finite_pressure_assumption="Assume N_closed=K*f throughout sampled vapor range; not a bound")
        fit=fit_isotherm(opened["fugacity_bar"],opened["loading_molecules_fu"]) if langmuir else None
        def delta(p):
            fug=reservoir_point(t,p,backend=spec["reservoir"],fluid=spec["fluid"])["fugacity_bar"]
            a=(-R*t*fit["n_sat_molecules_fu"]*np.log1p(fit["b_per_bar"]*fug) if fit else
               state_adsorption_omega(opened,[fug],t)[0])
            return float(df+PV*p*delta_v+a-state_adsorption_omega(closed,[fug],t)[0])
        grid=np.geomspace(lo,hi,401);y=[delta(p) for p in grid];crossings=[]
        for a,b,ya,yb in zip(grid[:-1],grid[1:],y[:-1],y[1:]):
            if ya*yb<0:
                crossings.append(dict(pressure_bar=float(brentq(delta,a,b,xtol=1e-11)),
                                      direction="closed_to_open" if ya>0 else "open_to_closed"))
        return dict(K_molecules_fu_bar=k,delta_F_host_kJ_mol_fu=df,crossings=crossings,
                    maximum_closed_loading_scenario_molecules_fu=k*maxf,
                    maximum_closed_stabilization_kJ_mol_fu=R*t*k*maxf,
                    delta_Omega_at_max_pressure_kJ_mol_fu=y[-1])
    residual=[dict(scale_multiplier=m,**evaluate(m*scale)) for m in (0,1,10,100)]
    host_sensitivity=[dict(name=s["name"],**evaluate(0,s["delta_F_host_proxy_kJ_mol_fu"])) for s in host["parameter_scenarios"]]
    loading=[]
    for sign in (-1,1):
        opened=deepcopy(states[1]);opened["loading_molecules_fu"]=np.maximum(0,np.array(opened["loading_molecules_fu"])+sign*np.array(opened["loading_native_or_seed_sensitivity"])).tolist()
        loading.append(dict(perturbation=sign,**evaluate(0,opened=opened)))
    return dict(classification="diagnostic_sensitivity_only_not_a_validated_gate_pressure",
        temperature_K=t,pressure_range_bar=[lo,hi],max_fugacity_bar=maxf,
        source_files_sha256={str(p):sha256(p) for p in paths.values()},states=states,
        widom_seed_K_molecules_fu_bar=[r["widom_henry_molecules_fu_bar"] for r in widom["rows"]],
        K_scale_molecules_fu_bar=scale,residual_adsorption_scenarios=residual,
        host_proxy_parameter_scenarios=host_sensitivity,
        open_loading_perturbation_scenarios=loading,
        open_langmuir_integration_model_scenario=evaluate(0,langmuir=True),
        experimental_gateopen_pressure_bar=None,
        limitations=[spec["open_isotherm_assumption"],
            "Henry evidence does not exclude cooperative adsorption or finite-sampling rare-site bias",
            "Zero occupancy during Widom is enforced by the algorithm, not a measurement of zero adsorption",
            "Residual, host and loading scenarios are NOT confidence intervals",
            "The open pressure grid is sparse; a Langmuir alternative is a model sensitivity, not a validated replacement",
            "Host entropy is uncalibrated and uses a different potential from static energy; ML phonon instability remains unresolved",
            "P_eq is conditional equilibrium, not experimental adsorption/desorption opening/closing pressure"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec",required=True,type=Path);parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists():raise ValueError("use a fresh output path")
    result=run(json.loads(args.spec.read_text()),args.spec.resolve().parent)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result["residual_adsorption_scenarios"],indent=2))
