#!/usr/bin/env python3
"""Resumable, fail-closed two-phase CIF-to-osmotic-pressure workflow.

Inputs are imported local/uploaded paths in a JSON specification. Each stage
retains immutable attempts, logs, hash-bound completed checkpoints and failures.
See references/end_to_end.md for schema and scientific boundaries.
"""
import argparse
import contextlib
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen.preflight import sha256,structure_gate,load_gas_profile,pressure_gate
    from flexible_mof_gateopen.host_models import harmonic_state_differences
    from flexible_mof_gateopen.advanced_thermo import phase_diagram,state_adsorption_omega,closed_adsorption_sensitivity
    from flexible_mof_gateopen.adsorption_fit import fit_isotherm
    from flexible_mof_gateopen.reservoir import reservoir_point
except ModuleNotFoundError:
    from preflight import sha256,structure_gate,load_gas_profile,pressure_gate
    from host_models import harmonic_state_differences
    from advanced_thermo import phase_diagram,state_adsorption_omega,closed_adsorption_sensitivity
    from adsorption_fit import fit_isotherm
    from reservoir import reservoir_point
from prepare_guest_rdkit import prepare,refine_tblite,validate_geometry
from prepare_raspa_inputs import CO2_JSON,pressure_tag
from parse_raspa_output import parse_point
from quick_host_free_energy import run as quick_host_run

SCRIPTS=Path(__file__).resolve().parent


def dump(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


class Workflow:
    def __init__(self,spec_path,out):
        self.spec_path=Path(spec_path).resolve(); self.base=self.spec_path.parent
        self.spec=json.loads(self.spec_path.read_text()); self.out=Path(out).resolve()
        self.out.mkdir(parents=True,exist_ok=True)
        self.results={}; self.failures={}
        # Paths, code and file bytes are part of the run identity, not only JSON.
        dependencies={str(self.spec_path):sha256(self.spec_path)}
        def discover(value):
            if isinstance(value,dict):
                for v in value.values(): discover(v)
            elif isinstance(value,list):
                for v in value: discover(v)
            elif isinstance(value,str):
                p=self.base/value
                if p.is_file(): dependencies[str(p.resolve())]=sha256(p)
        discover(self.spec)
        for p in SCRIPTS.glob("*.py"): dependencies[str(p)]=sha256(p)
        src=SCRIPTS.parent/"src/flexible_mof_gateopen"
        for p in src.glob("*.py"): dependencies[str(p)]=sha256(p)
        fingerprint=hashlib.sha256(json.dumps(dependencies,sort_keys=True).encode()).hexdigest()
        identity=self.out/"identity.json"
        if identity.exists() and json.loads(identity.read_text())["fingerprint"]!=fingerprint:
            raise ValueError("inputs/code changed; use a new run directory instead of stale resume")
        dump(identity,dict(fingerprint=fingerprint,dependencies=dependencies))

    def path(self,key): return (self.base/key).resolve()

    def command(self,args,directory):
        with (directory/"subprocess.log").open("a") as f:
            subprocess.run([str(v) for v in args],stdout=f,stderr=subprocess.STDOUT,check=True,
                           timeout=self.spec.get("subprocess_timeout_seconds",3600))

    def stage(self,name,fn):
        root=self.out/name; root.mkdir(exist_ok=True)
        checkpoint=root/"completed.json"
        if checkpoint.exists():
            record=json.loads(checkpoint.read_text())
            if any(not Path(p).is_file() or sha256(p)!=h for p,h in record["artifacts"].items()):
                raise ValueError(f"{name}: checkpoint artifact changed/missing")
            print(f"{name}: cached verified",flush=True)
            self.results[name]=record["result"]; return record["result"]
        attempt=root/f"attempt_{1+len(list(root.glob('attempt_*'))):03d}";attempt.mkdir()
        print(f"{name}: running ({attempt})",flush=True)
        try:
            with (attempt/"stage.log").open("w") as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
                result=fn(attempt)
            artifacts={str(p):sha256(p) for p in attempt.rglob("*") if p.is_file()}
            dump(checkpoint,dict(result=result,artifacts=artifacts))
            self.results[name]=result
            print(f"{name}: completed",flush=True)
            return result
        except Exception as error:
            failure=dict(error=type(error).__name__+": "+str(error),traceback=traceback.format_exc(),attempt=str(attempt))
            dump(attempt/"failure.json",failure); self.failures[name]=failure
            print(f"{name}: FAILED: {error}",flush=True)
            return None

    def clean(self,phase,directory):
        s=self.spec["phases"][phase]
        source=self.path(s["cif"])
        shutil.copy2(source,directory/"uploaded_original.cif")
        if s.get("preparation")=="public_mil53_hfree_4x4x4":
            from prepare_public_mil53 import reconstruct
            cif,report=reconstruct(directory/"uploaded_original.cif",directory/"clean",self.path(self.spec["clean_cif_auditor"]))
        else:
            # Generic chemically reviewed clean-cif output. Never set 'passed'
            # automatically for unknown solvent, disorder or activation state.
            report=self.path(s["cleaning_report"]);structure_gate(source,report)
            cif=directory/"host.cif";shutil.copy2(source,cif)
            r=json.loads(report.read_text())
            for field in ("source_audit","candidate_audit"):
                audit=self.path(str(report.parent/r[field]["path"]))
                destination=directory/(field+".json");shutil.copy2(audit,destination)
                r[field]["path"]=destination.name
            report=directory/"host.cleaning.json";dump(report,r)
        structure_gate(cif,report)
        return dict(cif=str(cif),cleaning_report=str(report),formula_units=float(s["formula_units"]),name=phase)

    def guest(self,directory):
        from rdkit import Chem
        g=self.spec["guest"]
        if g.get("geometry_sdf"):
            selected=directory/"guest_supplied.sdf"
            shutil.copy2(self.path(g["geometry_sdf"]),selected)
            if not g.get("geometry_source"):raise ValueError("supplied geometry requires geometry_source")
            if g.get("refinement") not in (None,"none"):raise ValueError("supplied geometry cannot also request refinement")
        else:
            prepare(g["smiles"],directory/"optimized",g.get("conformers",20),g.get("seed",20260913))
            selected=directory/"optimized/guest.sdf"
        if not g.get("geometry_sdf") and g.get("refinement")=="gfn2-xtb":
            self.command([sys.executable,SCRIPTS/"refine_guest_tblite.py","--output-dir",directory/"optimized"],directory)
            selected=directory/"optimized/guest_refined.sdf"
        elif g.get("refinement") not in (None,"none"):raise ValueError("unsupported guest refinement")
        mol=next(iter(Chem.SDMolSupplier(str(selected),removeHs=False)))
        if mol is None or Chem.MolToSmiles(Chem.RemoveHs(mol))!=Chem.MolToSmiles(Chem.MolFromSmiles(g["smiles"])):
            raise ValueError("guest geometry identity differs from requested SMILES")
        dump(directory/"geometry_validation.json",validate_geometry(mol))
        ff=json.loads(self.path(self.spec["force_field"]).read_text())
        types=g.get("atom_types")
        if types is None and g["name"]=="CO2" and Chem.MolToSmiles(Chem.RemoveHs(mol))=="O=C=O":
            types=[a.GetSymbol().capitalize()+"_co2" for a in mol.GetAtoms()]
        if types is None or len(types)!=mol.GetNumAtoms(): raise ValueError("explicit RDKit-order atom_types required")
        pseudo={a["name"]:a for a in ff["PseudoAtoms"]}
        if any(t not in pseudo or pseudo[t]["element"]!=a.GetSymbol() for t,a in zip(types,mol.GetAtoms())):
            raise ValueError("guest atom type/element mismatch")
        template=json.loads(self.path(g["molecule_json"]).read_text()) if g.get("molecule_json") else dict(CO2_JSON) if g["name"]=="CO2" else None
        if template is None: raise ValueError("custom gas requires molecule_json with sourced thermophysical metadata")
        if template.get("Type")!="rigid": raise ValueError("current optimized-geometry bridge supports rigid guests only")
        xyz=mol.GetConformer().GetPositions();xyz-=xyz.mean(axis=0)
        template["pseudoAtoms"]=[[t,x.tolist()] for t,x in zip(types,xyz)]
        model=directory/(g["name"]+".json");dump(model,template)
        dump(directory/"geometry_model_provenance.json",dict(force_field_source=self.spec["force_field_source"],
             guest_source=g["source"],geometry="RDKit conformer plus requested refinement, explicit atom mapping",selected_geometry=str(selected),
             warning="Geometry changed from original rigid model; adsorption parameter transfer needs validation"))
        return dict(molecule_json=str(model))

    def charge(self,state,directory):
        if self.spec.get("charges",{}).get("method")=="source":
            structure_gate(state["cif"],state["cleaning_report"],charges=True)
            if not self.spec["charges"].get("source"):raise ValueError("source charges require provenance")
            return state
        args=[sys.executable,SCRIPTS/"run_pacmof.py",state["cif"],"--output-dir",directory/"charged",
              "--cleaning-reports",state["cleaning_report"]]
        if self.spec.get("pacmof_root"): args += ["--pacmof-root",self.path(self.spec["pacmof_root"])]
        self.command(args,directory)
        target=directory/"charged"/(Path(state["cif"]).stem+"_pacmof_charged.cif")
        report=target.with_suffix(".cleaning.json");structure_gate(target,report,charges=True)
        return dict(state,cif=str(target),cleaning_report=str(report))

    def host(self,states,directory):
        method=self.spec["host"].get("method","phonons")
        if method=="elastic-network":
            quick_spec=dict(states=states,calculator=self.spec["host"]["calculator"],
                            temperature_K=self.spec["temperature_K"],cpu_threads=self.spec.get("cpu_threads",2),
                            mesh=self.spec["host"].get("proxy_mesh",[3,3,3]),
                            parameter_sets=self.spec["host"].get("proxy_parameter_sets"))
            result=quick_host_run(quick_spec,self.base)
            dump(directory/"host_thermodynamics.json",result)
            return result
        if method!="phonons":raise ValueError("host.method must be elastic-network or phonons")
        from advanced_free_energy import run as advanced_run
        spec=dict(states=states,calculator=self.spec["host"]["calculator"],phonons=self.spec["host"]["phonons"])
        spec["phonons"]=dict(spec["phonons"],temperatures=[self.spec["temperature_K"]])
        result=advanced_run("phonons",spec,self.base)
        dump(directory/"host_thermodynamics.json",result)
        return result

    def raspa(self,state,guest,directory):
        from ase.io import read
        r=self.spec["raspa"]; a=read(state["cif"])
        heights=1/np.linalg.norm(np.linalg.inv(a.cell.array),axis=0)
        cells=np.maximum(1,np.ceil(2*r["cutoff_A"]/heights).astype(int)).tolist()
        rows=[]
        for seed in r["seeds"]:
            args=[sys.executable,SCRIPTS/"prepare_raspa_inputs.py","--cif",state["cif"],
                  "--cleaning-report",state["cleaning_report"],"--phase",state["name"],
                  "--force-field",self.path(self.spec["force_field"]),"--molecule-json",guest["molecule_json"],
                  "--temperature",self.spec["temperature_K"],"--pressures",*self.spec["pressure_bar"],
                  "--unit-cells",*cells,"--cutoff",r["cutoff_A"],"--reservoir",self.spec.get("reservoir","ideal"),
                  "--output-dir",directory/f"seed_{seed}","--cycles",r["cycles"],"--seed",seed,
                  "--initialization-cycles",r["initialization_cycles"],"--equilibration-cycles",r["equilibration_cycles"]]
            if self.spec.get("gas_profile"): args += ["--gas-profile",self.path(self.spec["gas_profile"])]
            self.command(args,directory)
            for p in self.spec["pressure_bar"]:
                point=directory/f"seed_{seed}"/f"P_{pressure_tag(p)}bar"
                env=os.environ.copy()
                binary=self.path(r["binary"])
                env["RASPA_DIR"]=str(self.path(r["raspa_dir"])) if r.get("raspa_dir") else str(binary.parent.parent)
                with (point/"run.log").open("w") as log:
                    subprocess.run([str(binary)],cwd=point,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,
                                   timeout=r.get("point_timeout_seconds",1800))
                record=parse_point(point,self.spec["guest"]["name"],state["formula_units"],cells,seed)
                reservoir=json.loads((point/"preflight.json").read_text())["reservoir_point"]
                rows.append(dict(record,pressure_bar=p,fugacity_bar=reservoir["fugacity_bar"]))
                dump(directory/"partial_isotherm.json",rows)
        return dict(name=state["name"],rows=rows,unit_cells=cells,volume_A3_fu=a.get_volume()/state["formula_units"])

    def analyze(self,host,ads,directory):
        p=self.spec["pressure_bar"];t=self.spec["temperature_K"]
        states=[];fits={};blocks=[]
        for name in ("closed","open"):
            data=ads[name];n=[];u=[];f=[];error=[]
            model=self.spec.get("adsorption_models",{}).get(name)
            for pressure in p:
                rs=[r for r in data["rows"] if r["pressure_bar"]==pressure]
                ns=np.array([r["loading_molecules_fu"] for r in rs]);n.append(float(ns.mean()))
                u.append(float(np.mean([r["host_guest_energy_kJ_mol_fu"] for r in rs])))
                f.append(rs[0]["fugacity_bar"])
                # RASPA confidence limits are not 1-sigma standard errors.
                error.append(float(max([r["loading_confidence_molecules_fu"] for r in rs]+[np.ptp(ns)/2])))
                if not model and not all(r["accessibility_passed"] for r in rs): blocks.append(f"{name}@{pressure}: insertion/accessibility unresolved")
            fits[name]=fit_isotherm(f,n,u)
            henry=np.array(n[:2])/np.array(f[:2])
            # Conservative numeric check plus explicit run-specific reviewer evidence.
            stable=bool(min(henry)>0 and max(henry)/min(henry)<1.25)
            evidence=self.spec.get("henry_review",{}).get(name)
            if not model and (not stable or not evidence): blocks.append(f"{name}: first two Henry ratios/evidence not confirmed")
            states.append(dict(name=name,volume_A3_fu=data["volume_A3_fu"],fugacity_bar=f,
                loading_molecules_fu=n,loading_confidence_molecules_fu=error,henry_confirmed=stable and bool(evidence),
                host_source="same-model harmonic host calculation; not fitted to gate pressure"))
            if model:
                states[-1]["adsorption_model"]=model
                try:
                    state_adsorption_omega(states[-1],f,t)
                    if model.get("kind")=="nonadsorbing" and any(x>0 for x in n):
                        blocks.append(f"{name}: nonadsorbing model contradicts positive sampled loading")
                except (ValueError,KeyError,TypeError) as exc:
                    blocks.append(f"{name}: invalid analytical adsorption scenario: {exc}")
        dump(directory/"isotherm_and_energy_fits.json",fits)
        df=None
        proxy=bool(host and host.get("method")=="static_energy_plus_elastic_network_entropy_proxy")
        if host is None: blocks.append("host free energy unavailable")
        elif proxy:
            df=host["parameter_scenarios"][0]["delta_F_host_proxy_kJ_mol_fu"]
            for s in states:s["host_source"]="Static ML energy plus a DIFFERENT elastic-network entropy model; unvalidated fast screening approximation"
        else:
            difference=host["harmonic_state_differences"][0]
            if difference["status"]!="harmonic_estimate": blocks.append(difference["status"])
            else: df=difference["thermal"][0]["delta_F_host_harmonic_kJ_mol_fu"]
        result=dict(status="blocked" if blocks else "conditional_harmonic_screen",blocking_reasons=sorted(set(blocks)),
                    delta_F_open_minus_closed_kJ_mol_fu=df,states=states,fits=fits,
                    experimental_gateopen_pressure_bar=None,
                    limitations=["Harmonic model/phonon, pressure-grid and seed convergence are not established by a smoke run",
                                 "No kinetic barrier or hysteresis: only equilibrium crossings can be predicted",
                                 "Loading confidence and inter-seed spread are retained, not mislabeled as a validated pressure confidence interval"])
        if proxy:
            if not blocks:result["status"]="conditional_elastic_network_screen"
            result["host_free_energy_method"]="elastic_network_proxy"
            result["host_parameter_sensitivity_range_kJ_mol_fu"]=host["delta_F_proxy_scenario_range_kJ_mol_fu"]
            result["limitations"]+=host["warnings"]
        if any(s.get("adsorption_model") for s in states):
            result["limitations"].append("Analytical zero/Henry adsorption branches are explicitly assumed scenarios; Widom alone cannot bound cooperative finite-pressure adsorption")
        if not blocks:
            states[0]["host_F_kJ_mol_fu"]=0;states[1]["host_F_kJ_mol_fu"]=df
            profile=load_gas_profile(self.path(self.spec["gas_profile"]) if self.spec.get("gas_profile") else None)
            # Evaluate the actual integrated isotherms between sampled nodes;
            # a dense evaluation grid does not create new adsorption evidence.
            evaluation_p=np.unique(np.r_[p,np.geomspace(p[0],p[-1],401)]).tolist()
            reservoir=[reservoir_point(t,x,backend=self.spec.get("reservoir","ideal"),fluid=profile.get("coolprop_fluid",profile["name"])) for x in evaluation_p]
            result["equilibrium"]=phase_diagram(dict(temperature_K=t,pressure_bar=evaluation_p,states=states),reservoir)
            result["equilibrium_opening_pressures_bar"]=[r["pressure_bar"] for r in result["equilibrium"]["crossings"]
                if r.get("globally_stable") and r.get("stable_before")=="closed" and r.get("stable_after")=="open"]
            model=states[0].get("adsorption_model",{})
            if model.get("kind")=="nonadsorbing" and model.get("sensitivity_K_molecules_fu_bar"):
                result["closed_adsorption_sensitivity"]=closed_adsorption_sensitivity(
                    dict(temperature_K=t,pressure_bar=evaluation_p,states=states),reservoir,
                    model["sensitivity_K_molecules_fu_bar"])
            if proxy:
                from copy import deepcopy
                sensitivity=[]
                for scenario in host["parameter_scenarios"]:
                    varied=deepcopy(states)
                    varied[1]["host_F_kJ_mol_fu"]=scenario["delta_F_host_proxy_kJ_mol_fu"]
                    phase=phase_diagram(dict(temperature_K=t,pressure_bar=evaluation_p,states=varied),reservoir)
                    sensitivity.append(dict(name=scenario["name"],delta_F_host_kJ_mol_fu=varied[1]["host_F_kJ_mol_fu"],
                        equilibrium_opening_pressures_bar=[r["pressure_bar"] for r in phase["crossings"]
                            if r.get("globally_stable") and r.get("stable_before")=="closed" and r.get("stable_after")=="open"]))
                result["host_parameter_pressure_sensitivity"]=dict(scenarios=sensitivity,interpretation="parameter scenarios, not confidence interval")
        dump(directory/"result.json",result)
        return result

    def run(self):
        import torch
        torch.set_num_threads(self.spec.get("cpu_threads",2))
        if set(self.spec["phases"])!={"closed","open"}: raise ValueError("exactly closed/open phases required")
        pressures=np.asarray(self.spec["pressure_bar"],dtype=float)
        if len(pressures)<3 or np.any(np.diff(pressures)<=0): raise ValueError("at least 3 increasing pressure points")
        seeds=self.spec["raspa"]["seeds"]
        if not seeds or len(seeds)!=len(set(seeds)): raise ValueError("unique independent seeds required")
        profile=load_gas_profile(self.path(self.spec["gas_profile"]) if self.spec.get("gas_profile") else None)
        if profile["name"]!=self.spec["guest"]["name"]: raise ValueError("gas profile/name mismatch")
        pressure_gate(profile,self.spec["temperature_K"],pressures)
        for p in pressures: reservoir_point(self.spec["temperature_K"],p,backend=self.spec.get("reservoir","ideal"),fluid=profile.get("coolprop_fluid",profile["name"]))
        guest=self.stage("guest",self.guest)
        clean={k:self.stage("clean_"+k,lambda d,k=k:self.clean(k,d)) for k in ("closed","open")}
        states=[];ads={}
        for name,s in clean.items():
            if s is None: continue
            if self.spec["host"].get("relax",False):
                if self.spec["phases"][name].get("chemistry_adapter")!="mil53al":
                    self.failures["relax_"+name]="No chemistry re-audit adapter: upload relaxed, reviewed CIF and set relax=false"
                    continue
                def relaxation(d,s=s):
                    from relax_host import relax
                    calc=dict(self.spec["host"]["calculator"])
                    if calc.get("model") and self.path(calc["model"]).is_file():calc["model"]=str(self.path(calc["model"]))
                    relax(Path(s["cif"]),Path(s["cleaning_report"]),d/"relaxed",calc,
                          self.path(self.spec["clean_cif_auditor"]),self.spec["host"].get("fmax",.02),self.spec["host"].get("steps",500))
                    return dict(s,cif=str(d/"relaxed/host.cif"),cleaning_report=str(d/"relaxed/host.cleaning.json"))
                s=self.stage("relax_"+name,relaxation)
                if s is None: continue
            states.append(s)
            charged=self.stage("charge_"+name,lambda d,s=s:self.charge(s,d))
            if charged and guest: ads[name]=self.stage("raspa_"+name,lambda d,s=charged:self.raspa(s,guest,d))
        host=None
        if len(states)==2:
            # Reject inequivalent compositions before expensive host calculations.
            from ase.io import read
            comps=[{k:v/s["formula_units"] for k,v in Counter(read(s["cif"]).get_chemical_symbols()).items()} for s in states]
            if comps[0]!=comps[1]: raise ValueError("phase composition/f.u. mismatch")
            host=self.stage("host_free_energy",lambda d:self.host(states,d))
        analysis=None
        if set(ads)=={"closed","open"} and all(ads.values()):
            analysis=self.stage("analysis",lambda d:self.analyze(host,ads,d))
        summary=dict(status="completed_conditional_screen" if analysis and analysis["status"]!="blocked" else "incomplete_prediction",
                     completed_stages=list(self.results),failures=self.failures,analysis=analysis)
        dump(self.out/"workflow_summary.json",summary)
        return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec",required=True,type=Path);parser.add_argument("--output-dir",required=True,type=Path)
    args=parser.parse_args()
    summary=Workflow(args.spec,args.output_dir).run()
    print(json.dumps(dict(status=summary["status"],report=str(args.output_dir/"workflow_summary.json")),indent=2))
    sys.exit(0 if summary["status"]=="completed_conditional_screen" else 2)
