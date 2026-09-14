#!/usr/bin/env python3
"""Fast static-force-field + elastic-network entropy screening, without MD/DFT."""
import argparse
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen.preflight import structure_gate,sha256
    from flexible_mof_gateopen.host_models import calculator,energy_record,energy_difference
    from flexible_mof_gateopen.elastic_network import compare_states
except ModuleNotFoundError:
    from preflight import structure_gate,sha256
    from host_models import calculator,energy_record,energy_difference
    from elastic_network import compare_states


def run(spec,base):
    from ase.io import read
    start=time.monotonic()
    states=spec["states"]
    if [s["name"] for s in states]!=["closed","open"]:raise ValueError("ordered states closed, open required")
    audits=[structure_gate(base/s["cif"],base/s["cleaning_report"]) for s in states]
    atoms=[read(str(base/s["cif"])) for s in states]
    z=[s["formula_units"] for s in states]
    if spec.get("static_energy_json"):
        cached=json.loads((base/spec["static_energy_json"]).read_text())
        records=cached["states"];metadata=cached["calculator"]
        if len(records)!=2 or any(r["name"]!=s["name"] or r["audit"]["candidate_sha256"]!=a["candidate_sha256"] or r["formula_units"]!=u for r,s,a,u in zip(records,states,audits,z)):
            raise ValueError("cached static energies do not match exact CIFs and normalization")
    else:
        import torch
        torch.set_num_threads(spec.get("cpu_threads",2))
        model=dict(spec["calculator"])
        if model.get("model") and (base/model["model"]).is_file():model["model"]=str((base/model["model"]).resolve())
        calc,metadata=calculator(**model)
        records=[]
        for a,u,s,audit in zip(atoms,z,states,audits):
            a.calc=calc
            records.append(dict(name=s["name"],audit=audit,**energy_record(a,u)))
    records=energy_difference(records)
    result=compare_states(atoms,[r["energy_kJ_mol_fu"] for r in records],z,spec["temperature_K"],
                          mesh=spec.get("mesh",[3,3,3]),parameter_sets=spec.get("parameter_sets"))
    return dict(result,calculator=metadata,states=records,elapsed_seconds=time.monotonic()-start)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spec",required=True,type=Path);p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    if args.out.exists():p.error("use a new output path")
    result=run(json.loads(args.spec.read_text()),args.spec.resolve().parent)
    result["spec_sha256"]=sha256(args.spec)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for r in result["parameter_scenarios"]:
        print(r["name"],"dS",r["delta_S_proxy_J_mol_fu_K"],"dF",r["delta_F_host_proxy_kJ_mol_fu"])
    print("elapsed seconds:",result["elapsed_seconds"],args.out)
