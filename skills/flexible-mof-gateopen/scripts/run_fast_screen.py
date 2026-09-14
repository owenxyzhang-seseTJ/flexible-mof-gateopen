#!/usr/bin/env python3
"""Agent-native fast-screen entry. No UI, MD, DFT or finite-displacement phonons."""
import argparse
import json
from pathlib import Path
import sys
from run_workflow import Workflow


def fast_spec(spec,root):
    """Resolve bundle defaults without replacing explicit user chemistry choices."""
    spec=json.loads(json.dumps(spec))
    host=spec.setdefault("host",{})
    if host.get("method","elastic-network")!="elastic-network" or host.get("relax",False):
        raise ValueError("fast mode requires elastic-network and pre-reviewed fixed geometries; no host relaxation/phonons")
    host.update(method="elastic-network",relax=False)
    host.setdefault("calculator",dict(backend="mace",model=str(root/"assets/models/mofs_v2.model"),head="pbe_d3",dispersion=False))
    if host["calculator"].get("backend")!="mace":
        raise ValueError("this packaged fast route supports MACE only")
    if host["calculator"].get("head")!="pbe_d3" or host["calculator"].get("dispersion"):
        raise ValueError("bundled MACE model requires head=pbe_d3 and no added dispersion")
    spec.setdefault("pacmof_root",str(root/"vendor/pacmof"))
    if spec.get("charges",{}).get("method","pacmof") not in {"pacmof","source"}:
        raise ValueError("charges.method must be pacmof or source")
    spec.setdefault("reservoir","coolprop")
    r=spec.setdefault("raspa",{})
    r.setdefault("binary",str(root/"vendor/raspa3/bin/raspa3"))
    for k,v in dict(cycles=20000,initialization_cycles=3000,equilibration_cycles=5000,
                    seeds=[417,941],cutoff_A=12,point_timeout_seconds=1800).items():r.setdefault(k,v)
    if spec["guest"].get("refinement") not in (None,"none"):
        raise ValueError("fast bundle uses RDKit or a sourced supplied geometry, without electronic refinement")
    if any(p.get("preparation") for p in spec["phases"].values()):
        raise ValueError("fast bundle requires clean-cif reviewed inputs, not a material-specific reconstruction adapter")
    model=spec.get("adsorption_models",{}).get("closed")
    if model and model.get("kind")=="nonadsorbing" and not model.get("sensitivity_K_molecules_fu_bar"):
        raise ValueError("zero closed scenario requires residual adsorption sensitivity coefficients")
    return spec


def run(spec_path,out):
    spec_path=Path(spec_path).resolve();out=Path(out).resolve()
    root=Path(__file__).resolve().parent.parent
    if (root/"MANIFEST.json").exists():
        from verify_bundle import verify
        checked=verify(root)
        if not checked["integrity_passed"] or not checked["runtime_dependencies_present_and_sklearn_compatible"]:
            raise ValueError("bundle integrity/runtime gate failed: "+json.dumps(checked))
    spec=fast_spec(json.loads(spec_path.read_text()),root)
    # Resolve existing data paths while still in the user's spec directory.
    def paths(value):
        if isinstance(value,dict):return {k:paths(v) for k,v in value.items()}
        if isinstance(value,list):return [paths(v) for v in value]
        if isinstance(value,str) and (spec_path.parent/value).exists():return str((spec_path.parent/value).resolve())
        return value
    spec=paths(spec)
    out.mkdir(parents=True,exist_ok=True)
    resolved=out/"resolved_spec.json"
    text=json.dumps(spec,indent=2)+"\n"
    if resolved.exists() and resolved.read_text()!=text:raise ValueError("changed configuration requires new output directory")
    if not resolved.exists():resolved.write_text(text)
    summary=Workflow(resolved,out/"run").run()
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spec",required=True,type=Path);p.add_argument("--output-dir",required=True,type=Path)
    a=p.parse_args();r=run(a.spec,a.output_dir)
    print(json.dumps(dict(status=r["status"],summary=str(a.output_dir/"run/workflow_summary.json")),indent=2))
    sys.exit(0 if r["status"]=="completed_conditional_screen" else 2)
