#!/usr/bin/env python3
"""Fixed-cell host relaxation, preserving phase volume and reviewed connectivity."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read, write
from ase.optimize import FIRE
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen.host_models import calculator, energy_record
    from flexible_mof_gateopen.preflight import structure_gate
except ModuleNotFoundError:
    from host_models import calculator, energy_record
    from preflight import structure_gate
from prepare_public_mil53 import chemistry, certify, dump


def relax(cif, report, out, calculator_spec, auditor, fmax=.02, steps=500):
    structure_gate(cif,report)
    a=read(str(cif)); before=chemistry(a)
    out.mkdir(parents=True,exist_ok=False)
    a.calc,meta=calculator(**calculator_spec)
    initial=energy_record(a,before["formula_units"])
    opt=FIRE(a,logfile=str(out/"relax.log"),trajectory=str(out/"relax.traj"),maxstep=.05)
    converged=bool(opt.run(fmax=fmax,steps=steps))
    a.wrap(); target=out/"host.cif"; write(str(target),a)
    final=energy_record(a,before["formula_units"])
    result=dict(calculator=meta,initial=initial,final=final,converged=converged,
                steps=opt.nsteps,cell_policy="fixed cell; Helmholtz F at supplied phase volume, not relaxed Gibbs F")
    dump(out/"relaxation.json",result)
    if not converged: raise ValueError("host relaxation did not converge; see relaxation.json")
    after=chemistry(a)
    if before["bond_indices"]!=after["bond_indices"]:
        raise ValueError("host connectivity changed; needs new chemistry/phase review")
    certify(cif,target,auditor,"Fixed-cell FIRE position optimization; same material-specific periodic bond graph, composition and cell; no charge transfer from old geometry",result)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cif",required=True,type=Path);p.add_argument("--cleaning-report",required=True,type=Path)
    p.add_argument("--output-dir",required=True,type=Path);p.add_argument("--model",required=True)
    p.add_argument("--clean-cif-auditor",required=True,type=Path)
    p.add_argument("--head",default=None)
    p.add_argument("--steps",type=int,default=500);p.add_argument("--fmax",type=float,default=.02)
    args=p.parse_args()
    # Limit CPU oversubscription for small periodic cells.
    import torch
    torch.set_num_threads(2)
    relax(args.cif,args.cleaning_report,args.output_dir,dict(backend="mace",model=args.model,head=args.head),args.clean_cif_auditor,args.fmax,args.steps)
