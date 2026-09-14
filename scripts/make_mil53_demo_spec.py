#!/usr/bin/env python3
"""Create a local runnable specification from downloaded public MIL53 inputs.

No executable/checkpoint is downloaded or hidden in the JSON. Paths supplied
by the caller and public-source hashes are recorded in the output directory.
"""
import argparse
import json
from pathlib import Path
import sys
from prepare_public_mil53 import dump

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("closed-source","open-source","model","pacmof-root","raspa-bin","clean-cif-auditor","output-dir"):
        p.add_argument("--"+key,required=True,type=Path)
    p.add_argument("--cycles",type=int,default=1000)
    p.add_argument("--seeds",type=int,nargs="+",default=[417,941])
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    pseudo=[];interactions=[]
    rows=[("Al",True,"Al",26.9815385,0,254.0,4.008),
          ("C",True,"C",12.011,0,52.84,3.431),("O",True,"O",15.999,0,30.19,3.118),
          ("H",True,"H",1.008,0,22.14,2.571),
          ("C_co2",False,"C",12.011,.6512,29.933,2.745),
          ("O_co2",False,"O",15.999,-.3256,85.671,3.017)]
    for name,framework,element,mass,charge,epsilon,sigma in rows:
        pseudo.append(dict(name=name,framework=framework,element=element,mass=mass,charge=charge,print_to_output=True,print_as=element))
        interactions.append(dict(name=name,type="lennard-jones",parameters=[epsilon,sigma]))
    dump(a.output_dir/"force_field.json",dict(PseudoAtoms=pseudo,SelfInteractions=interactions,MixingRule="Lorentz-Berthelot",TruncationMethod="shifted",TailCorrections=False))
    spec=dict(temperature_K=268.,pressure_bar=[.001,.003,.03,.3,3.,20.],reservoir="coolprop",cpu_threads=2,
        clean_cif_auditor=str(a.clean_cif_auditor.resolve()),pacmof_root=str(a.pacmof_root.resolve()),
        phases={name:dict(cif=str(source.resolve()),formula_units=4,preparation="public_mil53_hfree_4x4x4",chemistry_adapter="mil53al")
                for name,source in (("closed",a.closed_source),("open",a.open_source))},
        guest=dict(name="CO2",smiles="O=C=O",refinement="gfn2-xtb",source="Garcia-Sanchez et al., J Phys Chem C 2009 113 8814; geometry refined by GFN2-xTB from RDKit, unvalidated screening transfer",conformers=10),
        force_field="force_field.json",force_field_source="UFF generic host LJ (Rappe et al., DOI 10.1021/ja00051a040), PACMOF charges and Garcia-Sanchez CO2 LJ/charges; NOT a validated MIL53 adsorption parameterization",
        host=dict(calculator=dict(backend="mace",model=str(a.model.resolve()),head="pbe_d3",dispersion=False),relax=True,fmax=.02,steps=500,
                  phonons=dict(supercell=[1,1,1],mesh=[2,2,2],displacement=.01,imaginary_tolerance_THz=.0001,force_tolerance=.03)),
        raspa=dict(binary=str(a.raspa_bin.resolve()),cutoff_A=8.,cycles=a.cycles,initialization_cycles=200,equilibration_cycles=300,seeds=a.seeds,point_timeout_seconds=1800),
        henry_review={},
        classification="real_backend_smoke_test; no automatic Henry review; 1x1x1 phonons and 8A GCMC cutoff require convergence")
    dump(a.output_dir/"workflow.json",spec)
    print(a.output_dir/"workflow.json")
