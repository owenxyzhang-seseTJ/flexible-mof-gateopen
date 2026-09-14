#!/usr/bin/env python3
"""Reproducible analytic and optional real-backend checks, NOT a MOF benchmark."""
import argparse
import json
import math
from pathlib import Path
import sys
import contextlib

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen import advanced_thermo as a, host_models as h, framework_ti as ti
    from flexible_mof_gateopen.reservoir import reservoir_point
except ModuleNotFoundError as error:
    if error.name != "flexible_mof_gateopen": raise
    import advanced_thermo as a, host_models as h, framework_ti as ti
    from reservoir import reservoir_point


def validate(ml=False):
    np=a.np
    q=np.linspace(-2,2,201); tilt=np.linspace(0,2,41)
    land=a.landscape_scan(q,np.exp(tilt),[(q*q-1)**2-x*q for x in tilt])
    numeric=math.log(land["spline_spinodals"][0]["pressure_bar"])
    exact=8/(3*math.sqrt(3))
    assert abs(numeric-exact)<.003
    mc=ti.mc_ti(lambda x: float(x@x),3,spring_eV_A2=1,bound_A=2,step_A=.3,
                lambdas=np.linspace(0,1,21),warmup=500,samples=2500,stride=3,seed=32)
    exact_ti=1.5*ti.KB*300*math.log(2)
    assert abs(mc["delta_F_TI_eV"]-exact_ti)<.0025
    reservoir=reservoir_point(268,19,backend="coolprop")
    report=dict(classification="analytic_and_Si_software_checks_not_MOF_prediction",
                reservoir_CO2_268K_19bar=reservoir,
                double_well_spinodal=dict(numeric_tilt=numeric,exact_tilt=exact,relative_error=abs(numeric/exact-1)),
                constrained_TI_harmonic=dict(numeric_delta_F_eV=mc["delta_F_TI_eV"],exact_delta_F_eV=exact_ti,
                                             approximate_block_SE_eV=mc["approximate_block_SE_eV"],seed=32),
                qha_parabola=a.qha_minimize(np.linspace(90,110,21),(np.linspace(90,110,21)-100)**2))
    if ml:
        from ase.build import bulk
        records=[]
        for backend in ("chgnet","mace"):
            calc,metadata=h.calculator(backend)
            phases=[]
            for lattice in (5.43,5.48):
                atoms=bulk("Si","diamond",a=lattice);atoms.calc=calc
                phases.append(dict(lattice_A=lattice,**h.energy_record(atoms,2)))
            entry=dict(calculator=metadata,static_strain_check=h.energy_difference(phases))
            if backend=="chgnet":
                atoms=bulk("Si","diamond",a=5.43);atoms.calc=calc
                entry["phonons"]=h.phonon_free_energy(atoms,2,[0,300],supercell=[2,2,2],mesh=[4,4,4])
                assert entry["phonons"]["status"]=="harmonic_screen"
                # Short wiring check only, not a converged atomistic TI result.
                entry["framework_MC_wiring_only"]=ti.framework_mc_ti(atoms,2,
                    lambdas=[0,.5,1],warmup=5,samples=20,stride=1,blocks=2,seed=32)
            records.append(entry)
        report["actual_backend_checks"]=records
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ml",action="store_true")
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.out.exists(): p.error("output exists")
    with contextlib.redirect_stdout(sys.stderr): result=validate(args.ml)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("x") as f: f.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(args.out)


if __name__=="__main__": main()
