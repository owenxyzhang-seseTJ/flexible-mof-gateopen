#!/usr/bin/env python3
"""Advanced ML, harmonic/QHA, constrained MC-TI and osmotic-landscape workflows.

Use --help and references/advanced_workflows.md for strict JSON schemas.
Paths inside a specification are relative to that JSON file.
"""
import argparse
import contextlib
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
try:
    from flexible_mof_gateopen import advanced_thermo as thermo
    from flexible_mof_gateopen import host_models, framework_ti
    from flexible_mof_gateopen.preflight import structure_gate, load_gas_profile, pressure_gate
    from flexible_mof_gateopen.reservoir import reservoir_point
except ModuleNotFoundError as error:
    if error.name != "flexible_mof_gateopen": raise
    import advanced_thermo as thermo
    import host_models, framework_ti
    from preflight import structure_gate, load_gas_profile, pressure_gate
    from reservoir import reservoir_point


def run(command,spec,base):
    def resolved(path): return (base/path).resolve()
    def reservoir(section):
        profile=load_gas_profile(resolved(section["gas_profile"]) if section.get("gas_profile") else None)
        pressure_gate(profile,section["temperature_K"],section["pressure_bar"])
        backend=section.get("reservoir_backend","ideal")
        # Explicit profile mapping for custom gases, not guessed aliases.
        fluid=profile.get("coolprop_fluid",profile["name"])
        return [reservoir_point(section["temperature_K"],p,backend=backend,fluid=fluid)
                for p in section["pressure_bar"]]
    if command=="reservoir":
        return dict(points=reservoir(spec))
    if command=="phase-diagram":
        return dict(slices=[thermo.phase_diagram(s,reservoir(s)) for s in spec["slices"]])
    if command=="qha":
        rows=[]
        for s in spec["slices"]:
            if not s.get("all_phonons_stable") or not s.get("source"):
                raise ValueError("QHA requires stable phonons at every volume and provenance")
            rows.append(dict(phase=s["phase"],temperature_K=s["temperature_K"],
                             pressure_bar=s.get("pressure_bar",0),
                             **thermo.qha_minimize(s["volume_A3_fu"],s["F_total_kJ_mol_fu"],s.get("pressure_bar",0))))
        return dict(slices=rows,warning="QHA only within each stable basin; do not fit across breathing phases")
    if command=="landscape":
        if not spec.get("landscape_source") or not spec.get("coordinate_measure"):
            raise ValueError("landscape source and coordinate measure/Jacobian convention required")
        if "states" in spec:
            # Construct Omega(q,P) from actual host F(q), V(q), and N(q,f).
            q=[s["q"] for s in spec["states"]]
            result=thermo.phase_diagram(spec,reservoir(spec))
            curves=thermo.np.array(result["omega_kJ_mol_fu"]).T.tolist()
            return dict(thermo.landscape_scan(q,spec["pressure_bar"],curves),osmotic=result,
                        coordinate_measure=spec["coordinate_measure"],source=spec["landscape_source"])
        if not spec.get("analytic_test",False):
            # Direct imported Omega must still respect vapor pressure bounds.
            reservoir(spec)
        return dict(thermo.landscape_scan(spec["q"],spec["pressure_bar"],spec["omega_kJ_mol_fu"]),
                    source=spec["landscape_source"],coordinate_measure=spec["coordinate_measure"])
    if command in {"host-energy","phonons","mc-ti"}:
        from ase.io import read
        # Audit every state BEFORE loading weights or calculating anything.
        states=spec.get("states",[])
        if not states: raise ValueError("states required")
        audits=[structure_gate(resolved(s["cif"]),resolved(s["cleaning_report"])) for s in states]
        atoms_list=[read(str(resolved(s["cif"]))) for s in states]
        compositions=[]
        for s,a in zip(states,atoms_list):
            z=float(s["formula_units"])
            if not math.isfinite(z) or z<=0: raise ValueError("positive formula_units required")
            compositions.append({k:v/z for k,v in Counter(a.get_chemical_symbols()).items()})
        reference=compositions[0]
        if any(set(c)!=set(reference) or any(not math.isclose(c[k],reference[k],abs_tol=1e-8) for k in c)
               for c in compositions):
            raise ValueError("all states must share the same composition per formula unit")
        ml=spec.get("calculator",{})
        if ml.get("model") and (base/ml["model"]).is_file():
            ml=dict(ml,model=str(resolved(ml["model"])))
        calc,metadata=host_models.calculator(**ml)
        results=[]
        for s,a,audit in zip(states,atoms_list,audits):
            a.calc=calc
            if command=="host-energy":
                r=host_models.energy_record(a,s["formula_units"])
            elif command=="phonons":
                r=host_models.phonon_free_energy(a,s["formula_units"],**spec["phonons"])
            else:
                r=framework_ti.framework_mc_ti(a,s["formula_units"],**spec["mc"])
            results.append(dict(name=s["name"],audit=audit,**r))
        if command=="host-energy": results=host_models.energy_difference(results)
        differences=host_models.harmonic_state_differences(results) if command=="phonons" else None
        return dict(calculator=metadata,states=results,
                    harmonic_state_differences=differences,
                    classification="software_calculation_material_accuracy_unvalidated",
                    warning="Static delta E is not finite-T host delta F; validate dispersion, relative phase energies and basin identity")
    raise ValueError("unknown command")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",choices=["host-energy","phonons","qha","mc-ti","reservoir","phase-diagram","landscape"])
    p.add_argument("--spec",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.out.exists(): p.error("output exists; choose a new file")
    try:
        spec=json.loads(args.spec.read_text())
        # Keep third-party ML progress messages away from machine-readable outputs.
        with contextlib.redirect_stdout(sys.stderr):
            result=run(args.command,spec,args.spec.resolve().parent)
        result.update(spec_sha256=hashlib.sha256(args.spec.read_bytes()).hexdigest(),command=args.command)
        text=json.dumps(result,indent=2,allow_nan=False)+"\n"
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open("x") as f: f.write(text)
        print(args.out)
    except (ValueError,KeyError,ImportError) as error:
        p.error(str(error))


if __name__=="__main__": main()
