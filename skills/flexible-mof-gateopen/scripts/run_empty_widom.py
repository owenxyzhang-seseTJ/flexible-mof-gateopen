#!/usr/bin/env python3
"""Repeat empty-host Widom in fresh directories from hash-verified RASPA inputs.

No geometry edits, relaxation, MD or DFT. All guest moves except Widom are
disabled. A dilute Henry coefficient is not a bound on finite-pressure uptake.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from parse_raspa_output import parse_point


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def run(source, checkpoint, binary, out, seeds, cycles, formula_units, timeout=600):
    source, checkpoint, binary, out = map(lambda x: Path(x).resolve(), (source, checkpoint, binary, out))
    if len(seeds)<2 or len(set(seeds))!=len(seeds) or cycles<100:
        raise ValueError("at least two unique seeds and 100 production cycles required")
    if not math.isfinite(formula_units) or formula_units<=0:
        raise ValueError("positive formula units per CIF required")
    identity=json.loads(checkpoint.read_text())["artifacts"]
    spec=json.loads((source/"simulation.json").read_text())
    if spec["SimulationType"]!="MonteCarlo" or len(spec["Systems"])!=1 or len(spec["Components"])!=1:
        raise ValueError("single fixed-host single-component MonteCarlo template required")
    system=spec["Systems"][0]; component=spec["Components"][0]
    if system["Type"]!="Framework" or system.get("UseChargesFrom")!="CIF_File":
        raise ValueError("charged framework template required")
    files=["simulation.json","preflight.json","force_field.json",system["Name"]+".cif",component["Name"]+".json"]
    hashes={}
    for name in files:
        path=source/name
        if identity.get(str(path))!=digest(path):
            raise ValueError(f"checkpoint input hash mismatch: {path}")
        hashes[str(path)]=digest(path)
    preflight=json.loads((source/"preflight.json").read_text())
    if preflight["structure"]["status"]!="passed" or preflight["structure"]["candidate_sha256"]!=digest(source/(system["Name"]+".cif")):
        raise ValueError("prepared CIF differs from passed clean-cif preflight")
    spec.update(NumberOfCycles=cycles,NumberOfInitializationCycles=0,
                NumberOfEquilibrationCycles=0,PrintEvery=max(1,cycles//10))
    for key in list(component):
        if key.endswith("Probability"):
            component[key]=0.
    component.update(WidomProbability=1.,CreateNumberOfMolecules=0)
    spec.pop("RestartFromBinary",None)
    out.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env["RASPA_DIR"]=str(binary.parent.parent)
    rows=[];start=time.monotonic()
    for seed in seeds:
        directory=out/f"seed_{seed}";directory.mkdir()
        for name in files:
            if name!="simulation.json":shutil.copy2(source/name,directory/name)
        spec["RandomSeed"]=seed;dump(directory/"simulation.json",spec)
        print(f"empty-host Widom seed {seed}: {cycles} cycles",flush=True)
        with (directory/"run.log").open("w") as log:
            subprocess.run([str(binary)],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,
                           check=True,timeout=timeout)
        row=parse_point(directory,component["Name"],formula_units,system["NumberOfUnitCells"],seed)
        if row["loading_molecules_fu"]!=0 or row["widom_henry_molecules_fu_bar"] is None:
            raise ValueError("empty-host occupancy or Henry-output audit failed")
        if not sum(row["raw"]["move_diagnostic"]["widom_constructed"] or []):
            raise ValueError("no constructed Widom trials: audit setup before estimating Henry")
        values=row["widom_henry_molecules_fu_bar"]
        if not all(math.isfinite(x) and x>=0 for x in values):
            raise ValueError("invalid Widom statistics")
        rows.append(row);dump(out/"partial.json",rows)
        print(f"seed {seed}: K={values[0]:.8g}, native confidence={values[1]:.8g}",flush=True)
    # Merely a reproducible scale for subsequent sensitivity, not an upper bound.
    scale=max(sum(r["widom_henry_molecules_fu_bar"]) for r in rows)
    result=dict(method="empty_framework_Widom_only",temperature_K=system["ExternalTemperature"],
                gas=component["Name"],formula_units_per_CIF=formula_units,cycles_per_seed=cycles,
                source_inputs_sha256=hashes,checkpoint_sha256=digest(checkpoint),
                binary_sha256=digest(binary),rows=rows,elapsed_seconds=time.monotonic()-start,
                sensitivity_K_scale_molecules_fu_bar=scale,
                interpretation="max(seed mean + native confidence) is a heuristic sensitivity scale, NOT a certified upper bound or standard deviation",
                limitations=["Finite random sampling may miss rare favorable cavities/orientations",
                             "Henry is a zero-density property; cooperative filling cannot be excluded",
                             "Zero occupancy is enforced by the Widom-only algorithm, not evidence of physical zero adsorption"])
    dump(out/"result.json",result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("source-point","checkpoint","binary","out"):
        parser.add_argument("--"+name,required=True,type=Path)
    parser.add_argument("--seeds",type=int,nargs="+",default=[1729,3253,5801,7919])
    parser.add_argument("--cycles",type=int,default=10000)
    parser.add_argument("--formula-units",type=float,required=True)
    parser.add_argument("--timeout",type=int,default=600)
    args=parser.parse_args()
    run(args.source_point,args.checkpoint,args.binary,args.out,args.seeds,args.cycles,args.formula_units,args.timeout)
