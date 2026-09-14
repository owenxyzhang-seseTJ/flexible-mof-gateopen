#!/usr/bin/env python3
"""Audit published MIL-53(Al) phonon entropy tables; not a new MOF calculation.

Source: Elena et al., 10.1038/s41524-025-01611-8, CC BY 4.0.
Download https://raw.githubusercontent.com/ddmms/data/main/mace-mof-0/phonons/MIL-53.zip
or GitHub blob 7c6ab0f7d40131240406b4c6b431aacbc45b4ca5 (raw API).
"""
import argparse
from collections import Counter
import hashlib
import io
import json
import lzma
from pathlib import Path
import zipfile

import numpy as np
import yaml

EXPECTED_SHA256="005e469f713b90f472672898b6bfc15329ee9cff12d7ce02e9bef50395eede2d"


def audit(path):
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest!=EXPECTED_SHA256:
        raise ValueError("archive version differs from audited source; inspect changes before using this script")
    rows=[]
    with zipfile.ZipFile(path) as archive:
        for name in ("0b","mof0","mof0-v2"):
            prefix=f"MIL-53/{name}/MIL-53"
            thermal=archive.read(prefix+"-thermal.dat")
            band=yaml.load(lzma.decompress(archive.read(prefix+"-auto_bands.yml.xz")),
                           Loader=getattr(yaml,"CSafeLoader",yaml.SafeLoader))
            composition=dict(Counter(x["symbol"] for x in band["points"]))
            if composition!={"Al":4,"O":20,"C":32,"H":20} or band["natom"]!=76:
                raise ValueError("unexpected primitive-cell composition")
            # 4 [Al(OH)(C8H4O4)] units in the phonon primitive cell.
            z=4
            frequencies=np.array([[b["frequency"] for b in q["band"]] for q in band["phonon"]])
            data=np.loadtxt(io.BytesIO(thermal))
            idx=np.flatnonzero(data[:,0]==300)
            if len(idx)!=1: raise ValueError("300 K row missing/duplicated")
            i=int(idx[0]);t,cv,f,s=data[i]
            slope=(data[i+1,2]-data[i-1,2])/(data[i+1,0]-data[i-1,0])
            predicted_s=-1000*slope
            # Header labels F as J/mol, but the numerical column follows
            # Phonopy's kJ/mol convention. Verify locally rather than trust it.
            if abs(predicted_s/s-1)>.001:
                raise ValueError("F/S unit cross-check failed")
            rows.append(dict(model_directory=name,temperature_K=float(t),
                primitive_composition=composition,formula_units=z,
                S_J_mol_primitive_K=float(s),S_J_mol_fu_K=float(s/z),
                Cv_J_mol_fu_K=float(cv/z),F_vib_kJ_mol_fu=float(f/z),
                min_band_frequency_THz=float(frequencies.min()),
                negative_band_entries_below_minus_1e_minus4=int((frequencies< -1e-4).sum()),
                harmonic_stability_on_sampled_path=bool(frequencies.min()>= -1e-4),
                source_header=thermal.decode().splitlines()[0],
                S_from_F_finite_difference_J_mol_primitive_K=float(predicted_s),
                entropy_consistency_relative_error=float(abs(predicted_s/s-1)),
                table_sha256=hashlib.sha256(thermal).hexdigest(),
                warning="Sampled band path is not a full Brillouin-zone stability proof; F is vibrational only"))
    return dict(source_doi="10.1038/s41524-025-01611-8",source_license="CC BY 4.0",
                source_url="https://github.com/ddmms/data/tree/main/mace-mof-0/phonons",
                archive_sha256=digest,results=rows,
                classification="reanalysis_of_published_single_structure_ML_data_not_open_closed_entropy_difference",
                unit_warning="Published F header says J/mol, but F values are consistent with kJ/mol and S with J/mol/K",
                limitations=["No paired open/closed phase assignment or independent DFT entropy benchmark established here",
                             "The 0b table has substantial imaginary band modes and is not a stable harmonic reference",
                             "No fresh CIF cleaning, ML inference, phonon sampling or convergence was performed"])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if args.out.exists():p.error("output exists")
    result=audit(args.archive)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("x") as f:f.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
