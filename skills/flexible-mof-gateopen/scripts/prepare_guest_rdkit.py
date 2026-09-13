#!/usr/bin/env python3
"""Seeded gas-phase guest conformers; geometry only, not a RASPA force field."""
import argparse
import json
import math
from pathlib import Path


def prepare(smiles, output, conformers=20, seed=20260913):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("provide valid single-component SMILES")
    if Chem.GetFormalCharge(mol) != 0:
        raise ValueError("this vapor-guest helper requires a neutral molecule")
    if conformers < 1:
        raise ValueError("positive conformer count required")
    output = Path(output)
    if output.exists():
        raise ValueError("use a new output directory")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=conformers, params=params))
    if not ids:
        raise ValueError("conformer embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        method = "MMFF94s"
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=2000, mmffVariant=method)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        method = "UFF"
        results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=2000)
    else:
        raise ValueError("MMFF and UFF parameters unavailable")
    eligible = [(energy, cid) for cid, (status, energy) in zip(ids, results) if status == 0 and math.isfinite(energy)]
    if not eligible:
        raise ValueError("no converged conformer")
    energy, selected = min(eligible)
    output.mkdir(parents=True)
    with Chem.SDWriter(str(output / "guest.sdf")) as writer:
        writer.write(mol, confId=selected)
    Chem.MolToXYZFile(mol, str(output / "guest.xyz"), confId=selected)
    report = {"smiles": smiles, "rdkit_version": rdBase.rdkitVersion, "seed": seed,
              "force_field": method, "selected_conformer": selected, "energy_kcal_mol": energy,
              "conformers": [{"id": cid, "status": int(s), "energy_kcal_mol": e} for cid, (s, e) in zip(ids, results)],
              "limitations": "Gas-phase geometry only. Supply validated RASPA atom types, charges, LJ and intramolecular parameters separately; lowest sampled energy is not proof of global minimum."}
    (output / "optimization.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--smiles", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--conformers", type=int, default=20)
    p.add_argument("--seed", type=int, default=20260913)
    a = p.parse_args()
    print(json.dumps(prepare(a.smiles, a.output_dir, a.conformers, a.seed), indent=2))
