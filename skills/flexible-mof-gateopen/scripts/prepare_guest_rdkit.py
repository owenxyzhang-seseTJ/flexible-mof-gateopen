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


def refine_tblite(output, fmax=.005, steps=300):
    """GFN2-xTB optimization of the selected RDKit conformer; optional dependency."""
    import numpy as np
    from ase import Atoms
    from ase.io import write
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    from rdkit import Chem
    from importlib.metadata import version
    output=Path(output)
    mol=next(iter(Chem.SDMolSupplier(str(output/"guest.sdf"),removeHs=False)))
    atoms=Atoms([a.GetSymbol() for a in mol.GetAtoms()],positions=mol.GetConformer().GetPositions())
    atoms.calc=TBLite(method="GFN2-xTB",verbosity=0)
    opt=BFGS(atoms,logfile=str(output/"gfn2_xtb.log"),trajectory=str(output/"gfn2_xtb.traj"),maxstep=.05)
    converged=bool(opt.run(fmax=fmax,steps=steps))
    record=dict(method="GFN2-xTB",tblite_version=version("tblite"),converged=converged,steps=opt.nsteps,
                energy_eV=float(atoms.get_potential_energy()),max_force_eV_A=float(np.linalg.norm(atoms.get_forces(),axis=1).max()),
                source="Bannwarth et al. DOI 10.1021/acs.jctc.8b01176",periodic=False)
    (output/"gfn2_xtb.json").write_text(json.dumps(record,indent=2)+"\n")
    if not converged:raise ValueError("GFN2-xTB guest optimization did not converge")
    conf=mol.GetConformer()
    for i,pos in enumerate(atoms.positions):conf.SetAtomPosition(i,pos.tolist())
    with Chem.SDWriter(str(output/"guest_refined.sdf")) as writer:writer.write(mol)
    write(str(output/"guest_refined.xyz"),atoms)
    return output/"guest_refined.sdf"


def validate_geometry(mol):
    """Chemical geometry gate independent of the optimizer's convergence flag."""
    import numpy as np
    from rdkit import Chem
    from ase.data import covalent_radii,atomic_numbers
    xyz=mol.GetConformer().GetPositions()
    distances=[]
    for bond in mol.GetBonds():
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        a,b=mol.GetAtomWithIdx(i).GetSymbol(),mol.GetAtomWithIdx(j).GetSymbol()
        length=float(np.linalg.norm(xyz[i]-xyz[j]));order=bond.GetBondTypeAsDouble()
        if {a,b}=={"C","O"} and order==2:lo,hi=1.10,1.32
        else:
            radius=covalent_radii[atomic_numbers[a]]+covalent_radii[atomic_numbers[b]]
            lo,hi=.65*radius,1.20*radius
        if not lo<=length<=hi:raise ValueError(f"implausible optimized bond {i}-{j} {a}-{b}: {length:.4f} A")
        distances.append(dict(atoms=[i,j],elements=[a,b],order=order,length_A=length))
    if Chem.MolToSmiles(Chem.RemoveHs(mol))=="O=C=O":
        if any(not 1.10<=d["length_A"]<=1.22 for d in distances):
            raise ValueError("CO2 C=O geometry outside 1.10-1.22 A; MMFF/UFF convergence alone is insufficient")
        from rdkit.Chem import rdMolTransforms
        oxygen=[a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()=="O"]
        carbon=[a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()=="C"][0]
        if rdMolTransforms.GetAngleDeg(mol.GetConformer(),oxygen[0],carbon,oxygen[1])<175:
            raise ValueError("CO2 must be approximately linear")
    return dict(status="passed",bonds=distances,warning="Bond checks do not prove global conformer or force-field accuracy")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--smiles", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--conformers", type=int, default=20)
    p.add_argument("--seed", type=int, default=20260913)
    a = p.parse_args()
    print(json.dumps(prepare(a.smiles, a.output_dir, a.conformers, a.seed), indent=2))
