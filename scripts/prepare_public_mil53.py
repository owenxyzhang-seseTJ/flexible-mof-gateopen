#!/usr/bin/env python3
"""Reconstruct missing H in the public H-free MIL-53(Al) pair; NOT a refinement.

Only this known Al4C32O20 translational supercell is supported. Other materials
must use clean-cif and an independently reviewed chemistry record.
"""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys

import gemmi
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
try:
    from flexible_mof_gateopen.preflight import sha256, structure_gate
except ModuleNotFoundError:
    from preflight import sha256, structure_gate


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def load_auditor(path):
    spec = importlib.util.spec_from_file_location("clean_cif_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def audit(path, auditor):
    pm, structure = auditor.parse_pymatgen(path)
    ar, atoms = auditor.parse_ase(path)
    if not pm["success"] or not ar["success"]:
        raise ValueError("independent structure parsers failed")
    return dict(raw=auditor.raw_metadata(path), pymatgen=pm, ase=ar,
                geometry=auditor.structure_geometry(structure, 30),
                spglib=auditor.recovered_symmetry(atoms, .01))


def chemistry(atoms):
    """Material-specific periodic graph and valence checks, not universal cutoffs."""
    sym = atoms.get_chemical_symbols()
    z = sym.count("Al")
    if Counter(sym) != Counter({"Al": z, "C": 8*z, "O": 5*z, "H": 5*z}) or not z:
        raise ValueError("expected activated Al(OH)(BDC), including all H")
    ii, jj, dd = neighbor_list("ijd", atoms, 2.4)
    bonds, adjacency = [], [[] for _ in atoms]
    limits = {frozenset(("C", "C")): (1.15, 1.70),
              frozenset(("C", "O")): (1.12, 1.50),
              frozenset(("Al", "O")): (1.65, 2.30),
              frozenset(("H", "C")): (.85, 1.20),
              frozenset(("H", "O")): (.80, 1.15)}
    nonbonded = []
    for i, j, d in zip(ii, jj, dd):
        key = frozenset((sym[i], sym[j]))
        bonded = key in limits and limits[key][0] <= d <= limits[key][1]
        if bonded:
            adjacency[i].append(j)
            if i < j: bonds.append((int(i), int(j), float(d)))
        elif i < j:
            nonbonded.append((float(d), int(i), int(j)))
            threshold = 1.15 if "H" in key else 1.45
            if d < threshold:
                raise ValueError(f"unexplained short contact {i}:{sym[i]} {j}:{sym[j]} {d}")
    for i, s in enumerate(sym):
        ns = Counter(sym[j] for j in adjacency[i])
        ok = (s == "Al" and ns == {"O": 6}) or (s == "H" and sum(ns.values()) == 1)
        ok |= s == "C" and ns in ({"C": 2, "H": 1}, {"C": 3}, {"C": 1, "O": 2})
        ok |= s == "O" and ns in ({"C": 1, "Al": 1}, {"Al": 2, "H": 1})
        if not ok: raise ValueError(f"MIL53 valence/coordination failed at {i}:{s}: {ns}")
    # Every vertex must be in a framework component, with no detached solvent.
    visited, stack = set(), [0]
    while stack:
        i = stack.pop()
        if i not in visited:
            visited.add(i); stack.extend(adjacency[i])
    if len(visited) != len(atoms): raise ValueError("disconnected guest/framework fragments")
    al_o = [d for i,j,d in bonds if "Al" in (sym[i],sym[j])]
    return dict(formula_units=z, formula=atoms.get_chemical_formula(),
                Al_O_range_A=[min(al_o), max(al_o)],
                Al_coordination=6, hydrogen_count=sym.count("H"),
                shortest_nonbonded=sorted(nonbonded)[:15],
                bond_indices=sorted([[i,j] for i,j,d in bonds]))


def certify(source, candidate, auditor_path, transformation, extra=None):
    auditor = load_auditor(auditor_path)
    evidence = chemistry(read(str(candidate)))
    sa = candidate.with_suffix(".source_audit.json")
    ca = candidate.with_suffix(".candidate_audit.json")
    dump(sa, audit(source, auditor)); dump(ca, audit(candidate, auditor))
    report = dict(schema_version=1, status="passed", source_sha256=sha256(source),
                  candidate_sha256=sha256(candidate), transformation_record=transformation,
                  source_audit=dict(path=sa.name, sha256=sha256(sa)),
                  candidate_audit=dict(path=ca.name, sha256=sha256(ca)),
                  disorder=dict(status="absent", reason="Input is an ordered computational P1 representative; not an experimental disorder ensemble"),
                  solvent=dict(status="absent", reason="Periodic single framework graph; exact Al(OH)(BDC) formula; no detached solvent. Activation is a model assumption."),
                  chemistry={k: dict(status="passed", evidence=v) for k,v in {
                      "connectivity": "All atoms in one periodic framework graph; carboxyl/aromatic valence checked",
                      "hydrogen_hosts": f"Each of {evidence['hydrogen_count']} H has one C/O host",
                      "metal_coordination": f"Every Al has 6 O: {evidence['Al_O_range_A']} Angstrom",
                      "nonbonded_contacts": "Material-specific bond exclusions; all residual contacts checked against H/heavy thresholds; numeric pairs in chemistry_details",
                      "composition_charge_balance": "Al(III)+OH(-)+BDC(2-)=0; AlC8H5O5 per formula unit"
                  }.items()}, chemistry_details=evidence, provenance=extra or {})
    target = candidate.with_suffix(".cleaning.json")
    dump(target, report)
    structure_gate(candidate, target)
    return target


def reconstruct(source, output, auditor_path):
    output.mkdir(parents=True, exist_ok=False)
    block = gemmi.cif.read_file(str(source)).sole_block()
    # Inspect raw fields before allowing either parser to expand anything.
    sym = list(block.find_values("_atom_site_type_symbol"))
    xyz = np.array([[gemmi.cif.as_number(x) for x in block.find_values(f"_atom_site_fract_{a}")] for a in "xyz"]).T
    occ = [gemmi.cif.as_number(x) for x in block.find_values("_atom_site_occupancy")]
    if occ and not np.allclose(occ, 1): raise ValueError("disorder requires explicit review")
    if Counter(sym) != Counter({"Al":256,"C":2048,"O":1280}):
        raise ValueError("not the supported public 4x4x4 Al4C32O20 model")
    for tag in ("_atom_site_disorder_group", "_atom_site_disorder_assembly"):
        if any(x not in {".","?","0"} for x in block.find_values(tag)):
            raise ValueError("unexpected disorder fields")
    from ase.cell import Cell
    cellpar = [gemmi.cif.as_number(block.find_value("_cell_"+s)) for s in
               ("length_a","length_b","length_c","angle_alpha","angle_beta","angle_gamma")]
    cellpar[:3] = [v/4 for v in cellpar[:3]]
    cell = Cell.fromcellpar(cellpar)
    mask = np.all((xyz >= -1e-8) & (xyz < .25-1e-8), axis=1)
    selected = np.flatnonzero(mask)
    if len(selected)!=56: raise ValueError("unexpected reduced cell")
    coords=(4*xyz[selected])%1
    # Verify every source atom maps to exactly one reduced site, 64 repetitions.
    counts=np.zeros(56, dtype=int)
    for s,x in zip(sym, (4*xyz)%1):
        d=x-coords; d-=np.round(d)
        hits=[i for i in np.flatnonzero(np.linalg.norm(d@cell.array,axis=1)<.005) if sym[selected[i]]==s]
        if len(hits)!=1: raise ValueError("supercell is not an exact chemical repetition")
        counts[hits[0]]+=1
    if not np.all(counts==64): raise ValueError("unequal translation multiplicities")
    a=Atoms([sym[i] for i in selected], cell=cell, scaled_positions=coords, pbc=True)
    reduced=output/"reduced_source.cif"
    write(str(reduced), a)
    symbols=a.get_chemical_symbols()
    hydrogens=[]; mapping=[]
    for i,s in enumerate(symbols):
        vectors=a.get_distances(i, range(len(a)), mic=True, vector=True)
        dist=np.linalg.norm(vectors, axis=1)
        carbon=[j for j,t in enumerate(symbols) if t=="C" and 1.1<dist[j]<1.7]
        if s=="C" and len(carbon)==2 and not any(t=="O" and dist[j]<1.5 for j,t in enumerate(symbols)):
            near=carbon; length=1.08
        elif s=="O" and not carbon:
            near=[j for j,t in enumerate(symbols) if t=="Al" and 1.65<dist[j]<2.3]; length=.98
            if len(near)!=2: raise ValueError("unresolved hydroxyl geometry")
        else: continue
        direction=-np.sum(vectors[near]/dist[near,None], axis=0)
        if np.linalg.norm(direction)<.05: raise ValueError("ambiguous hydrogen direction")
        hydrogens.append(a.positions[i]+length*direction/np.linalg.norm(direction)); mapping.append(i)
    if len(hydrogens)!=20: raise ValueError("expected 16 aromatic H and 4 hydroxyl H")
    a.extend(Atoms("H20", positions=hydrogens)); a.wrap()
    candidate=output/"host.cif"; write(str(candidate), a)
    record=certify(reduced, candidate, auditor_path,
        "Verified 64-fold translations; reduced 4x4x4 cell; added 16 aromatic and 4 bridging-OH H along local bond bisectors. Hypothetical activated model; no non-H deletions beyond repetitions.",
        dict(public_source_sha256=sha256(source),public_source=str(source.resolve()),
             selected_raw_indices=selected.tolist(),hydrogen_parent_indices=mapping,
             original_atom_count=3584, original_H_count=0))
    return candidate,record


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--clean-cif-auditor",type=Path,required=True)
    args=p.parse_args()
    print(reconstruct(args.source,args.output_dir,args.clean_cif_auditor))
