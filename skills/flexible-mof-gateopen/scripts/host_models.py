"""Optional ASE ML calculators and finite-displacement harmonic free energies."""
from collections import Counter
import hashlib
from importlib.metadata import version
import math
from pathlib import Path
import numpy as np

EV = 96.48533212331002


def calculator(backend="chgnet", model=None, device="cpu", dispersion=False, head=None):
    """Load published checkpoints; no user-supplied executable plugins."""
    if backend == "chgnet":
        if head is not None: raise ValueError("head applies only to MACE")
        if dispersion:
            raise ValueError("CHGNet dispersion correction not configured; no silent double counting")
        from chgnet.model.model import CHGNet
        from chgnet.model.dynamics import CHGNetCalculator
        net = CHGNet.load(model_name=model or "0.3.0", use_device=device)
        calc = CHGNetCalculator(model=net, use_device=device, on_isolated_atoms="error")
        package="chgnet"
    elif backend == "mace":
        from mace.calculators import mace_mp
        if head == "pbe_d3" and dispersion:
            raise ValueError("pbe_d3 head already includes dispersion; do not double count D3")
        calc = mace_mp(model=model or "small", device=device, default_dtype="float64", dispersion=dispersion,
                       **({"head": head} if head else {}))
        package="mace-torch"
    else:
        raise ValueError("ML backend must be chgnet or mace")
    metadata=dict(backend=backend,model=model or ("0.3.0" if backend=="chgnet" else "small"),
                  package_version=version(package),device=device,dispersion=dispersion,head=head)
    # Hash actual loaded tensor state, including named parameters and buffers.
    nets = [net] if backend == "chgnet" else getattr(calc,"models",[])
    if not nets and hasattr(calc,"calcs"):
        nets=getattr(calc.calcs[0],"models",[])
    digest=hashlib.sha256()
    for network in nets:
        for name,tensor in sorted(network.state_dict().items()):
            digest.update(name.encode()); digest.update(str(tensor.shape).encode())
            digest.update(tensor.detach().cpu().contiguous().numpy().tobytes())
    metadata["weights_sha256"]=digest.hexdigest() if nets else None
    if model and Path(model).is_file():
        metadata["checkpoint_sha256"]=hashlib.sha256(Path(model).read_bytes()).hexdigest()
    return calc,metadata


def energy_record(atoms, formula_units):
    z=float(formula_units)
    if not math.isfinite(z) or z<=0:
        raise ValueError("positive formula units per input cell required")
    energy=float(atoms.get_potential_energy())
    forces=np.asarray(atoms.get_forces())
    stress=np.asarray(atoms.get_stress())
    if not np.isfinite([energy,*forces.ravel(),*stress.ravel()]).all():
        raise ValueError("ML returned nonfinite energy, forces or stress")
    return dict(energy_eV_cell=energy,energy_kJ_mol_fu=energy*EV/z,
                formula_units=z,composition_per_fu={k:v/z for k,v in Counter(atoms.get_chemical_symbols()).items()},
                volume_A3_fu=float(atoms.get_volume()/z),
                max_force_eV_A=float(np.linalg.norm(forces,axis=1).max()),
                stress_eV_A3=stress.tolist(),atom_count=len(atoms))


def energy_difference(records):
    ref=records[0]
    for r in records:
        a,b=r["composition_per_fu"],ref["composition_per_fu"]
        if set(a)!=set(b) or any(not math.isclose(a[k],b[k],abs_tol=1e-8) for k in a):
            raise ValueError("phase compositions per formula unit differ; compare same activation state")
    return [dict(r,delta_E_kJ_mol_fu=r["energy_kJ_mol_fu"]-ref["energy_kJ_mol_fu"]) for r in records]


def harmonic_state_differences(records, *, max_negative_frequency_THz=1e-4):
    """Combine static energy and vibrational thermodynamics from ONE calculator run.

    Reference is the first state (put closed first). This is a harmonic
    estimate, not a correction for configurational or strong anharmonic entropy.
    """
    if len(records)<2: return []
    names=[r["name"] for r in records]
    if len(set(names))!=len(names): raise ValueError("state names must be unique")
    if not math.isfinite(max_negative_frequency_THz) or max_negative_frequency_THz<0:
        raise ValueError("invalid negative-frequency threshold")
    energy_difference([r["base"] for r in records])  # same composition per f.u.
    ref=records[0]
    def valid(r):
        return bool(r.get("thermal")) and r["min_frequency_THz"]>= -max_negative_frequency_THz
    results=[]
    for state in records[1:]:
        if not valid(ref) or not valid(state):
            results.append(dict(state=state["name"],reference=ref["name"],status="blocked_unstable_phonons",
                                thermal=None,negative_frequency_threshold_THz=max_negative_frequency_THz))
            continue
        rt={r["temperature_K"]:r for r in ref["thermal"]}
        st={r["temperature_K"]:r for r in state["thermal"]}
        if set(rt)!=set(st) or len(rt)!=len(ref["thermal"]) or len(st)!=len(state["thermal"]):
            raise ValueError("unique matching temperature grids required")
        de=state["base"]["energy_kJ_mol_fu"]-ref["base"]["energy_kJ_mol_fu"]
        thermal=[]
        for t in sorted(rt):
            dfv=st[t]["F_vib_kJ_mol_fu"]-rt[t]["F_vib_kJ_mol_fu"]
            ds=st[t]["S_J_mol_fu_K"]-rt[t]["S_J_mol_fu_K"]
            minus_ts=-t*ds/1000
            duv=dfv-minus_ts  # U_vib includes ZPE and thermal excitation
            df=de+dfv
            total_difference=st[t]["F_total_kJ_mol_fu"]-rt[t]["F_total_kJ_mol_fu"]
            if not np.isfinite([de,dfv,ds,df,duv]).all() or not math.isclose(df,total_difference,abs_tol=1e-7):
                raise ValueError("inconsistent static/vibrational free-energy bookkeeping")
            thermal.append(dict(temperature_K=t,delta_E_static_kJ_mol_fu=de,
                                delta_U_vib_kJ_mol_fu=duv,delta_S_vib_J_mol_fu_K=ds,
                                minus_T_delta_S_kJ_mol_fu=minus_ts,delta_F_vib_kJ_mol_fu=dfv,
                                delta_F_host_harmonic_kJ_mol_fu=df))
        results.append(dict(state=state["name"],reference=ref["name"],status="harmonic_estimate",
                            thermal=thermal,negative_frequency_threshold_THz=max_negative_frequency_THz,
                            warning="Same calculator; each state uses its own volume. No PV, anharmonic or configurational correction"))
    return results


def phonon_free_energy(atoms, formula_units, temperatures, supercell=(2,2,2), mesh=(8,8,8),
                       displacement=.01, imaginary_tolerance_THz=.1, force_tolerance=.03):
    """Phonopy primitive=identity: thermal output normalized to the INPUT cell."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from ase import Atoms
    if any(int(x)!=x or x<1 for x in (*supercell,*mesh)) or len(supercell)!=3 or len(mesh)!=3:
        raise ValueError("positive integer 3-vectors required")
    if displacement<=0 or imaginary_tolerance_THz<0 or force_tolerance<=0:
        raise ValueError("invalid displacement/tolerance")
    ts=np.asarray(temperatures,dtype=float)
    if ts.ndim!=1 or not len(ts) or np.any(ts<0) or not np.isfinite(ts).all():
        raise ValueError("finite nonnegative temperatures required")
    base=energy_record(atoms,formula_units)
    if base["max_force_eV_A"]>force_tolerance:
        raise ValueError("structure is not stationary; optimize positions at fixed volume and re-audit first")
    unit=PhonopyAtoms(symbols=atoms.get_chemical_symbols(),cell=atoms.cell.array,
                      scaled_positions=atoms.get_scaled_positions())
    ph=Phonopy(unit,supercell_matrix=np.diag(supercell),primitive_matrix=np.eye(3))
    ph.generate_displacements(distance=displacement)
    forces=[]
    drift=[]
    for cell in ph.supercells_with_displacements:
        a=Atoms(symbols=cell.symbols,cell=cell.cell,scaled_positions=cell.scaled_positions,pbc=True)
        a.calc=atoms.calc
        f=a.get_forces()
        if not np.isfinite(f).all(): raise ValueError("nonfinite displacement forces")
        drift.append(float(np.linalg.norm(f.mean(axis=0))))
        forces.append(f-f.mean(axis=0))
    ph.forces=forces
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.run_mesh(mesh,is_gamma_center=True)
    frequencies=np.asarray(ph.mesh.frequencies)
    minimum=float(frequencies.min())
    if not np.isfinite(frequencies).all(): raise ValueError("nonfinite phonon frequencies")
    result=dict(base=base,supercell=list(supercell),mesh=list(mesh),displacement_A=displacement,
                min_frequency_THz=minimum,imaginary_tolerance_THz=imaginary_tolerance_THz,
                nonpositive_frequency_entries=int(np.count_nonzero(frequencies<=1e-5)),
                soft_negative_entries=int(np.count_nonzero(frequencies< -1e-5)),
                displacements=len(forces),max_force_drift_eV_A=max(drift),
                mesh_qpoints=ph.mesh.qpoints.tolist(),mesh_frequencies_THz=frequencies.tolist(),
                phonopy_version=version("phonopy"),normalization="input cell / formula_units; primitive matrix identity")
    if minimum < -imaginary_tolerance_THz:
        return dict(result,status="unstable_phonons",thermal=None,
                    warning="No harmonic F supplied: genuine imaginary modes require resolving instability/anharmonic sampling")
    ph.run_thermal_properties(temperatures=ts,cutoff_frequency=1e-5,pretend_real=False)
    thermal=ph.thermal_properties.thermal_properties
    # Public property is (temperatures, free energy, entropy, Cv), kJ/mol of primitive cells.
    temps,free,entropy,cv=thermal
    return dict(result,status="harmonic_screen",thermal=[dict(temperature_K=float(t),
                 F_vib_kJ_mol_fu=float(f/formula_units),
                 F_total_kJ_mol_fu=float(f/formula_units+base["energy_kJ_mol_fu"]),
                 S_J_mol_fu_K=float(s/formula_units),Cv_J_mol_fu_K=float(c/formula_units))
                 for t,f,s,c in zip(temps,free,entropy,cv)],
                 warning="ZPE included; <=1e-5 THz modes excluded. Converge soft modes, displacement, mesh and supercell; not MOF validation")
