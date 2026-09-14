"""Very fast periodic elastic-network entropy proxy; no MD or electronic DFT.

This is a NEW approximate vibrational model, not a repair of unstable phonons
from another force field. Phenomenological transverse springs represent angular
resistance. Its entropy has no established MOF relative-free-energy accuracy.
"""
from itertools import product
import math
import numpy as np

R = .008314462618  # kJ mol^-1 K^-1
EV = 96.48533212331002
H_OVER_K_THz = 47.9924307337  # K/THz
EV_AMU_A2_TO_THZ = 15.633304239856194
DEFAULT = dict(covalent_k=15., hydrogen_k=40., metal_k=3., transverse_fraction=.03,
               packing_k=.03, packing_cutoff_A=4.5)


def build_edges(atoms, parameters=None):
    from ase.data import covalent_radii
    from ase.neighborlist import neighbor_list
    if set(parameters or {})-set(DEFAULT):raise ValueError("unknown elastic-network parameter")
    p=dict(DEFAULT,**(parameters or {}))
    if any(not math.isfinite(float(v)) or v<0 for v in p.values()) or min(p["covalent_k"],p["hydrogen_k"],p["metal_k"],p["transverse_fraction"],p["packing_cutoff_A"])<=0:
        raise ValueError("positive bond springs/transverse fraction/cutoff; nonnegative packing spring required")
    symbols=atoms.get_chemical_symbols();z=atoms.numbers
    metals={"Al","Zn","Cu","Co","Ni","Fe","Cr","Zr","Hf","Ti","Mg","Ca","Sc","V","Mn"}
    donors={"O","N","S","F","Cl"}
    ii,jj,ss,vec,dist=neighbor_list("ijSDd",atoms,p["packing_cutoff_A"])
    bonds=[];nonbonded=[];adj=[[] for _ in atoms]
    for i,j,s,v,d in zip(ii,jj,ss,vec,dist):
        i,j=int(i),int(j);s=tuple(int(x) for x in s)
        a,b=symbols[i],symbols[j]
        metal=(a in metals and b in donors) or (b in metals and a in donors)
        if "H" in (a,b) and d<1.30 and a!=b:k=p["hydrogen_k"];kind="X-H"
        elif metal and d<covalent_radii[z[i]]+covalent_radii[z[j]]+.35:k=p["metal_k"];kind="metal-donor"
        elif a not in metals and b not in metals and "H" not in (a,b) and d<covalent_radii[z[i]]+covalent_radii[z[j]]+.25:k=p["covalent_k"];kind="covalent"
        else:nonbonded.append((i,j,s,v,d));continue
        if d<.75:raise ValueError("short-contact structure is not suitable for proxy entropy")
        adj[i].append((j,s));bonds.append((i,j,s,v/d,k,kind))
    if any(not neighbors for neighbors in adj):raise ValueError("unconnected atom in bonded network")
    # Exclude true periodic 1-2 and 1-3 pairs from packing springs.
    excluded=set()
    for i,neighbors in enumerate(adj):
        for j,s in neighbors:
            excluded.add((i,j,s))
            for k,t in adj[j]:excluded.add((i,k,tuple(a+b for a,b in zip(s,t))))
    edges=[]
    for i,j,s,n,k,kind in bonds:
        block=k*(np.outer(n,n)+p["transverse_fraction"]*(np.eye(3)-np.outer(n,n)))
        edges.append((i,j,s,block,kind))
    for i,j,s,v,d in nonbonded:
        if p["packing_k"]==0 or "H" in (symbols[i],symbols[j]) or (i,j,s) in excluded:continue
        switch=.5*(1+math.cos(math.pi*d/p["packing_cutoff_A"]))
        k=p["packing_k"]*math.exp(-((d-3.)/.8)**2)*switch
        edges.append((i,j,s,k*np.outer(v/d,v/d),"packing"))
    return edges,p


def network_thermodynamics(atoms, formula_units, temperature_K, *, mesh=(3,3,3), parameters=None):
    if not math.isfinite(temperature_K) or temperature_K<=0 or not math.isfinite(formula_units) or formula_units<=0:
        raise ValueError("positive finite T and formula_units required")
    if len(mesh)!=3 or any(int(n)!=n or n<1 for n in mesh):raise ValueError("positive integer q mesh required")
    edges,p=build_edges(atoms,parameters)
    n=len(atoms);mass=atoms.get_masses();sqrtmass=np.sqrt(mass)
    qpoints=list(product(*[np.arange(m)/m for m in mesh]));spectra=[]
    for q in qpoints:
        d=np.zeros((3*n,3*n),dtype=complex)
        for i,j,s,k,kind in edges:
            x=slice(3*i,3*i+3);y=slice(3*j,3*j+3)
            d[x,x]+=k/mass[i]
            d[x,y]-=k*np.exp(2j*np.pi*np.dot(q,s))/(sqrtmass[i]*sqrtmass[j])
        if not np.allclose(d,d.conj().T,atol=1e-8):raise ValueError("periodic network is not Hermitian")
        eig=np.linalg.eigvalsh(d)
        if eig.min() < -1e-7:raise ValueError("positive-spring network unexpectedly indefinite")
        zero=eig<1e-9
        expected=3 if np.allclose(q,0) else 0
        if int(zero.sum())!=expected:raise ValueError("disconnected/floppy network beyond the three Gamma translations")
        # Numerical zero removal only for known translations, not negative ML phonons.
        freq=np.sqrt(eig[~zero])*EV_AMU_A2_TO_THZ
        spectra.append(freq)
    freq=np.concatenate(spectra)
    x=H_OVER_K_THz*freq/temperature_K
    thermal_population=np.exp(-x)/(-np.expm1(-x))
    norm=len(qpoints)*formula_units
    entropy=R*1000*np.sum(x*thermal_population-np.log(-np.expm1(-x)))/norm
    energy=R*temperature_K*np.sum(x*(.5+thermal_population))/norm
    free=energy-temperature_K*entropy/1000
    return dict(method="phenomenological_periodic_elastic_network",parameters=p,mesh=list(mesh),
                S_proxy_J_mol_fu_K=float(entropy),U_vib_proxy_kJ_mol_fu=float(energy),
                F_vib_proxy_kJ_mol_fu=float(free),ZPE_proxy_kJ_mol_fu=float(R*H_OVER_K_THz*np.sum(freq)*.5/norm),
                smallest_nonzero_frequency_THz=float(freq.min()),largest_frequency_THz=float(freq.max()),
                directed_bond_edges=sum(e[-1]!="packing" for e in edges),directed_packing_edges=sum(e[-1]=="packing" for e in edges),
                excluded_translations=3,q_points=len(qpoints),formula_units=formula_units,
                warning="Surrogate entropy only; positive network modes do NOT establish stability under the ML force field")


def compare_states(atoms, energies, formula_units, temperature_K, *, mesh=(3,3,3), parameter_sets=None):
    from collections import Counter
    if len(atoms)!=2 or len(energies)!=2 or len(formula_units)!=2:raise ValueError("closed and open states required in that order")
    comp=[{k:v/z for k,v in Counter(a.get_chemical_symbols()).items()} for a,z in zip(atoms,formula_units)]
    if comp[0]!=comp[1]:raise ValueError("phase compositions per formula unit differ")
    if not np.isfinite(energies).all():raise ValueError("finite static energies in kJ/mol-f.u. required")
    choices=parameter_sets or [dict(name="baseline",transverse_fraction=.03,packing_k=.03),
        dict(name="soft_hinges",transverse_fraction=.01,packing_k=.03),
        dict(name="stiff_hinges",transverse_fraction=.1,packing_k=.03),
        dict(name="no_packing",transverse_fraction=.03,packing_k=0.),
        dict(name="stronger_packing",transverse_fraction=.03,packing_k=.1)]
    rows=[];de=float(energies[1]-energies[0])
    for choice in choices:
        parameters={k:v for k,v in choice.items() if k!="name"}
        r=[network_thermodynamics(a,z,temperature_K,mesh=mesh,parameters=parameters) for a,z in zip(atoms,formula_units)]
        ds=r[1]["S_proxy_J_mol_fu_K"]-r[0]["S_proxy_J_mol_fu_K"]
        du=r[1]["U_vib_proxy_kJ_mol_fu"]-r[0]["U_vib_proxy_kJ_mol_fu"]
        df=de+du-temperature_K*ds/1000
        rows.append(dict(name=choice.get("name","custom"),closed=r[0],open=r[1],
                         delta_E_static_kJ_mol_fu=de,delta_S_proxy_J_mol_fu_K=ds,
                         delta_U_vib_proxy_kJ_mol_fu=du,minus_T_delta_S_proxy_kJ_mol_fu=-temperature_K*ds/1000,
                         delta_F_host_proxy_kJ_mol_fu=df))
    fs=[r["delta_F_host_proxy_kJ_mol_fu"] for r in rows]
    return dict(method="static_energy_plus_elastic_network_entropy_proxy",temperature_K=temperature_K,
                classification="fast_unvalidated_multifidelity_screen_not_converged_host_free_energy",
                parameter_scenarios=rows,delta_F_proxy_scenario_range_kJ_mol_fu=[min(fs),max(fs)],
                sign_consistent_across_scenarios=bool(min(fs)>0 or max(fs)<0),
                warnings=["Static energy and entropy come from DIFFERENT models; this is not a thermodynamically consistent ML-potential free energy",
                          "Scenario range is parameter sensitivity, not a confidence interval or known error bound",
                          "No MD, DFT, quantum electronic optimization, finite-displacement ML phonons or temperature sampling",
                          "A positive proxy spectrum does not establish stability under the ML potential or include unknown configurational/mixing entropy"])
