"""Descriptive fits; a single-T isotherm does not identify a microscopic force field."""
import numpy as np
from scipy.optimize import least_squares


def fit_isotherm(fugacity, loading, energy=None):
    f,n=np.asarray(fugacity,dtype=float),np.asarray(loading,dtype=float)
    if f.ndim!=1 or len(f)<3 or f.shape!=n.shape or not np.isfinite([f,n]).all() or np.any(f<=0) or np.any(n<0) or np.any(np.diff(f)<=0):
        raise ValueError("at least three increasing positive fugacities and finite nonnegative loadings required")
    if not np.any(n): return dict(status="unidentifiable_zero_loading")
    def model(x):
        ns,b=np.exp(x); return ns*b*f/(1+b*f)
    opt=least_squares(lambda x:model(x)-n,np.log([max(n)*1.1,1/np.median(f)]),bounds=([-20,-30],[20,30]))
    ns,b=np.exp(opt.x); pred=model(opt.x)
    result=dict(status="descriptive_fit",model="single_site_Langmuir",n_sat_molecules_fu=float(ns),
                b_per_bar=float(b),rmse_molecules_fu=float(np.sqrt(np.mean((pred-n)**2))),
                predicted_loading=pred.tolist(),
                capacity_identifiable=bool(b*f[-1]>3),
                warning="Single T b is an affinity, not a unique interaction energy; pressure prediction integrates raw N(f), not this fit")
    if energy is not None:
        u=np.asarray(energy,dtype=float)
        if u.shape!=n.shape or not np.isfinite(u).all(): raise ValueError("energy shape/finite check failed")
        # Extensive host-guest U tends to zero at zero loading. Fit U=a*N+b*N².
        design=np.column_stack([n,n*n])
        coeff,_,rank,_=np.linalg.lstsq(design,u,rcond=None)
        result["interaction_energy_fit"]=dict(model="U_host_guest=a*N+b*N^2; extensive energy per framework f.u.",
            a_kJ_mol_adsorbate=float(coeff[0]),b=float(coeff[1]),rank=int(rank),
            rmse_kJ_mol_fu=float(np.sqrt(np.mean((design@coeff-u)**2))),
            warning="Descriptive ensemble energy vs loading; not LJ parameter calibration, adsorption free energy, or isosteric heat")
    return result
