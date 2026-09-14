"""Fixed-cell, pinned-atom, bounded Cartesian configurational MC thermodynamic integration.

This is a deliberately explicit constrained ensemble, not a full unconstrained
Frenkel-Ladd absolute crystal free energy. No ideal kinetic/COM correction is implied.
"""
import math
import numpy as np

KB = 8.617333262145e-5  # eV/K
EV = 96.48533212331002


def mc_ti(target_energy, dimension, *, temperature_K=300., spring_eV_A2=2.,
          bound_A=.3, lambdas=None, warmup=1000, samples=5000, stride=5,
          step_A=.05, blocks=10, seed=0):
    """Return F_conf with measure product dx/(1 A); scalar energies in eV.

    U_lambda=(1-lambda)*k*x^2/2 + lambda*U_target. Symmetric single-coordinate
    proposals, hard-wall cube, exact truncated Einstein reference partition.
    The callable receives a flat displacement vector; caller fixes anchor/cell.
    """
    lam=np.linspace(0,1,11) if lambdas is None else np.asarray(lambdas,dtype=float)
    if lam.ndim!=1 or len(lam)<2 or not np.isfinite(lam).all() or lam[0]!=0 or lam[-1]!=1 or np.any(np.diff(lam)<=0):
        raise ValueError("lambda grid must strictly span 0 to 1")
    if dimension<1 or int(dimension)!=dimension or samples<20 or blocks<2 or samples<blocks*2 or warmup<0 or stride<1:
        raise ValueError("invalid dimension or sampling controls")
    if not all(math.isfinite(v) and v>0 for v in (temperature_K,spring_eV_A2,bound_A,step_A)):
        raise ValueError("positive finite MC parameters required")
    beta=1/(KB*temperature_K)
    integral=math.sqrt(2*math.pi/(beta*spring_eV_A2))*math.erf(bound_A*math.sqrt(beta*spring_eV_A2/2))
    reference=-dimension/beta*math.log(integral)
    rng=np.random.default_rng(seed)
    x=np.zeros(dimension)
    def evaluate(x):
        e=float(target_energy(x))
        if not math.isfinite(e): raise ValueError("nonfinite target energy")
        return .5*spring_eV_A2*float(x@x),e
    u0,u1=evaluate(x)
    records=[]
    for l in lam:
        values=[]; accepted=0; production_accepted=0
        for iteration in range(warmup+samples*stride):
            proposal=x.copy(); j=rng.integers(dimension)
            proposal[j]+=rng.uniform(-step_A,step_A)
            if abs(proposal[j]) <= bound_A:
                v0,v1=evaluate(proposal)
                de=(1-l)*(v0-u0)+l*(v1-u1)
                if de<=0 or rng.random()<math.exp(-beta*de):
                    x=proposal;u0,u1=v0,v1;accepted+=1
                    if iteration>=warmup: production_accepted+=1
            if iteration>=warmup and (iteration-warmup)%stride==0:
                values.append(u1-u0)
        means=np.array([np.mean(b) for b in np.array_split(values,blocks)])
        records.append(dict(lambda_value=float(l),mean_dU_eV=float(np.mean(values)),
                            block_SE_eV=float(means.std(ddof=1)/math.sqrt(blocks)),
                            acceptance=production_accepted/(samples*stride),
                            samples=len(values)))
    means=np.array([r["mean_dU_eV"] for r in records])
    weights=np.zeros(len(lam)); weights[:-1]+=np.diff(lam)/2; weights[1:]+=np.diff(lam)/2
    delta=float(weights@means)
    error=float(np.linalg.norm(weights*np.array([r["block_SE_eV"] for r in records])))
    return dict(F_reference_eV=reference,delta_F_TI_eV=delta,F_conf_eV=reference+delta,
                approximate_block_SE_eV=error,windows=records,seed=seed,dimension=dimension,
                bound_A=bound_A,spring_eV_A2=spring_eV_A2,temperature_K=temperature_K,
                measure="product dx/(1 Angstrom), pinned atom, fixed cell, labelled atoms, hard-wall Cartesian cube",
                classification="constrained_configurational_TI_not_unconstrained_host_F",
                warnings=["Block SE assumes decorrelation, including between lambda windows; verify block-size plateaus and independent seeds.",
                          "Converge lambda quadrature, warmup, proposals and bound; no COM/kinetic/permutation correction.",
                          "Hard-wall basin depends on bound; it is not an automatically identified physical phase basin."])


def framework_mc_ti(atoms, formula_units, **kwargs):
    if len(atoms)<2 or not math.isfinite(formula_units) or formula_units<=0:
        raise ValueError("at least two atoms and positive finite formula units required")
    origin=atoms.positions.copy()
    def target(x):
        positions=origin.copy();positions[1:]+=x.reshape((-1,3))
        atoms.set_positions(positions)
        return atoms.get_potential_energy()
    try:
        result=mc_ti(target,3*(len(atoms)-1),**kwargs)
    finally:
        atoms.set_positions(origin)
    return dict(result,F_conf_kJ_mol_fu=result["F_conf_eV"]*EV/formula_units,
                formula_units=formula_units,anchor_atom_index=0)
