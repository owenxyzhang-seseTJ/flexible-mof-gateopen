"""Auditable numerical post-processing; supplied surfaces are not predictions."""
import math
from copy import deepcopy
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline, RectBivariateSpline
from scipy.optimize import brentq, minimize_scalar, root

R = .008314462618
PV = 6.02214076e-5  # bar A^3 -> kJ/mol
EV = 96.48533212331002


def vector(value, name, minimum=2):
    x = np.asarray(value, dtype=float)
    if x.ndim != 1 or len(x) < minimum or not np.isfinite(x).all():
        raise ValueError(f"{name} must be a finite 1D array (length >= {minimum})")
    return x


def increasing(value, name, positive=False, minimum=2):
    x = vector(value, name, minimum)
    if np.any(np.diff(x) <= 0) or (positive and np.any(x <= 0)):
        raise ValueError(f"{name} must be strictly increasing" + (" and positive" if positive else ""))
    return x


def adsorption_omega(fugacity_bar, loading_per_fu, targets, temperature_K, *, henry_confirmed=False):
    """Piecewise linear N(log f), analytically integrated including Henry tail."""
    f = increasing(fugacity_bar, "fugacity", True)
    n = vector(loading_per_fu, "loading")
    ft = vector(targets, "target fugacity", 1)
    if len(n) != len(f) or np.any(n < 0) or np.any(ft <= 0) or np.any(ft > f[-1]):
        raise ValueError("loading/grid mismatch, negative loading or fugacity extrapolation")
    if not henry_confirmed:
        raise ValueError("confirm the first isotherm point is in the Henry regime")
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("positive temperature required")
    x = np.log(f)
    cumulative = np.r_[n[0], n[0] + np.cumsum((n[1:]+n[:-1])*.5*np.diff(x))]
    result = []
    for t in ft:
        if t <= f[0]:
            integral = n[0]*t/f[0]
        else:
            i = min(np.searchsorted(f, t)-1, len(f)-2)
            dx = math.log(t/f[i])
            slope = (n[i+1]-n[i])/(x[i+1]-x[i])
            integral = cumulative[i]+n[i]*dx+.5*slope*dx*dx
        result.append(-R*temperature_K*integral)
    return np.asarray(result)


def all_crossings(pressure_bar, delta):
    """All isolated roots of a PCHIP-interpolated nodal difference; plateaus explicit."""
    p = increasing(pressure_bar, "pressure", True)
    y = vector(delta, "delta")
    if len(y) != len(p):
        raise ValueError("grid mismatch")
    fun = PchipInterpolator(np.log(p), y, extrapolate=False)
    roots, plateaus = [], []
    for i in range(len(p)-1):
        if y[i] == y[i+1] == 0:
            plateaus.append([float(p[i]), float(p[i+1])])
        elif y[i]*y[i+1] < 0:
            roots.append(float(np.exp(brentq(fun, math.log(p[i]), math.log(p[i+1])))))
    for i in np.flatnonzero(y == 0):
        if (i == 0 or y[i-1] != 0) and (i == len(y)-1 or y[i+1] != 0):
            roots.append(float(p[i]))
    return dict(pressure_bar=sorted(set(roots)), coexistence_intervals_bar=plateaus,
                interpolation="PCHIP of nodal delta vs log pressure; refine grid to converge")


def state_adsorption_omega(state, targets, temperature_K):
    """Explicit zero/linear scenarios; never infer physical exclusion from zeros.

    A Henry scenario is an assumed N=K*f over its declared range, NOT a
    finite-pressure bound inferred from a dilute Widom calculation.
    """
    model = state.get("adsorption_model")
    if model is None:
        return adsorption_omega(state["fugacity_bar"], state["loading_molecules_fu"],
                                targets, temperature_K,
                                henry_confirmed=state.get("henry_confirmed", False))
    kind = model.get("kind")
    if kind not in {"nonadsorbing", "henry_scenario"}:
        raise ValueError("unknown adsorption_model kind")
    if model.get("interpretation") != "conditional_scenario" or not str(model.get("evidence", "")).strip():
        raise ValueError("analytical adsorption model requires evidence and conditional_scenario interpretation")
    ft = vector(targets, "target fugacity", 1)
    maximum = float(model["max_fugacity_bar"])
    if not math.isfinite(maximum) or maximum <= 0 or np.any(ft <= 0) or np.any(ft > maximum*(1+1e-12)):
        raise ValueError("analytical adsorption model fugacity range/extrapolation error")
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("positive temperature required")
    k = 0. if kind == "nonadsorbing" else float(model["K_molecules_fu_bar"])
    if not math.isfinite(k) or k < 0:
        raise ValueError("finite nonnegative Henry coefficient required")
    if kind == "henry_scenario" and not model.get("finite_pressure_assumption"):
        raise ValueError("Henry scenario requires a finite_pressure_assumption, not just dilute Widom evidence")
    # Integral_0^f K*f' dln(f') = K*f exactly, independent of pressure grid.
    return -R*temperature_K*k*ft


def phase_diagram(spec, reservoir):
    """Multiple T slices; no T interpolation or metastable branch as global phase."""
    t = float(spec["temperature_K"])
    p = increasing(spec["pressure_bar"], "pressure", True)
    f = np.asarray([r["fugacity_bar"] for r in reservoir])
    if len(f) != len(p):
        raise ValueError("reservoir grid mismatch")
    states = spec["states"]
    names = [s["name"] for s in states]
    if len(states) < 2 or len(set(names)) != len(names):
        raise ValueError("at least two uniquely named states required")
    curves = []
    for s in states:
        host, vol = float(s["host_F_kJ_mol_fu"]), float(s["volume_A3_fu"])
        if not math.isfinite(host) or not math.isfinite(vol) or vol <= 0 or not s.get("host_source"):
            raise ValueError("finite sourced host F and positive volume required")
        # Dataset axis is actual simulation fugacity, not physical pressure.
        ads = state_adsorption_omega(s, f, t)
        curves.append(host + PV*p*vol + ads)
    curves = np.array(curves)
    interpolants = [PchipInterpolator(np.log(p), curve) for curve in curves]
    crossings = []
    for i in range(len(states)):
        for j in range(i+1, len(states)):
            # Roots of the same interpolated potentials used to classify stability.
            # Difference of piecewise cubics may have >1 root in one interval.
            from scipy.interpolate import PPoly
            difference = PPoly(interpolants[j].c-interpolants[i].c, np.log(p), extrapolate=False)
            roots = difference.roots(extrapolate=False)
            plateaus = [[float(p[k]),float(p[k+1])] for k in range(len(p)-1)
                        if np.all(difference.c[:,k] == 0)]
            for x in sorted(set(float(v) for v in roots if np.isfinite(v) and np.log(p[0]) <= v <= np.log(p[-1]))):
                energies = np.array([fun(x) for fun in interpolants])
                eps = min(1e-5, (np.log(p[-1])-np.log(p[0]))*1e-5)
                left = max(np.log(p[0]), x-eps)
                right = min(np.log(p[-1]), x+eps)
                crossings.append(dict(pair=[names[i], names[j]], pressure_bar=float(np.exp(x)),
                                      globally_stable=bool(energies[i] <= min(energies)+1e-7),
                                      stable_before=names[int(np.argmin([fun(left) for fun in interpolants]))],
                                      stable_after=names[int(np.argmin([fun(right) for fun in interpolants]))]))
            if plateaus:
                crossings.append(dict(pair=[names[i],names[j]], coexistence_intervals_bar=plateaus))
    return dict(temperature_K=t, pressure_bar=p.tolist(), reservoir=reservoir,
                state_names=names, omega_kJ_mol_fu=curves.tolist(),
                stable_state=[names[i] for i in np.argmin(curves, axis=0)], crossings=crossings,
                classification="conditional_equilibrium_phase_diagram",
                mechanics="P_external_equals_reservoir; P*V included once",
                analytical_adsorption_scenarios={s["name"]:s["adsorption_model"] for s in states if s.get("adsorption_model")},
                warnings=["Grid convergence required; host F must exclude PV.",
                          "Explicit analytical adsorption models are conditional scenarios, not measured zeros or certified bounds.",
                          "Roots are equilibrium, not kinetic opening/closing pressures."])


def closed_adsorption_sensitivity(spec, reservoir, coefficients, closed_name="closed"):
    """Compare zero adsorption with explicit residual linear-uptake scenarios."""
    values=vector(coefficients,"Henry sensitivity coefficients",1)
    if np.any(values<0):raise ValueError("nonnegative Henry sensitivity coefficients required")
    base=deepcopy(spec)
    closed=next((s for s in base["states"] if s["name"]==closed_name),None)
    if closed is None or closed.get("adsorption_model",{}).get("kind")!="nonadsorbing":
        raise ValueError("explicit nonadsorbing reference state required")
    results=[]
    for k in sorted(set([0.]+values.tolist())):
        model=closed["adsorption_model"]
        model.update(kind="nonadsorbing" if k==0 else "henry_scenario",K_molecules_fu_bar=k,
                     finite_pressure_assumption="Assume N_closed=K*f throughout the declared range for sensitivity only; no claim of a Widom-derived bound")
        diagram=phase_diagram(base,reservoir)
        results.append(dict(K_molecules_fu_bar=k,
            maximum_closed_loading_scenario_molecules_fu=k*max(r["fugacity_bar"] for r in reservoir),
            maximum_closed_adsorption_stabilization_kJ_mol_fu=R*spec["temperature_K"]*k*max(r["fugacity_bar"] for r in reservoir),
            crossings=diagram["crossings"],
            equilibrium_opening_pressures_bar=[r["pressure_bar"] for r in diagram["crossings"]
                if r.get("globally_stable") and r.get("stable_before")==closed_name and r.get("stable_after")!=closed_name]))
    return dict(classification="conditional_residual_adsorption_sensitivity_not_confidence_interval",
                scenarios=results,warning="Widom cannot rule out cooperative uptake or unsampled cavities; analytical linear scenarios are not rigorous bounds")


def qha_minimize(volume_A3_fu, helmholtz_kJ_mol_fu, pressure_bar=0):
    """Basin-specific shape-preserving F(V)+PV minimization (not an EOS fit)."""
    v = increasing(volume_A3_fu, "volume", True, minimum=5)
    f = vector(helmholtz_kJ_mol_fu, "Helmholtz F", 5)
    if len(v) != len(f) or not math.isfinite(pressure_bar) or pressure_bar < 0:
        raise ValueError("invalid QHA data/pressure")
    fun = PchipInterpolator(v, f)
    candidates = [(float(x),float(fun(x)+PV*pressure_bar*x)) for x in v]
    for a,b in zip(v[:-1],v[1:]):
        opt = minimize_scalar(lambda x: fun(x)+PV*pressure_bar*x, bounds=(a,b), method="bounded")
        candidates.append((float(opt.x),float(opt.fun)))
    x,y = min(candidates,key=lambda pair: pair[1])
    return dict(volume_A3_fu=x, G_kJ_mol_fu=y,
                boundary_minimum=bool(np.isclose(x,v[0]) or np.isclose(x,v[-1])),
                method="basin-specific PCHIP QHA minimization; no extrapolation")


def landscape_extrema(q, energy):
    """Stationary points of a supplied 1D cubic-spline constrained free energy."""
    q = increasing(q,"q",minimum=5)
    y = vector(energy,"landscape",5)
    if len(q) != len(y):
        raise ValueError("landscape length mismatch")
    fun = CubicSpline(q,y,extrapolate=False)
    roots = fun.derivative().roots(extrapolate=False)
    extrema=[]
    for x in roots:
        if not np.isfinite(x) or x <= q[0]+1e-9 or x >= q[-1]-1e-9:
            continue
        curvature=float(fun(x,2))
        extrema.append(dict(q=float(x),F_kJ_mol_fu=float(fun(x)),curvature=curvature,
                            kind="minimum" if curvature>1e-8 else "maximum" if curvature< -1e-8 else "marginal"))
    barriers=[]
    for k, e in enumerate(extrema):
        if e["kind"] != "maximum":
            continue
        for direction in (-1,1):
            i=k+direction
            if 0 <= i < len(extrema) and extrema[i]["kind"]=="minimum":
                barriers.append(dict(from_q=extrema[i]["q"],saddle_q=e["q"],
                                     barrier_kJ_mol_fu=e["F_kJ_mol_fu"]-extrema[i]["F_kJ_mol_fu"]))
    boundary = []
    if fun(q[0],1) > 0: boundary.append(float(q[0]))
    if fun(q[-1],1) < 0: boundary.append(float(q[-1]))
    return dict(extrema=extrema,barriers=barriers,boundary_minima_q=boundary,
                warning="1D constrained path only; spline and q-grid convergence required")


def landscape_scan(q, pressure_bar, surfaces):
    p=increasing(pressure_bar,"pressure",True)
    if len(surfaces)!=len(p): raise ValueError("one landscape per pressure required")
    slices=[dict(pressure_bar=float(x),**landscape_extrema(q,y)) for x,y in zip(p,surfaces)]
    brackets=[]
    for a,b in zip(slices[:-1],slices[1:]):
        na=sum(e["kind"]=="minimum" for e in a["extrema"])
        nb=sum(e["kind"]=="minimum" for e in b["extrema"])
        if na!=nb:
            brackets.append(dict(pressure_bracket_bar=[a["pressure_bar"],b["pressure_bar"]],
                                 minima_before=na,minima_after=nb,
                                 classification="candidate_spinodal_bracket_not_certified",
                                 boundary_warning=bool(a["boundary_minima_q"] or b["boundary_minima_q"])))
    spinodals=[]
    if len(p)>=4:
        q=increasing(q,"q",minimum=5)
        y=np.asarray(surfaces,dtype=float)
        if y.shape != (len(p),len(q)) or not np.isfinite(y).all():
            raise ValueError("landscape must be finite pressure-by-q matrix")
        x=np.log(p)
        spline=RectBivariateSpline(x,q,y,kx=3,ky=3,s=0)
        # Scale equations to energy units; solve dF/dq=d2F/dq2=0.
        span=q[-1]-q[0]
        scale=max(float(np.ptp(y)),1.)
        def equations(z):
            return [float(spline.ev(z[0],z[1],dy=1))*span/scale,
                    float(spline.ev(z[0],z[1],dy=2))*span**2/scale]
        for bracket in brackets:
            lo,hi=np.log(bracket["pressure_bracket_bar"])
            for seed in np.linspace(q[0],q[-1],min(21,len(q)))[1:-1]:
                sol=root(equations,[(lo+hi)/2,seed])
                a,b=sol.x
                if not sol.success or np.linalg.norm(equations(sol.x))>1e-7:
                    continue
                if not (lo-1e-8<=a<=hi+1e-8 and q[0]+1e-6*span<b<q[-1]-1e-6*span):
                    continue
                if any(abs(a-math.log(s["pressure_bar"]))<1e-5 and abs(b-s["q"])<1e-5*span for s in spinodals):
                    continue
                h=span*1e-4
                third=float((spline.ev(a,b+h,dy=2)-spline.ev(a,b-h,dy=2))/(2*h))
                coupling=float(spline.ev(a,b,dx=1,dy=1))
                if abs(third*span**3/scale)<1e-6 or abs(coupling*span/scale)<1e-6:
                    continue  # degenerate critical point, not a generic saddle-node
                spinodals.append(dict(pressure_bar=float(np.exp(a)),q=float(b),
                                      dF_dq=float(spline.ev(a,b,dy=1)),d2F_dq2=float(spline.ev(a,b,dy=2)),
                                      d3F_dq3=third,d2F_dlogP_dq=coupling,
                                      classification="spinodal_of_supplied_1D_spline_not_material_validation"))
    return dict(slices=slices,candidate_spinodals=brackets,spline_spinodals=spinodals,
                warning="Confirm stationary-point annihilation F'=F''=0; count changes alone are not proof.")
