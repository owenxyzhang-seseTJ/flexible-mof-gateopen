"""Pure stable vapor/supercritical reservoirs. Never relabel ideal loadings."""
import math


def reservoir_point(temperature_K, pressure_bar, *, backend="ideal", fluid="CO2"):
    t, p = float(temperature_K), float(pressure_bar)
    if not all(math.isfinite(x) and x > 0 for x in (t, p)):
        raise ValueError("positive finite temperature and pressure required")
    if backend == "ideal":
        return dict(temperature_K=t, pressure_bar=p, fugacity_bar=p,
                    fugacity_coefficient=1., phase="ideal_assumption", backend=backend, fluid=fluid)
    if backend != "coolprop":
        raise ValueError("backend must be ideal or coolprop")
    import CoolProp
    import CoolProp.CoolProp as CP
    if "&" in fluid or "::" in fluid:
        raise ValueError("only single-component HEOS fluids supported")
    state = CoolProp.AbstractState("HEOS", fluid)
    if t <= state.Ttriple():
        raise ValueError("below-triple-point gas/solid equilibrium is not implemented")
    # No imposed gas branch: reject liquid and two-phase reservoirs.
    psat = None
    if state.Ttriple() < t < state.T_critical():
        state.update(CP.QT_INPUTS, 1., t)
        psat = state.p() / 1e5
        if p > psat * (1 + 1e-10):
            raise ValueError(f"pressure exceeds stable-vapor Psat={psat:g} bar")
    saturated_vapor = psat is not None and abs(p / psat - 1) <= 1e-10
    if saturated_vapor:
        state.update(CP.QT_INPUTS, 1., t)  # saturated vapor endpoint
        vapor_density = state.rhomolar()
        # HEOS labels Q=1 as two-phase and disables fugacity there. Evaluate
        # the SAME saturated-vapor density on its single-phase EOS branch.
        state.specify_phase(CP.iphase_gas)
        state.update(CP.DmolarT_INPUTS, vapor_density, t)
    else:
        state.update(CP.PT_INPUTS, p * 1e5, t)
    allowed = {CP.iphase_gas, CP.iphase_supercritical, CP.iphase_supercritical_gas}
    if state.phase() not in allowed and not saturated_vapor:
        raise ValueError("reservoir is not stable vapor/supercritical gas; liquid/two-phase unsupported")
    phi = float(state.fugacity_coefficient(0))
    if not math.isfinite(phi) or phi <= 0:
        raise ValueError("invalid EOS fugacity")
    return dict(temperature_K=t, pressure_bar=p, fugacity_bar=phi*p,
                fugacity_coefficient=phi, phase=int(state.phase()), backend="CoolProp/HEOS",
                version=CoolProp.__version__, fluid=fluid, psat_bar=psat)
