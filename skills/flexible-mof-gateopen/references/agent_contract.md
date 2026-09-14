# Fast-screen input/output contract

## Runtime

Use Python 3.11+ and requirements.txt; runtime-lock.json records actual versions.
Verify MANIFEST.json before loading models. MACE-MOF0-v2 uses `pbe_d3` with no
extra D3. PACMOF RF-DDEC requires sklearn 1.8.0. The included RASPA binary is
macOS arm64-specific with external dynamic libraries; source is included.

## Required inputs

Top-level keys: `temperature_K`, increasing positive `pressure_bar` (≥3 points),
`phases.closed`, `phases.open`, `guest`, `force_field`, `force_field_source`.
Each phase has `cif`, `cleaning_report`, `formula_units` (Z per CIF cell).
Paths are absolute or relative to the JSON. Use fresh output directories.
Read vendor/clean-cif instructions and [cleaning_gate.md](cleaning_gate.md);
source/target audits and chemistry evidence must be real, not automatic flags.

| Optional setting | Default/meaning |
|---|---|
| host.method | elastic-network; other methods rejected by fast entry |
| host.relax | false; pre-reviewed endpoint geometry |
| host.calculator | bundled MACE, pbe_d3, dispersion=false |
| host.proxy_mesh | [3,3,3]; compare finer q mesh |
| host.proxy_parameter_sets | five spring/packing scenarios |
| charges.method | pacmof; source requires charges.source |
| pacmof_root | bundled vendor/pacmof |
| reservoir | coolprop; ideal requires explicit scientific justification |
| raspa.binary | bundled binary; override on other platforms |
| raspa.cycles | 20000, not proof of convergence |
| raspa.initialization_cycles | 3000 |
| raspa.equilibration_cycles | 5000 |
| raspa.seeds | [417,941] |
| raspa.cutoff_A | 12; automatic box heights ≥2×cutoff |
| raspa.point_timeout_seconds | 1800 |
| henry_review | per-phase run-specific evidence strings |

Guest: `name`, neutral one-component `smiles`, `source` (adsorption model).
Default geometry route is RDKit ETKDG/MMFF94s/UFF. Alternatively supply
`geometry_sdf` and `geometry_source`: explicit-H, valid single-molecule SDF,
identity checked against SMILES. No electronic refinement in this bundle.
Use a sourced linear CO₂ geometry if the generic optimizer produces implausible
C=O bonds. Custom guests need `molecule_json` (rigid) and `atom_types` in exact
RDKit/SDF atom order including H. Geometry does not supply LJ/charges.

VOCs/custom gases require `gas_profile`: name/kind/source, optional coolprop_fluid,
and sourced Psat(T) or Antoine coefficients with temperature/pressure units and
validity interval. See [vapor_and_guest.md](vapor_and_guest.md). No mixtures,
electrolytes or flexible-guest ensemble sampling in this contract.

## Closed nonadsorbing branch

An explicit model replaces only the closed insertion/Henry gate:

```json
{"adsorption_models":{"closed":{
  "kind":"nonadsorbing", "interpretation":"conditional_scenario",
  "evidence":"Actual run-specific exclusion/Widom review and explicit assumption",
  "max_fugacity_bar":10.0,
  "sensitivity_K_molecules_fu_bar":[0.000001,0.00001,0.0001]
}}}
```

Numbers are schema illustrations, not default material data. The maximum
fugacity must cover the requested range. Residual coefficients assume N=Kf over
that range; they are not statistical bounds. Positive sampled closed loading
blocks a conflicting zero model. Open and host gates remain active.
`run_empty_widom.py --help` supplies the command to audit dilute Henry behavior
using verified RASPA inputs and independent seeds with all population-changing
moves disabled. Zero occupancy is enforced by that algorithm, not measured
exclusion; Henry cannot rule out finite-pressure cooperative filling.

## Numerical gates and interpretation

For sampled phases the first two mean N/f values must agree within 25% plus
henry_review evidence. This heuristic is necessary, not sufficient convergence.
Review independent seeds, longer runs and actual new pressure points near roots.
Native confidence limits are not relabeled standard deviations. Normalize by
Z×supercell product for extensive box quantities; use absolute, not excess,
loading. Integrate at the original simulation fugacity; include PΔV once.

Host spring/q-grid sensitivity, adsorption loading sensitivity and residual
uptake sensitivity diagnose fragile rankings; none is automatically a calibrated
confidence interval. The full runner propagates the five host parameter scenarios
to pressure and reports declared residual scenarios. To propagate other host terms or loading perturbations,
re-evaluate phase_diagram with the same source N(f); do not fit to a target root.
An elastic network's positive spectrum is not stability under the ML potential.

## Output contract

- resolved_spec.json: effective defaults and resolved input paths.
- run/identity.json: input/code fingerprint.
- run/<stage>/completed.json: result and artifact hashes; failure.json on error.
- host_free_energy/.../host_thermodynamics.json: ΔE, proxy ΔU/ΔS/ΔF, five parameter
  scenarios, per-state values and model provenance.
- raspa_closed|raspa_open/...: native inputs/outputs, actual f, absolute N/f.u.,
  energies, confidence, seeds and Widom diagnostics.
- analysis/.../result.json: blockers, fits, states, equilibrium crossings,
  residual-adsorption sensitivity when declared.
- run/workflow_summary.json: final machine-readable summary.

Exit 0 = completed conditional screen; exit 2 = incomplete prediction. A missing
root means no predicted opening in the sampled range only. Report P_eq, not
experimental P_open/P_close. Fits of isotherms/interaction energy do not replace
the raw adsorption integral; a single-T affinity is not an adsorption enthalpy.
