---
name: flexible-mof-gateopen
description: Agent-native fast screening of flexible-MOF equilibrium opening pressures from closed/open CIFs with clean-cif, PACMOF, static MACE, entropy proxy and RASPA GCMC. Supports nonadsorbing closed phases and saturation-limited VOCs; no MD, DFT or web interface.
---

# Flexible-MOF fast screening

Read [README.md](README.md) and [agent_contract.md](references/agent_contract.md).
Full theory: [principles.pdf](assets/principles.pdf), with searchable
[Markdown source](references/principles.zh-CN.md) and the reproducible
[XeLaTeX source](references/principles.tex).

1. Run `python scripts/verify_bundle.py`. Use hash-checked trusted model files;
   pickle/Torch weights can execute code. Install requirements.txt if needed;
   never suppress the PACMOF scikit-learn version check. The RASPA executable is
   macOS arm64-specific; other systems must build the supplied source.
2. Obtain closed/open CIFs, gas, T, Z per CIF cell, sourced adsorption force field
   and pressure range. Preserve originals. Use clean-cif: read the bundled
   `vendor/clean-cif/SKILL.md` and its cleaning decisions. Resolve disorder and
   pore solvent with chemistry review; preserve framework ligands, structural
   water and counterions. Require same-composition/activation P1 structures and
   hash-bound dual-parser/chemistry audits. See [cleaning_gate.md](references/cleaning_gate.md).
3. Use sourced neutral CIF charges or PACMOF RF-DDEC. Preserve coordinates/cell
   and audit neutrality/order. PACMOF predicts charges, not host free energy.
4. Prepare the guest with RDKit or a sourced valid SDF. Keep atom types, LJ,
   charges and rigid-geometry mapping explicit. Optimizer convergence is not
   geometry/model validation. No electronic refinement, MD, DFT, host relaxation
   or finite-displacement phonons is executed by the fast entry point.
5. Run `python scripts/run_fast_screen.py --spec /absolute/job.json --output-dir /absolute/new_run`.
   This combines static MACE-MOF0-v2 (`pbe_d3`, no extra D3), periodic elastic-
   network entropy, independent-seed GCMC, actual-fugacity integration and PΔV.
   Identical checkpointed runs may resume; changed inputs require a new directory.
   Use quick_host_free_energy.py alone for a host-only screen.
6. Review Henry behavior, seeds/cycles and root-region pressure resolution. VOCs
   must remain ≤ sourced Psat(T). Never extrapolate isotherms or relabel ideal
   simulations using real-gas fugacity after the fact.
7. Closed phases need not adsorb. Use an explicit zero-uptake branch plus residual
   sensitivity; run_empty_widom.py provides dilute evidence when needed. Zero
   accepted swaps is not proof of exclusion; a small pore window is not proof
   that interior cavities cannot adsorb. Widom does not rule out cooperative
   filling. Open-state sampling gates remain active.
8. Report all stable closed-to-open equilibrium crossings or sampled-range
   no-root outcomes, with host/adsorption/residual sensitivity. Refine actual
   simulations when needed; dense integration nodes are not new GCMC evidence.

## Scientific interpretation

Use open minus closed per identical formula unit:
`ΔF_screen=ΔE_static+ΔU_vib,proxy−TΔS_proxy`;
`ΔΩ=ΔF_screen+PΔV−RT∫(N_open−N_closed)dln(f)`.

The entropy proxy is uncalibrated and uses a different model from static MACE.
It is not self-consistent host F or proof of phase stability. Parameter scenarios
are not confidence intervals. Report P_eq, not experimental P_open/P_close;
endpoints do not determine barriers, kinetics or hysteresis. Fitting ΔF to a
gate pressure and reproducing that same point is calibration, not validation.

Deliver reviewed/charged CIFs, gas/model provenance, native isotherms and energy
fits, ΔE/ΔU/ΔS/ΔF, P_eq or explicit blocker/no-root status, and limitations.
