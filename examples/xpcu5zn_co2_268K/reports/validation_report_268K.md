# X-pcu-5-Zn 268 K validation report

This report audits the existing fixed-host RASPA result; it does not turn a rigid-host calculation into a flexible-framework kinetic simulation.

- Checks: {'pass': 29, 'warn': 1, 'fail': 0}
- Ideal-fugacity required host penalty at 19 bar: **67.936 ± 0.940 kJ per framework formula unit**.
- PR post-processing correction: φ(19 bar)=0.8610; calibrated 67.9-kJ scenario root **22.34 bar**.
- SRK post-processing correction: φ(19 bar)=0.8706; calibrated 67.9-kJ scenario root **22.03 bar**.
- Independent low-pressure grid: **69.237 ± 0.930 kJ per formula unit** after replacing the main grid through 1 bar; the 67.9-kJ scenario root is **17.59 bar**.

## What is verified

- Source SI page 39 contains the reported β→open 16 bar / 11 bar and γ→open 19 bar / 11 bar values.
- Gemmi and ASE agree on charged-CIF atom counts and cell volumes; PACMOF charge sums are effectively neutral and coordinates are unchanged after charge assignment.
- Every main CSV point re-parses from its RASPA text output; production inputs use Ewald and `UseChargesFrom: CIF_File`.
- Formula-unit normalization is independently reproduced from the raw CSV.
- Independent short and long RASPA repeats are recorded under `validation_replicates_268K*` and `validation_long_268K_*`.
- The gamma/closed long repeat records zero absolute loading together with no constructed translation moves; this is retained as an accessibility/setup flag, not hidden as a zero-uncertainty measurement.

## What remains unverified

- **host_free_energy:** No independently computed empty-host F_host(alpha)-F_host(gamma) is present; approximately 67.9 kJ/formula unit is a calibrated scenario, not an independent prediction.
- **flexibility:** All production and repeat RASPA calculations keep the framework rigid; no q-grid, spinodal, saddle, or barrier was sampled.
- **force_field:** Framework LJ parameters are explicitly UFF-like placeholders; PACMOF charges are an electrostatic input, not a host free-energy model.
- **fugacity:** Production inputs use FugacityCoefficient=1.0; PR/SRK corrections are post-processing sensitivities.
- **kinetics:** Monte Carlo steps are not physical time; P_open/P_close hysteresis and rates are not identified by this dataset.
- **intermediate_structures:** Only alpha and gamma ordered endpoints are available; no chemically validated intermediate X-pcu-5-Zn structures were supplied for a multi-rigid q-grid.
- **accessibility:** The gamma/closed 19-bar long repeat has zero absolute loading and no constructed translation moves; treat this as an explicit accessibility/setup audit flag rather than a zero-uncertainty observation.

## Interpretation

The defensible result is a conditional osmotic equilibrium crossing. With an explicitly assumed empty-host penalty of 67.9 kJ per framework formula unit, the ideal-fugacity analysis gives about 19.0 bar; PR/SRK fugacity post-processing shifts that scenario to roughly 22.0–22.3 bar. Neither value is an independent prediction of experimental adsorption/desorption hysteresis.

Detailed machine-readable evidence is in `validation_report_268K.json`.
