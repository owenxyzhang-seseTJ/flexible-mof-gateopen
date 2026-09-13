# Method boundaries

## v0.2 corrections to archived examples

- The ideal-reservoir X-pcu-5-Zn integral is 67.9285 kJ/mol of formula units at
  19 bar. Correct propagation combines weights of shared integration nodes;
  its propagated input-error magnitude is 1.2913 kJ/mol, replacing 0.9397.
  Input RASPA error bars retain their original convention; this is not
  automatically a one-sigma error or a complete model uncertainty.
- Previous 22-bar PR/SRK scenarios are withdrawn. The dataset was sampled at
  f=P. Integrating those same N(P) values against a new EOS f(P) is inconsistent.
  Obtain N at the actual fugacities (or remap the reservoir consistently), then
  integrate. CO2 at 268 K also requires checking stable vapor/liquid branches;
  monotonic fugacity alone is not a stable-phase test.
- The 67.9-kJ host term was calibrated from the 19-bar switching point.
  Agreement at that pressure is not independent prediction or validation.
- Existing reports and plots are historical records; their original pass counts
  do not certify the new cleaning/solvent gate, force field or host free energy.
  Gamma has zero accepted swaps but nonzero constructed insertion attempts.
  Its zero-loading uncertainty requires an accessibility audit.
- The helper neglects external mechanical P DeltaV. If needed, add it in units
  consistent with the host term before solving the osmotic crossing.

## What fixed-host GCMC provides

For a fixed framework geometry (q), RASPA GCMC provides an adsorbate
grand-canonical response.  An adsorption contribution can be reconstructed as

\[
\Delta\Omega_\mathrm{ads}(f)=-RT\int_0^f
  [N_\mathrm{open}(f')-N_\mathrm{closed}(f')]\,d\ln f'.
\]

Normalize each loading to the same framework formula unit before integrating.
Below the first measured pressure use the Henry approximation only after
checking the low-pressure regime, (N(f)\propto f). At
elevated pressure use fugacity rather than silently treating pressure as
fugacity.

## What it does not provide

PACMOF charges are electrostatic input data; they are not a framework free
energy.  A rigid-host GCMC scan does not sample linker, bond, cell, or pore
degrees of freedom and does not identify a barrier, spinodal, switching rate,
or experimental (P_\mathrm{open})/(P_\mathrm{close}).

The two-state equilibrium scenario is

\[
\Delta\Omega_\mathrm{os}(f)=\Delta F_\mathrm{host}(T)+
\Delta\Omega_\mathrm{ads}(f).
\]

The empty-host term must come from experiment, an explicitly calibrated model,
a validated classical potential plus entropy treatment, or an ML potential
that has been benchmarked for the material.  A static ML energy difference is
not a finite-temperature free energy.

## Recommended evidence ladder

1. Public structure provenance and CIF/Gemmi/ASE audit.
2. Explicit disorder/order policy and formula-unit normalization.
3. PACMOF charge file with charge-loop and neutrality audit.
4. RASPA input audit: `ChargeMethod: Ewald` and
   `UseChargesFrom: CIF_File`.
5. Raw loading/error table plus move/accessibility diagnostics.
6. Henry-corrected adsorption integral and quadrature sensitivity.
7. Explicit host-free-energy source and uncertainty.
8. Only with a validated (q)-grid or enhanced sampling: spinodal/barrier
   claims.
