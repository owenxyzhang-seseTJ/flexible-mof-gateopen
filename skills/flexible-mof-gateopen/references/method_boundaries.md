# Method boundaries

## What fixed-host GCMC provides

For a fixed framework geometry (q), RASPA GCMC provides an adsorbate
grand-canonical response.  An adsorption contribution can be reconstructed as

\[
\Delta\Omega_\mathrm{ads}(f)=-RT\int_0^f
  [N_\mathrm{open}(f')-N_\mathrm{closed}(f')]\,d\ln f'.
\]

Normalize each loading to the same framework formula unit before integrating.
At the first measured pressure use the Henry limit, (N(f)\propto f).  At
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
