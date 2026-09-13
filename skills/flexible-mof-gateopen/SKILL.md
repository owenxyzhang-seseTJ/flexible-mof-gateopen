---
name: flexible-mof-gateopen
description: Prepare clean flexible-MOF structures and guest molecules, run saturation-limited RASPA adsorption, and screen two-state osmotic equilibrium pressures with explicit host free energies.
---

# Flexible MOF adsorption

## 1. Prepare each framework state — mandatory

Use $clean-cif before assigning charges or preparing GCMC. Read its SKILL.md
and cleaning decisions. If unavailable, locate or install it before continuing
structure preparation. Analysis of existing isotherms remains usable.

Preserve the source, resolve disorder into a justified ordered configuration,
and identify/remove pore solvent. Distinguish coordinated solvent, structural
water and charge-balancing counterions; document the activation state. Never
delete framework fragments or set all disorder occupancies to one.

Require chemistry review, Pymatgen/ASE agreement, Spglib cross-check and
identity-only P1 output. Follow [cleaning_gate.md](references/cleaning_gate.md)
for solvent decisions and the hash-bound audit record. Run check_cleaning.py.
Unresolved disorder, solvent or chemistry blocks PACMOF and GCMC.

## 2. Prepare charges and guest

Use justified source charges or run run_pacmof.py with a cleaning report per
input. Verify neutrality, element order, coordinates, occupancies and cell
preservation; retain the charged CIF and updated record.
PACMOF predicts charges, not framework free energy.

For VOCs read [vapor_and_guest.md](references/vapor_and_guest.md).
prepare_guest_rdkit.py uses explicit H, ETKDGv3 and MMFF94s/UFF to select the
lowest-energy converged sampled conformer. Review stereochemistry and relevant
torsions. Supply validated RASPA atom types, charges and force-field parameters
separately; gas-phase optimization does not establish an adsorption model.

## 3. Run fixed-host adsorption

Specify gas, temperature and force-field provenance. VOCs require sourced
Psat(T) or Antoine coefficients with explicit units and a valid temperature
range. Enforce 0 < P <= Psat(T) in generation and analysis. Report P/Psat;
never extrapolate saturation data.

prepare_raspa_inputs.py requires a passed cleaning record, neutral charges,
explicit supercell and cutoff. Perpendicular supercell heights must be at
least twice the largest real-space cutoff. It writes Ewald/CIF-charge inputs
and rejects existing output directories. Custom guests require a gas profile
and matching molecule JSON. Read [raspa_input.md](references/raspa_input.md).

Launch when simulation is requested; --run requires a local RASPA executable.
Parse absolute loading and guest moves. Zero loading with no accepted
insertions requires an accessibility audit. Use independent seeds, sufficient
equilibration and longer repeats near a proposed crossing.

## 4. Screen equilibrium pressure

Normalize phases to the same formula unit:
DeltaOmega_ads = -RT integral (N_open - N_closed) d ln f.

Verify the Henry regime below the first point. Supply host free energy with
source, temperature, normalization and uncertainty. Account for external
P DeltaV or state its neglect. Static ML energy is a prior unless finite-
temperature entropy and material-specific accuracy have been validated.

quick_free_energy_screen.py gives required host penalty and a conditional
two-state crossing. Its CSV convention is ideal reservoir f=P; do not relabel
those loadings with EOS fugacities at unchanged pressures. Report nodal errors
and sensitivity to pressure grid, force field and host free energy.
Read [method_boundaries.md](references/method_boundaries.md).

Call the result P_eq. Experimental P_open/P_close, barriers and spinodals
require intermediate states and a justified landscape/criterion. Do not
validate a pressure against the same point used to fit the host term.

## Deliver

Provide source/derived CIFs, cleaning and charge records, guest/model provenance,
Psat and pressure grid, RASPA inputs/loadings, and the pressure result with its
assumptions. Keep binaries, credentials and restricted files out of public
packages. During skill maintenance, run tests, skill validation and the bundled
example; synchronize repository and installed skill scripts.
