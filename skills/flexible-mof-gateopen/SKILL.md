---
name: flexible-mof-gateopen
description: Build and audit a public, reproducible flexible-MOF gate-opening workflow from CIF provenance through PACMOF charges, RASPA fixed-host GCMC, quick osmotic free-energy screening, and explicitly bounded gate-pressure estimates.
---

# Flexible-MOF gate-opening workflow

Use this skill when the user wants to download or inspect public CIFs, resolve
disorder into a reproducible representative, add PACMOF charges, prepare or
run RASPA GCMC, integrate adsorption free energies, screen framework free
energy models, or estimate a flexible-MOF gate-opening pressure.

## Core scientific boundary

Keep three quantities separate in every output:

1. `Omega_ads`: adsorbate grand-potential contribution from fixed-host GCMC;
2. `F_host`: empty-framework free-energy difference, which must be an explicit
   input, calibration, or independently validated model;
3. kinetic/barrier quantities (`P_open`, `P_close`, spinodals, rates), which
   require a framework reaction coordinate and enhanced sampling or an
   experiment-calibrated criterion.

PACMOF charges are electrostatic input, not `F_host`.  A rigid GCMC result is
not a flexible-framework kinetic prediction.  Monte Carlo cycles are not
experimental time.

Read [references/method_boundaries.md](references/method_boundaries.md) when
writing a methods section or deciding whether a claim is supported.

## Workflow

### 1. Provenance and structure intake

- Prefer a public DOI/URL and preserve the source hash in a provenance record.
- Select the intended CIF data block before parsing.  Record space group,
  formula-unit count `Z`, cell volume, atom counts, disorder/occupancy policy,
  and whether the output is a derived ordered P1 representative.
- Never call an ordered representative a new structure solution or refinement.
- Audit with Gemmi (and ASE when available); check occupancies, formula, cell,
  atom count, and coordinate preservation.
- Do not upload restricted PDFs, private logs, credentials, RASPA binaries,
  restart binaries, or personal absolute paths to a public repository.

### 2. Charges

- Keep the source CIF unchanged and write a separate charged CIF.
- PACMOF is optional but appropriate when the original CIF lacks charges.
- Verify `_atom_site_charge` is present, the atom ordering/coordinates are
  unchanged, and the cell charge is neutral to a stated tolerance.
- Report the PACMOF model/version if it is available; do not reinterpret its
  charges as a framework energy model.

### 3. RASPA fixed-host GCMC

- Require an explicit gas model, temperature, pressure/fugacity convention,
  force-field provenance, cycle counts, seeds, and formula-unit normalization.
- For charged frameworks use `ChargeMethod: Ewald` and
  `UseChargesFrom: CIF_File`; audit the actual `simulation.json`, not just a
  template.
- Keep framework degrees of freedom fixed in the quick screen.  State that
  this samples guest translation/rotation/reinsertion/swap only.
- Parse absolute and excess loading, statistical error, and move diagnostics
  from raw text.  A zero loading with no constructed translations or no
  accepted guest moves is an accessibility/setup flag that needs an audit
  before quantitative use.
- Run at least one independent seed and a longer repeat near a proposed
  crossing when the result is intended as evidence.

Read [references/raspa_input.md](references/raspa_input.md) when preparing or
auditing RASPA directories.

### 4. Fast free-energy screening

For paired open/closed fixed-host isotherms, normalize each loading to one
framework formula unit and calculate

`DeltaOmega_ads(f) = -R*T*integral[(N_open-N_closed) d ln(f)]`.

Use the exact Henry-limit contribution below the first pressure point.  Do not
linearly extrapolate loading in `ln(P)` to an artificial zero-pressure point.
Use fugacity for non-ideal reservoir sensitivity and label post-processing EOS
corrections as such.

The fastest defensible ranking is:

1. cheap static/ML host-energy prior and pore-accessibility audit;
2. 3–5 log-spaced endpoint GCMC points for shortlisted candidates;
3. full endpoint scan and independent repeats for the top candidates;
4. multi-`q` rigid states only after the endpoint screen passes.

The skill's `quick_free_energy_screen.py` reports the required host penalty and
an optional scenario crossing.  A host penalty calibrated from a known gate
pressure is a conditional calibration, not an independent prediction.

### 5. Gate-pressure language

- `P_eq`: equilibrium crossing of two explicitly defined osmotic states.
- `P_open`/`P_close`: do not report as independently predicted unless a
  framework free-energy landscape, spinodal/barrier criterion, or validated
  experimental calibration is present.
- `spinodal`: local-minimum loss along a stated reaction coordinate.
- `barrier`: a free-energy difference, not a Monte Carlo time or rate.

Always include an uncertainty/sensitivity table for `F_host`, force field,
fugacity, pressure grid, and formula-unit normalization.

## Reusable scripts

The public mirror under `skills/flexible-mof-gateopen/` contains portable
copies of these helpers:

- `fetch_cif.py`: download a public file and verify SHA-256;
- `audit_cif.py`: Gemmi structure/charge audit;
- `run_pacmof.py`: invoke a local PACMOF installation without overwriting the
  source CIF;
- `prepare_raspa_inputs.py`: write portable fixed-host RASPA3 inputs;
- `parse_raspa_output.py`: parse loadings and accessibility diagnostics;
- `quick_free_energy_screen.py`: Henry-corrected two-state screening and
  pressure crossing.

These scripts default to read-only preparation or write-only derived outputs;
running RASPA requires an explicit `--run` and a user-supplied executable.

## Public-package quality gate

Before claiming completion, run:

```bash
python scripts/quick_validate.py /path/to/skill
python -m pytest -q
```

Then run a bundled example without network access and inspect the generated
JSON.  For a real RASPA case, archive the input JSON, charged-CIF hash, raw
loading table, validation report, and limitations; do not archive restart
binaries.  Use the public provenance schema in
`references/provenance_schema.md`.
