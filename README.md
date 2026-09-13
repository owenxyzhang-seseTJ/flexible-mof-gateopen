# flexible-mof-gateopen

An auditable, lightweight workflow for screening adsorption-induced gate opening in flexible metal-organic frameworks:

```text
public CIF + provenance
        -> ordered/disorder audit
        -> PACMOF charges (separate CIF)
        -> RASPA3 fixed-host GCMC
        -> Henry/fugacity-corrected adsorption free energy
        -> explicit host-free-energy scenario
        -> conditional osmotic equilibrium pressure
```

The package is deliberately honest about scope. A fixed-host GCMC scan gives an adsorbate grand-potential contribution. It does not by itself give the empty-framework free energy, a framework barrier, a spinodal, a switching rate, or an experimental `P_open`/`P_close`. Those quantities require a validated framework free-energy model and, for hysteresis, intermediate states or enhanced sampling.

## Included examples

### X-pcu-5-Zn + CO2, 268 K

[`examples/xpcu5zn_co2_268K/`](examples/xpcu5zn_co2_268K/) contains a public Figshare CIF, ordered alpha/gamma P1 representatives, PACMOF-charged CIFs, portable force-field/CO2 definitions, loading tables, independent repeats, the low-pressure grid, and a validation report. The source DOI is `10.1021/jacs.8b08642` and the public Figshare CIF is:

<https://ndownloader.figshare.com/files/13528844>

At 268 K, the fixed-host adsorption integral requires about `67.94 +/- 0.94 kJ per framework formula unit` to make the alpha/open and gamma/closed states cross at 19 bar under ideal fugacity. Taking `67.9` kJ per formula unit as an explicit calibration gives a conditional crossing of about `18.97 bar`; PR/SRK post-processing gives about `22.34/22.03 bar`. These are screening/scenario values, not independent kinetic predictions.

Replacing the main grid through 1 bar with the independent low-pressure grid gives `69.24 +/- 0.93 kJ per formula unit` and a `67.9-kJ` scenario root near `17.6 bar`; this grid sensitivity is retained rather than folded into a falsely precise single pressure.

The gamma/closed calculation has zero accepted CO2 loading, zero constructed translation moves, and zero accepted swap moves in the archived RASPA output. The example therefore labels this as an accessibility/setup audit flag as well as a possible non-porous state; it must not be hidden behind a zero statistical error.

A native local RASPA3 smoke test is summarized in [`examples/xpcu5zn_co2_268K/reports/raspa_input_smoke.json`](examples/xpcu5zn_co2_268K/reports/raspa_input_smoke.json). It verifies that the charged CIFs, `ChargeMethod: Ewald`, `UseChargesFrom: CIF_File`, guest moves, and text-output parser work with a macOS arm64 RASPA3 binary. The executable, raw logs, and restart files are intentionally not distributed.

### MIL-53(Al) + CO2

[`examples/mil53al_co2/`](examples/mil53al_co2/) contains the published-parameter two-state osmotic benchmark from Boutin et al. At 223 K its higher equilibrium branch is approximately `0.01696 bar`. This is a literature parameter benchmark, not a validated atomistic flexible-host simulation.

## Install and test

The core free-energy functions use only the Python standard library. Gemmi is needed for CIF auditing; RASPA3 and PACMOF are external tools and are not bundled.

```bash
python -m pip install -e '.[dev,cif]'
python -m pytest -q
```

Run the bundled screen without network access:

```bash
python scripts/quick_free_energy_screen.py \
  --csv examples/xpcu5zn_co2_268K/data/fixed_host_isotherms_268K.csv \
  --temperature 268 \
  --open-phase open_alpha --closed-phase closed_gamma \
  --open-formula-units 2 --closed-formula-units 4 \
  --target-pressure 19 --host-penalty 67.9
```

Audit the two charged structures:

```bash
python scripts/audit_cif.py examples/xpcu5zn_co2_268K/structures/alpha_pacmof_charged.cif
python scripts/audit_cif.py examples/xpcu5zn_co2_268K/structures/gamma_pacmof_charged.cif
```

## From a new public CIF to RASPA

1. Download with `scripts/fetch_cif.py` and record the URL, DOI, retrieval date, and SHA-256.
2. Select the correct CIF data block and resolve disorder explicitly. Keep the source file and write a new ordered representative.
3. Run `scripts/run_pacmof.py` if the source lacks charges; verify neutrality and coordinate preservation.
4. Use `scripts/prepare_raspa_inputs.py` to write fixed-host inputs. It enforces Ewald electrostatics and `UseChargesFrom: CIF_File`; add `--run` only when a local RASPA executable is intentionally supplied.
5. Parse each output with `scripts/parse_raspa_output.py`, inspect move diagnostics, and retain independent seeds near a crossing.
6. Run `scripts/quick_free_energy_screen.py`. Supply `F_host` explicitly or label it as a calibration/scenario. Use a multi-`q` workflow before making spinodal/barrier or hysteresis claims.

See [`references/method_boundaries.md`](references/method_boundaries.md) and the mirrored Codex skill at [`skills/flexible-mof-gateopen/SKILL.md`](skills/flexible-mof-gateopen/SKILL.md).

To install the skill locally for Codex, copy the directory
`skills/flexible-mof-gateopen/` into `~/.codex/skills/` (or the configured
`CODEX_HOME/skills/` directory), then run its `quick_validate.py` helper.

## Public-data policy

This repository intentionally excludes RASPA binaries, restart binaries, machine-specific absolute paths, private logs, credentials, and restricted supporting-information PDFs. Public-source provenance and derived-file hashes are kept in the example metadata.
