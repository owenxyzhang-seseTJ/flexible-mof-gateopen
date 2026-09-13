# Portable RASPA3 input contract

```bash
python scripts/prepare_raspa_inputs.py \
  --cif host_charged.cif --cleaning-report host_charged.cleaning.json \
  --force-field force_field.json --gas-profile ethanol.json \
  --molecule-json ethanol_raspa.json --temperature 298.15 \
  --pressures 0.001 0.01 --unit-cells 2 2 2 --cutoff 12 \
  --output-dir ethanol_runs
```

Example pressures/supercell must be checked against your sourced Psat and cell.
Omit gas-profile/molecule-json only for the bundled CO2 model. Add --run and
--raspa-bin /path/to/raspa3 to launch an authorized local calculation.
Run scripts/prepare_raspa_inputs.py --help for required arguments.

The helper writes fixed-host, single-component GCMC directories.
The framework CIF must contain `_atom_site_charge` values.  The input must
contain:

```json
{
  "ChargeMethod": "Ewald",
  "UseChargesFrom": "CIF_File"
}
```

The framework is rigid in this screening mode.  The default cycle count in the
template is a smoke-test value; increase it and run independent seeds for
quantitative uncertainty.  Keep the pressure/fugacity convention explicit in
the report.  Do not upload RASPA executables, restart binaries, private logs,
or machine-specific absolute paths.

Before treating a zero loading as meaningful, inspect the text output for
constructed and accepted insertion/swap/translation moves. A zero loading
with no constructed translations or no accepted swaps is an
accessibility/setup flag and should trigger a pore geometry, molecule-size,
and force-field audit.
