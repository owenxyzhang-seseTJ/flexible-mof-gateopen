# Portable RASPA3 input contract

The helper script writes a fixed-host, single-component CO2 GCMC directory.
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
