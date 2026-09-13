# X-pcu-5-Zn + CO2 at 268 K

This is an archived v0.1 dataset for testing analysis. It has not been certified
against the v0.2 mandatory cleaning/solvent gate. The source CIF contains alpha,
beta and gamma blocks; the historical scan used ordered alpha/open and
gamma/closed representatives. See [corrections](../../references/method_boundaries.md)
before using historical uncertainty or EOS scenarios.

| endpoint | P1 formula | atoms | formula units/cell | volume (A^3) |
|---|---:|---:|---:|---:|
| alpha/open | C64H44N4O16S8Zn4 | 140 | 2 | 2559.90 |
| gamma/closed | C128H88N8O32S16Zn8 | 280 | 4 | 3402.81 |

The PACMOF files preserve coordinates and occupancies and are neutral to machine precision. The PACMOF wrapper was also re-run from the public PACMOF training tables in a temporary environment; the regenerated alpha/gamma charges remained neutral and differed from the archived files by at most 0.0035/0.0111 e per atom. The serialized model binary is not published. RASPA used Garcia-Sanchez CO2, Ewald electrostatics, `UseChargesFrom: CIF_File`, and a clearly labelled UFF-like placeholder framework Lennard-Jones model.

The archived loading tables are sanitized to use relative `raspa_output/` paths; raw restart binaries and local output directories are intentionally not published.

The ordered representatives can be regenerated from the source block with:

```bash
python extract_ordered_xpcu5zn.py \
  --source structures/source_ja8b08642_si_001.cif \
  --output-dir /tmp/xpcu5zn_ordered
```

The most important result is a conditional two-state screen. At 19 bar the adsorption integral is `67.94 +/- 0.94 kJ per framework formula unit`. A `67.9 kJ` empty-host penalty is an explicit calibration to the reported 19 bar benchmark, not an independently calculated framework free energy.

The independent low-pressure refinement changes the 19-bar integral to `69.24 +/- 0.93 kJ per framework formula unit` and the same illustrative 67.9-kJ crossing to about `17.6 bar`. This is a useful numerical sensitivity, not a new physical prediction.

Read `reports/validation_report_268K.md` before reusing the number.
