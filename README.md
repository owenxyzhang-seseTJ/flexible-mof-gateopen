# Flexible MOF Gate Opening

[中文](README.zh-CN.md) | English

A Codex skill for preparing flexible-MOF adsorption calculations and screening equilibrium transition pressures.

**Clean CIF → charges → guest model → RASPA → adsorption + host free energy → P_eq**

## Install

```bash
git clone https://github.com/owenxyzhang-seseTJ/flexible-mof-gateopen.git
cd flexible-mof-gateopen
python -m pip install -e '.[dev,structure,guest]'
cp -R skills/flexible-mof-gateopen ~/.codex/skills/
```

Also install **clean-cif**. RASPA3 and PACMOF are separate installations. Core free-energy analysis needs only Python 3.10+.

## Use in Codex

> Use $flexible-mof-gateopen for my open/closed CIFs and ethanol adsorption at 298.15 K. First use $clean-cif to resolve disorder and solvent, assign PACMOF charges, prepare the guest with RDKit, and scan only up to sourced Psat(T). Report P_eq for my supplied host free-energy range.

Provide both structures, gas, temperature, a validated guest/force-field model, and host free energy or a clearly labeled range to test.

| Step | Requirement / output |
|---|---|
| Framework | Mandatory disorder/solvent review, P1 CIF, dual-parser audit and hash-bound record |
| Charges | Source charges or PACMOF; neutral cell and preserved geometry |
| Guest | RDKit conformers → MMFF94s/UFF → SDF, XYZ and convergence record |
| VOC pressure | Sourced Psat at the actual temperature; reject P > Psat |
| GCMC | Explicit supercell/cutoff, Ewald and CIF charges; matching model for custom guests |
| Free energy | Formula-unit-normalized absolute loading → required host penalty and conditional P_eq |

Preparation stops if the cleaning record is missing, stale or unresolved. RDKit geometry still requires separate RASPA force-field parameters.

## Try the bundled analysis

```bash
python scripts/quick_free_energy_screen.py \
  --csv examples/xpcu5zn_co2_268K/data/fixed_host_isotherms_268K.csv \
  --temperature 268 --open-formula-units 2 --closed-formula-units 4 \
  --target-pressure 19 --host-penalty 67.9
```

The archived X-pcu-5-Zn/CO2 data give a required host penalty of **67.93 kJ/mol of framework formula units** at 19 bar and **P_eq ≈ 18.96 bar** for an assumed 67.9 kJ/mol penalty. This legacy analysis example has not been certified by the new mandatory cleaning gate. See [numerical corrections](references/method_boundaries.md).

Experimental opening/closing pressures require a validated framework free-energy landscape and barrier or instability criterion. Two rigid structures alone do not determine hysteresis or switching rates.

## Details

- [Structure and solvent preparation](references/cleaning_gate.md)
- [VOC saturation and RDKit](references/vapor_and_guest.md)
- [RASPA commands](references/raspa_input.md)
- [Skill instructions](skills/flexible-mof-gateopen/SKILL.md)

```bash
python -m pytest -q
python scripts/validate_example.py
python skills/flexible-mof-gateopen/scripts/quick_validate.py skills/flexible-mof-gateopen
```
