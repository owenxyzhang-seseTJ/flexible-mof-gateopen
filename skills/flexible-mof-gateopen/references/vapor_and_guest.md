# Single-component vapor and guest inputs

For VOCs enforce `0 < P <= Psat(T)` both when preparing RASPA and analyzing an
isotherm. An over-limit request fails before creating simulation directories;
it is not silently truncated. Below-saturation pore condensation remains allowed.
For a mixture this single-component helper is insufficient: use partial
fugacities and an appropriate phase-equilibrium model.

A gas profile contains `name`, `kind` (`VOC` or `non-VOC`), and `source`.
Non-CO2 gases also require a RASPA molecule JSON and a matching force field.
For a VOC include a sourced saturation record:

```json
{
  "name": "VOC_NAME",
  "kind": "VOC",
  "source": "Guest identity and force-field reference",
  "saturation": {
    "method": "value",
    "temperature_K": 298.15,
    "pressure_bar": 0.1,
    "source": "REPLACE: measured or validated Psat reference at exactly this temperature"
  }
}
```

The number above is a schema illustration, not data for a named VOC. Alternatively
use `method: antoine`, numeric `A`, `B`, `C`, `temperature_unit: K` or `C`,
`pressure_unit: bar`, `Pa`, `kPa` or `mmHg`, and `valid_temperature_K: [Tmin,Tmax]`.
The implemented convention is `log10(P) = A - B/(T + C)` in those explicit units.
Never extrapolate coefficients; provide `critical_temperature_K` when available.
No saturation pressure is assumed for supercritical conditions.

Generate guest geometry with RDKit:

```bash
python scripts/prepare_guest_rdkit.py --smiles CCO --output-dir ethanol_geometry
```

The script adds H, uses seeded ETKDGv3 conformers, optimizes with MMFF94s or UFF,
and exports the lowest-energy converged sampled geometry as SDF/XYZ plus JSON.
Verify stereochemistry and relevant torsional conformers. This is a gas-phase
geometry preparation step; energies are not framework free energies. Supply
validated RASPA atom types, charges, LJ and any flexible intramolecular terms
separately. Do not silently transfer MMFF charges or use a CO2 definition for VOCs.
