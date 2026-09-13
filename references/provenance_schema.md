# Provenance record

Each public example should include a small JSON record with:

```json
{
  "material": "…",
  "source_url": "https://…",
  "source_doi": "10.…",
  "retrieved_utc": "YYYY-MM-DD",
  "source_sha256": "…",
  "derived_structure_policy": "…",
  "charge_method": "PACMOF RF-DDEC",
  "simulation_temperature_K": 268.0,
  "gas": "CO2",
  "limitations": ["…"]
}
```

Use URLs and hashes instead of embedding restricted supporting information or
personal absolute paths.  Derived ordered P1 CIFs must say that they are
ordered representatives rather than newly refined structures.
