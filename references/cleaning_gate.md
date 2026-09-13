# Required structure preparation

Use `$clean-cif` before PACMOF or GCMC for **every framework state**.
Preserve the source, select its data block, inspect raw loops, classify disorder,
identify solvent, and write a new computation-ready P1 model. Run the clean-cif
audit on the selected source block and candidate with Pymatgen, ASE and Spglib.

Choose complete, compositionally defensible disorder configurations. Never keep
all alternatives at occupancy one. Identify guests using periodic connectivity
and chemical identity: remove unbound pore solvent; distinguish coordinated
solvent, counterions and structural water. Retain framework ligands and required
counterions. Coordinated-solvent removal requires an explicit activation model
and a new coordination/charge audit. Report SQUEEZE/masked solvent even when no
solvent atom rows exist. If the activation state is unresolved, stop before GCMC.
No automatic energy minimization of the experimental framework is implied.

After chemistry review, create this record next to the two audit artifacts.
Replace every placeholder with actual evidence. `passed` is a review conclusion,
not a flag to set merely to run a calculation.

```json
{
  "schema_version": 1,
  "status": "passed",
  "source_sha256": "SHA256_OF_PRESERVED_SOURCE",
  "candidate_sha256": "SHA256_OF_CANDIDATE",
  "transformation_record": "Selected block, source labels removed/retained, symmetry operations, disorder selection, formula and Z before/after",
  "source_audit": {"path": "source_audit.json", "sha256": "SHA256_OF_AUDIT"},
  "candidate_audit": {"path": "candidate_audit.json", "sha256": "SHA256_OF_AUDIT"},
  "disorder": {"status": "ordered", "reason": "Actual configuration and justification"},
  "solvent": {"status": "removed", "reason": "Actual species and labels removed; retained coordinated species/counterions and reason"},
  "chemistry": {
    "connectivity": {"status": "passed", "evidence": "Periodic framework and complete ligands checked"},
    "hydrogen_hosts": {"status": "passed", "evidence": "Each H has its intended bonded host"},
    "metal_coordination": {"status": "passed", "evidence": "Actual coordination counts and distance ranges"},
    "nonbonded_contacts": {"status": "passed", "evidence": "Shortest nonbonded pairs reviewed"},
    "composition_charge_balance": {"status": "passed", "evidence": "Actual formula, counterions, charge balance and formula units"}
  }
}
```

Disorder statuses: `absent`, `ordered`. Solvent statuses: `absent`, `removed`,
`retained_with_reason`. Chemistry checks allow `not_applicable` only with a reason.
Use paths relative to the record for portable artifacts. The gate verifies
artifact hashes, raw occupancies/labels/P1, dual-parser composition and cell,
periodic near-duplicates, and Spglib availability. It does not automate the
chemical judgments recorded above. PACMOF rebinds the report to the charged CIF
only after checking element order, coordinates, occupancies and cell preservation.

```bash
python scripts/check_cleaning.py --cif host_clean.cif --report host.cleaning.json
```
