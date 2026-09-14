---
name: clean-cif
description: Clean, repair, symmetry-expand, and validate crystallographic CIF files. Use when Codex must expand an asymmetric unit to P1, remove periodic display copies or duplicated symmetry sites, resolve overlap caused by pre-expanded molecules or disorder, normalize labels and occupancies, check bonding and short contacts, prepare computation-ready crystal structures, or diagnose why a CIF fails in Pymatgen, ASE, Spglib, VESTA, Materials Studio, or other crystallographic software.
---

# Clean CIF

Produce a new, traceable CIF rather than overwriting the source. Treat symmetry,
occupancy, molecular completeness, and periodic geometry as separate questions.

## Required workflow

1. Preserve the source and inspect its raw CIF loops before parsing it as a
   structure. Record cell parameters, declared space group, symmetry operations,
   atom rows, labels, elements, coordinates, occupancies, and disorder fields.
2. Parse with Pymatgen and ASE independently. Do not treat one successful parser
   as proof that the structure is chemically valid.
3. Classify the input before expanding:
   - valid asymmetric unit;
   - partially or fully pre-expanded structure;
   - boundary/display images outside `[0, 1)`;
   - complete molecule crossing a special position;
   - crystallographic disorder or partial occupancy;
   - malformed CIF syntax or metadata only.
4. Build the P1 model using the applicable path in
   [references/cleaning-decisions.md](references/cleaning-decisions.md). Never
   delete atoms solely because two coordinates are close; confirm their elements,
   symmetry relationship, molecular role, occupancy, and connectivity.
5. Wrap final fractional coordinates into `[0, 1)`, assign globally unique atom
   labels, write occupancy explicitly, and declare only `P 1` with the identity
   operation. Keep any source-label mapping in an audit file rather than a custom
   CIF loop when viewer compatibility matters.
6. Run `scripts/audit_cif.py` on both the source and candidate output. Inspect the
   shortest contacts, formula, atom count, parser agreement, recovered symmetry,
   and coordination environment.
7. Perform chemistry-aware validation. Check expected covalent bonds, metal
   coordination, complete molecules across periodic boundaries, and nonbonded
   clashes. Do not "fix" experimental bond lengths by optimization unless the
   user explicitly requests geometry optimization.
8. Deliver the cleaned CIF, a concise audit report, and a reproducible script or
   transformation record. State every scientifically consequential choice, such
   as selecting one ordered disorder configuration or lowering symmetry.

## Expansion rules

- For a genuine asymmetric unit, apply all declared symmetry operations and
  deduplicate periodic equivalents of the same crystallographic site.
- For pre-expanded display atoms, identify periodic or symmetry-equivalent copies
  before applying additional operations.
- For a complete molecule that overlaps its symmetry mate, expand molecular
  fragments rather than isolated atoms. Determine the stabilizer, cosets, or a
  group-closed subgroup; require a clash-free, compositionally consistent model.
- For partial occupancy or disorder, preserve the experimental model unless a
  computation-ready ordered structure is requested. When ordering is necessary,
  generate and document a defensible configuration; do not silently set every
  occupancy to one.
- Preserve chemically meaningful periodic connectivity. Coordinates at `x=1`
  and `x=0` can be the same site, but atoms near opposite boundaries can also be
  distinct bonded images.

## Validation gate

Do not call a CIF clean until all applicable checks pass:

- source is preserved and output is a new file;
- Pymatgen and ASE agree on atom count and composition;
- no unexplained sub-angstrom duplicate contacts remain;
- formula change, if any, is explained by removed duplicates or disorder ordering;
- each expected atom has plausible connectivity and each H has one intended host;
- metal coordination counts and distances are chemically plausible;
- final CIF contains unique labels, valid occupancies, identity-only P1 symmetry,
  and conservative standard tags;
- Spglib recovery is reported as an informational cross-check, not used to
  overwrite the requested P1 declaration.

## Audit utility

Run:

```bash
python3 scripts/audit_cif.py input.cif --json audit.json --shortest 30
```

The utility is read-only. It reports raw-loop metadata, independent parser
results, minimum-image contacts, and Spglib symmetry. Use its output as evidence;
interpret contacts with the chemical rules in the reference file.

If the report marks Pymatgen, ASE, or Spglib unavailable, switch to a configured
Python environment containing that dependency and rerun. Do not claim dual-parser
validation from a partial report.
