# CIF cleaning decisions

## 1. Diagnose the representation

Inspect the raw atom loop before accepting parser expansion.

| Signal | Likely meaning | Action |
|---|---|---|
| Repeated labels, coordinates at `1`, negative coordinates | Display or boundary images | Compare modulo lattice translations; remove only confirmed copies |
| Pymatgen occupancy sum greater than one | Duplicate symmetry sites, disorder, or both | Inspect site orbits and occupancies; do not rescale blindly |
| One complete molecule plus non-P1 symmetry | Molecule may cross a special position or encode one disorder orientation | Expand whole fragments and test group-consistent ordered models |
| Several sites map within 0.01–0.05 Å | Exact/numerical duplicates | Merge only when element, orbit, and chemical role agree |
| Symmetry expansion gives 0.05–1.0 Å contacts | Pre-expansion or incompatible disorder configurations | Stop atom-wise expansion and reconstruct fragment/orbit logic |
| Occupancy below one or alternate-location fields | Crystallographic disorder | Preserve, enumerate ordered configurations, or request scientific direction |

## 2. Choose the transformation

### Valid asymmetric unit

Apply all symmetry operations. Wrap fractional coordinates and cluster periodic
equivalents by element and crystallographic site. Use a tight Cartesian tolerance,
normally 0.01 Å for exact generated coordinates and up to 0.05 Å only when source
rounding justifies it.

### Pre-expanded or viewer-oriented file

Identify which rows are source sites and which are periodic/symmetry images. Use
orbit membership and integer lattice translations, not labels alone. Retain one
representative per genuine site before expansion.

### Complete molecule on a special position

Build its intramolecular bond graph first. Apply candidate symmetry operations to
the complete fragment. Group transformations that reproduce the same fragment.
Select coset representatives or a group-closed subgroup only when it:

- gives the expected number of molecules and formula;
- preserves whole molecular connectivity;
- yields consistent coordination at every symmetry-related metal;
- removes impossible intermolecular contacts;
- can be described as an ordered configuration.

Record the chosen operation IDs and recovered subgroup. Other ordered
configurations can be scientifically distinct even when their local geometry is
similar.

### Disorder and partial occupancy

Do not convert fractional occupancy to full occupancy without justification. For
DFT or force-field calculations, generate one or more ordered supercell models
that satisfy composition, charge balance, local coordination, and minimum-image
separation. State that the output is an ordered model rather than the unique
experimental structure.

## 3. Geometry checks

Treat numeric cutoffs as flags, not universal chemistry laws.

- Any unrelated heavy-atom pair below about 1.0 Å is usually a duplicate or clash.
- Nonbonded H-H below about 1.2 Å is usually suspicious; molecular H2 is an
  exception.
- Typical placed C-H distances are about 0.9–1.1 Å.
- Common C-C, C-N, and C-O bonds are roughly 1.15–1.70 Å, depending on bond order.
- Metal-ligand distances are element- and oxidation-state-specific. Inspect the
  full coordination shell rather than applying one universal cutoff.
- Check bonds across every periodic boundary. A visually fragmented molecule may
  still be chemically complete through a lattice translation.
- Compare coordination counts across symmetry-related metal sites; inconsistent
  counts often reveal a missing or duplicated fragment.

Use covalent-radius sums as a first-pass bond heuristic, then verify expected
functional groups and coordination chemistry explicitly.

## 4. Output contract

Write a conservative P1 CIF with:

- cell lengths and angles;
- `_symmetry_space_group_name_H-M 'P 1'` and International Tables number 1;
- one identity symmetry operation `x,y,z`;
- unique atom labels, element symbols, fractional coordinates, and occupancies;
- all fractional coordinates wrapped into `[0, 1)` unless an unwrapped display
  artifact was explicitly requested.

Keep provenance and source-operation mappings in JSON or Markdown audit output.
Avoid nonstandard CIF tags when broad viewer compatibility is important.

## 5. Report the result

Include:

- source and output paths;
- rows removed or merged and why;
- symmetry operations or subgroup used;
- input/output atom counts and compositions;
- minimum periodic distance and suspicious-contact count;
- representative bond and metal-coordination ranges;
- Pymatgen and ASE parser results;
- Spglib recovered symmetry at the stated tolerance;
- limitations, especially disorder-ordering choices or missing occupancy data.

