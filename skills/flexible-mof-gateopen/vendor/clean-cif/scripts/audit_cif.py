#!/usr/bin/env python3
"""Read-only CIF audit using raw loops, Pymatgen, ASE, and Spglib."""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

def central_number(value):
    try:
        return float(re.sub(r"\([^)]*\)$", "", str(value)))
    except (TypeError, ValueError):
        return None


def raw_metadata(path: Path):
    result = {"readable": False}
    try:
        from pymatgen.io.cif import CifFile

        cif = CifFile.from_file(path)
        block_name, block = next(iter(cif.data.items()))
        data = block.data
        labels = list(data.get("_atom_site_label", []))
        elements = list(data.get("_atom_site_type_symbol", []))
        xs = list(data.get("_atom_site_fract_x", []))
        ys = list(data.get("_atom_site_fract_y", []))
        zs = list(data.get("_atom_site_fract_z", []))
        coordinates = [
            [central_number(x), central_number(y), central_number(z)]
            for x, y, z in zip(xs, ys, zs)
        ]
        outside = [
            index + 1
            for index, xyz in enumerate(coordinates)
            if None not in xyz and any(value < 0 or value >= 1 for value in xyz)
        ]
        result.update(
            {
                "readable": True,
                "block": block_name,
                "declared_space_group": data.get(
                    "_symmetry_space_group_name_H-M",
                    data.get("_space_group_name_H-M_alt"),
                ),
                "declared_number": data.get(
                    "_symmetry_Int_Tables_number",
                    data.get("_space_group_IT_number"),
                ),
                "symmetry_operation_count": len(
                    data.get(
                        "_symmetry_equiv_pos_as_xyz",
                        data.get("_space_group_symop_operation_xyz", []),
                    )
                ),
                "atom_rows": len(labels),
                "duplicate_labels": sorted(
                    label for label in set(labels) if labels.count(label) > 1
                ),
                "rows_outside_unit_interval": outside,
                "elements_in_raw_loop": sorted(set(elements)),
                "occupancy_field_present": "_atom_site_occupancy" in data,
            }
        )
    except Exception as exc:  # audit must report malformed inputs, not hide them
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def parse_pymatgen(path: Path):
    result = {"success": False}
    try:
        from pymatgen.io.cif import CifParser

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            structure = CifParser(str(path)).parse_structures(primitive=False)[0]
        result.update(
            {
                "success": True,
                "atom_count": len(structure),
                "formula": structure.composition.formula,
                "warnings": [str(item.message) for item in caught],
            }
        )
        return result, structure
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result, None


def parse_ase(path: Path):
    result = {"success": False}
    try:
        from ase.io import read

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            atoms = read(str(path))
        result.update(
            {
                "success": True,
                "atom_count": len(atoms),
                "formula": atoms.get_chemical_formula(),
                "warnings": [str(item.message) for item in caught],
            }
        )
        return result, atoms
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        return result, None


def structure_geometry(structure, shortest_count: int):
    pairs = []
    for i in range(len(structure)):
        for j in range(i + 1, len(structure)):
            pairs.append((float(structure.get_distance(i, j)), i, j))
    pairs.sort()
    shortest = []
    for distance, i, j in pairs[:shortest_count]:
        shortest.append(
            {
                "distance_A": round(distance, 6),
                "site_1": i + 1,
                "element_1": structure[i].specie.symbol,
                "site_2": j + 1,
                "element_2": structure[j].specie.symbol,
            }
        )
    return {
        "minimum_periodic_distance_A": round(pairs[0][0], 6) if pairs else None,
        "pairs_under_0.75_A": sum(distance < 0.75 for distance, _, _ in pairs),
        "pairs_under_1.00_A": sum(distance < 1.00 for distance, _, _ in pairs),
        "shortest_pairs": shortest,
    }


def recovered_symmetry(atoms, symprec: float):
    result = {"success": False, "symprec_A": symprec}
    try:
        import spglib

        dataset = spglib.get_symmetry_dataset(
            (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers),
            symprec=symprec,
        )
        if dataset is None:
            result["error"] = "Spglib returned no symmetry dataset"
        else:
            result.update(
                {
                    "success": True,
                    "international": dataset.international,
                    "number": int(dataset.number),
                    "operation_count": len(dataset.rotations),
                }
            )
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    parser.add_argument("--json", type=Path, help="Write the report to this path")
    parser.add_argument("--shortest", type=int, default=20)
    parser.add_argument("--symprec", type=float, default=0.01)
    args = parser.parse_args()

    path = args.cif.expanduser().resolve()
    report = {"file": str(path), "raw": raw_metadata(path)}
    report["pymatgen"], structure = parse_pymatgen(path)
    report["ase"], atoms = parse_ase(path)
    if structure is not None:
        report["geometry"] = structure_geometry(structure, args.shortest)
    if atoms is not None:
        report["spglib"] = recovered_symmetry(atoms, args.symprec)
    if report["pymatgen"].get("success") and report["ase"].get("success"):
        report["parser_atom_counts_agree"] = (
            report["pymatgen"]["atom_count"] == report["ase"]["atom_count"]
        )

    text = json.dumps(report, indent=2, ensure_ascii=False)
    print(text)
    if args.json:
        args.json.expanduser().resolve().write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
