"""Optional Gemmi-backed CIF and charge audit helpers."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any


def audit_cif(path: str | Path) -> dict[str, Any]:
    """Return atom counts, cell volume, occupancies, and charge neutrality.

    Gemmi is intentionally imported lazily so the free-energy utilities remain
    usable in a lightweight environment.
    """
    try:
        import gemmi
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError("CIF auditing requires gemmi") from exc
    block = gemmi.cif.read_file(str(path)).sole_block()
    loop = block.find_loop("_atom_site_label").get_loop()
    tags = list(loop.tags)
    index = {tag: i for i, tag in enumerate(tags)}
    rows = [loop.values[i : i + loop.width()] for i in range(0, len(loop.values), loop.width())]
    counts: dict[str, int] = {}
    occupancies: list[float] = []
    charges: list[float] = []
    for row in rows:
        element = row[index.get("_atom_site_type_symbol", index.get("_atom_site_label", 0))]
        counts[element] = counts.get(element, 0) + 1
        if "_atom_site_occupancy" in index:
            occupancies.append(float(row[index["_atom_site_occupancy"]]))
        if "_atom_site_charge" in index:
            charges.append(float(row[index["_atom_site_charge"]]))
    def scalar(tag: str) -> float:
        value = block.find_value(tag)
        if value in {"", "?", "."}:
            raise ValueError(f"missing CIF cell parameter {tag}")
        return float(value)

    a, b, c = (scalar(tag) for tag in ("_cell_length_a", "_cell_length_b", "_cell_length_c"))
    alpha, beta, gamma = (math.radians(scalar(tag)) for tag in ("_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))
    volume = a * b * c * math.sqrt(
        1.0
        - math.cos(alpha) ** 2
        - math.cos(beta) ** 2
        - math.cos(gamma) ** 2
        + 2.0 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
    )
    return {
        "path": str(Path(path)),
        "block": block.name,
        "atom_count": len(rows),
        "element_counts": dict(sorted(counts.items())),
        "cell_volume_A3": volume,
        "max_occupancy_deviation": max((abs(value - 1.0) for value in occupancies), default=None),
        "charge_tag_present": bool(charges),
        "charge_sum_e": sum(charges) if charges else None,
        "charge_neutral": abs(sum(charges)) < 1.0e-8 if charges else None,
    }
