#!/usr/bin/env python3
"""Audit a CIF with Gemmi and report formula, cell, occupancy, and charges."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from flexible_mof_gateopen.cif_audit import audit_cif
except ModuleNotFoundError:  # standalone copy inside a skill directory
    from cif_audit import audit_cif  # type: ignore[no-redef]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit_cif(args.cif), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
