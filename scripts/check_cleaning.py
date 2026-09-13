#!/usr/bin/env python3
"""Validate a reviewed clean-cif candidate before charging or GCMC."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
try:
    from flexible_mof_gateopen.preflight import structure_gate
except ModuleNotFoundError:
    from preflight import structure_gate

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cif", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--charges", action="store_true")
    a = p.parse_args()
    print(json.dumps(structure_gate(a.cif, a.report, charges=a.charges), indent=2))
