#!/usr/bin/env python3
"""Assign PACMOF charges without overwriting source CIFs.

PACMOF installations expose the same ``get_charges_single_serial`` API but are
not always installed into the active Python environment.  Supply either
``--pacmof-root`` (a directory containing ``pacmof.py``) or a normal importable
installation.  The script verifies that each generated CIF has an explicit
charge loop and is neutral to numerical tolerance.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from flexible_mof_gateopen.cif_audit import audit_cif
except ModuleNotFoundError:  # standalone copy inside a skill directory
    from cif_audit import audit_cif  # type: ignore[no-redef]


def load_module(root: Path | None):
    # PACMOF releases before the pymatgen namespace move import CrystalNN from
    # pymatgen.analysis.local_env.  Keep this compatibility alias local to the
    # wrapper; do not modify the user's installed pymatgen package.
    try:
        sys.modules.setdefault("pymatgen.analysis.local_env", importlib.import_module("pymatgen.core.local_env"))
    except ImportError:
        pass
    if root is None:
        return importlib.import_module("pacmof")
    root = root.resolve()
    sys.path.insert(0, str(root))
    module_path = root / "pacmof.py"
    if not module_path.is_file():
        raise FileNotFoundError(module_path)
    spec = importlib.util.spec_from_file_location("pacmof_local", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load PACMOF from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pacmof-root", type=Path, default=None)
    parser.add_argument("--suffix", default="_pacmof_charged")
    args = parser.parse_args()
    module = load_module(args.pacmof_root)
    if not hasattr(module, "get_charges_single_serial"):
        raise SystemExit("PACMOF module lacks get_charges_single_serial")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for input_path in args.inputs:
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
        atoms = module.get_charges_single_serial(
            str(input_path),
            create_cif=True,
            path_to_output_dir=str(args.output_dir),
            add_string=args.suffix,
        )
        generated = args.output_dir / f"{input_path.stem}{args.suffix}.cif"
        if not generated.is_file():
            candidates = sorted(args.output_dir.glob("*.cif"))
            if len(candidates) != 1:
                raise FileNotFoundError(f"could not resolve PACMOF output for {input_path}")
            generated = candidates[0]
        audit = audit_cif(generated)
        if not audit["charge_tag_present"] or not audit["charge_neutral"]:
            raise RuntimeError(f"charged CIF failed neutrality audit: {generated}: {audit}")
        total = sum(float(value) for value in atoms.info.get("_atom_site_charge", []))
        print(f"{generated}: atoms={audit['atom_count']} charge_sum={total:.12g}")


if __name__ == "__main__":
    main()
