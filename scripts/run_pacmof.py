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
import json
import warnings
from importlib.metadata import version

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

try:
    from flexible_mof_gateopen.cif_audit import audit_cif
    from flexible_mof_gateopen.preflight import structure_gate, sha256
except ModuleNotFoundError:  # standalone copy inside a skill directory
    from cif_audit import audit_cif  # type: ignore[no-redef]
    from preflight import structure_gate, sha256


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
    parser.add_argument("--cleaning-reports", nargs="+", type=Path, required=True, help="One reviewed clean-cif record per input, in matching order")
    args = parser.parse_args()
    if len(args.inputs) != len(args.cleaning_reports):
        raise ValueError("one cleaning report per input is required")
    if not args.suffix or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in args.suffix):
        raise ValueError("nonempty safe output suffix required")
    destinations = [args.output_dir / f"{p.stem}{args.suffix}.cif" for p in args.inputs]
    if len(set(destinations)) != len(destinations) or any(p.exists() for p in destinations):
        raise ValueError("PACMOF output exists or collides; use a new directory")
    for path, report in zip(args.inputs, args.cleaning_reports):
        structure_gate(path, report)
    module = load_module(args.pacmof_root)
    if not hasattr(module, "get_charges_single_serial"):
        raise SystemExit("PACMOF module lacks get_charges_single_serial")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for input_path, report_path in zip(args.inputs, args.cleaning_reports):
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
        from sklearn.exceptions import InconsistentVersionWarning
        with warnings.catch_warnings():
            warnings.simplefilter("error", InconsistentVersionWarning)
            atoms = module.get_charges_single_serial(str(input_path), create_cif=False)
        generated = args.output_dir / f"{input_path.stem}{args.suffix}.cif"
        import gemmi
        import numpy as np
        # PACMOF's legacy writer rounds coordinates; preserve the exact input
        # atom table and replace ONLY the inferred charge column instead.
        doc = gemmi.cif.read_file(str(input_path))
        b = doc.sole_block()
        labels = b.find_values("_atom_site_label")
        predicted = np.asarray(atoms.info.get("_atom_site_charge", []), dtype=float)
        if len(predicted) != len(labels) or not np.isfinite(predicted).all() or abs(predicted.sum()) > 1e-6:
            raise ValueError("PACMOF returned missing/nonfinite/nonneutral charges")
        if not b.find_values("_atom_site_charge"):
            b.find_loop("_atom_site_label").get_loop().add_columns(["_atom_site_charge"], "0")
        column = b.find_values("_atom_site_charge")
        for i, value in enumerate(predicted): column[i] = f"{value:.15g}"
        doc.write_file(str(generated))
        before = gemmi.cif.read_file(str(input_path)).sole_block()
        after = gemmi.cif.read_file(str(generated)).sole_block()
        for tag in ("_atom_site_type_symbol",):
            if list(before.find_values(tag)) != list(after.find_values(tag)):
                raise ValueError("PACMOF changed element ordering")
        for tag in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"):
            x = [gemmi.cif.as_number(v) for v in before.find_values(tag)]
            y = [gemmi.cif.as_number(v) for v in after.find_values(tag)]
            if len(x) != len(y) or not np.allclose(x, y, atol=1e-8, rtol=0):
                raise ValueError(f"PACMOF changed {tag}")
        for tag in ("_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"):
            if not np.isclose(gemmi.cif.as_number(before.find_value(tag)), gemmi.cif.as_number(after.find_value(tag)), atol=1e-7, rtol=0):
                raise ValueError("PACMOF changed cell")
        audit = audit_cif(generated)
        if not audit["charge_tag_present"] or not audit["charge_neutral"]:
            raise RuntimeError(f"charged CIF failed neutrality audit: {generated}: {audit}")
        # Preserve the reviewed chemistry and rebind only after the geometry checks.
        report = json.loads(report_path.read_text())
        for field in ("source_audit", "candidate_audit"):
            report[field]["path"] = str((report_path.parent / report[field]["path"]).resolve())
        report["precharge_sha256"] = sha256(input_path)
        model_path = Path(module.__file__).parent / "data/Model_RF_DDEC.pkl"
        report["charge_model"] = dict(method="PACMOF RF-DDEC", sklearn_version=version("scikit-learn"),
                                     model_sha256=sha256(model_path), code_sha256=sha256(module.__file__),
                                     coordinate_policy="unchanged original table; charge-only insertion")
        report["candidate_sha256"] = sha256(generated)
        charged_report = generated.with_suffix(".cleaning.json")
        charged_report.write_text(json.dumps(report, indent=2) + "\n")
        structure_gate(generated, charged_report, charges=True)
        total = sum(float(value) for value in atoms.info.get("_atom_site_charge", []))
        print(f"{generated}: atoms={audit['atom_count']} charge_sum={total:.12g}")


if __name__ == "__main__":
    main()
