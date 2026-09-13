"""Fail-closed structure and single-component vapor preparation gates."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path


def positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_gas_profile(path=None):
    if path is None:
        return {"name": "CO2", "kind": "non-VOC", "source": "bundled CO2 model"}
    profile = json.loads(Path(path).read_text())
    if profile.get("kind") not in {"VOC", "non-VOC"} or not profile.get("source"):
        raise ValueError("gas profile requires kind=VOC/non-VOC and source")
    name = profile.get("name", "")
    if not name or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in name):
        raise ValueError("gas name must be a safe molecule filename stem")
    return profile


def pressure_gate(profile, temperature, pressures):
    """Return provenance and P/Psat; never extrapolate a vapor-pressure fit."""
    temperature = positive(temperature, "temperature_K")
    pressures = [positive(p, "pressure_bar") for p in pressures]
    if len(set(pressures)) != len(pressures):
        raise ValueError("duplicate pressures")
    saturation = profile.get("saturation")
    if profile["kind"] == "VOC" and not saturation:
        raise ValueError("VOC requires sourced saturation pressure at this temperature")
    if not saturation:
        return {"gas": profile["name"], "temperature_K": temperature, "pressure_bar": pressures, "psat_bar": None}
    if not saturation.get("source"):
        raise ValueError("saturation.source is required")
    if "critical_temperature_K" in profile and temperature >= positive(profile["critical_temperature_K"], "critical temperature"):
        raise ValueError("saturation-limited workflow is undefined at/above critical temperature")
    method = saturation.get("method")
    if method == "value":
        if not math.isclose(temperature, float(saturation["temperature_K"]), rel_tol=0, abs_tol=1e-6):
            raise ValueError("Psat temperature does not match simulation temperature")
        psat = positive(saturation["pressure_bar"], "Psat")
    elif method == "antoine":
        lo, hi = map(float, saturation["valid_temperature_K"])
        if not (0 < lo <= temperature <= hi and math.isfinite(hi)):
            raise ValueError("temperature outside Antoine validity interval")
        if saturation.get("temperature_unit") not in {"K", "C"}:
            raise ValueError("Antoine temperature_unit must be K or C")
        scales = {"bar": 1., "Pa": 1e-5, "kPa": .01, "mmHg": .001333223684}
        if saturation.get("pressure_unit") not in scales:
            raise ValueError("unsupported Antoine pressure unit")
        a, b, c = (float(saturation[k]) for k in ("A", "B", "C"))
        t = temperature if saturation["temperature_unit"] == "K" else temperature - 273.15
        psat = positive(10 ** (a - b / (t + c)) * scales[saturation["pressure_unit"]], "Psat")
    else:
        raise ValueError("saturation method must be value or antoine")
    if any(p > psat * (1 + 1e-12) for p in pressures):
        raise ValueError(f"pressure exceeds {profile['name']} Psat({temperature:g} K)={psat:.9g} bar")
    return {"gas": profile["name"], "temperature_K": temperature, "pressure_bar": pressures,
            "psat_bar": psat, "relative_pressure": [p / psat for p in pressures], "saturation": saturation}


def structure_gate(cif, report_path, *, charges=False):
    """Bind a chemistry-reviewed clean-cif record to the exact candidate file.

    Chemical judgments (solvent identity, counterions, ordered configuration)
    are explicit audit inputs; numeric checks do not manufacture them.
    """
    import gemmi
    import numpy as np
    import spglib
    from ase.io import read
    from pymatgen.io.cif import CifParser

    cif, report_path = Path(cif), Path(report_path)
    report = json.loads(report_path.read_text())
    if report.get("schema_version") != 1 or report.get("status") != "passed":
        raise ValueError("a passed schema_version=1 clean-cif report is required")
    if report.get("candidate_sha256") != sha256(cif):
        raise ValueError("clean-cif report hash does not match candidate CIF")
    source_digest = report.get("source_sha256", "")
    if len(source_digest) != 64 or any(c not in "0123456789abcdef" for c in source_digest) or not report.get("transformation_record"):
        raise ValueError("source hash and transformation record are required")
    for name in ("disorder", "solvent"):
        record = report.get(name, {})
        choices = {"absent", "ordered"} if name == "disorder" else {"absent", "removed", "retained_with_reason"}
        if record.get("status") not in choices or not record.get("reason"):
            raise ValueError(f"unresolved {name}: explicit disposition and reason required")
    for check in ("connectivity", "hydrogen_hosts", "metal_coordination", "nonbonded_contacts", "composition_charge_balance"):
        record = report.get("chemistry", {}).get(check, {})
        if record.get("status") not in {"passed", "not_applicable"} or not record.get("evidence"):
            raise ValueError(f"chemistry review missing: {check}")
    for check in ("source_audit", "candidate_audit"):
        record = report.get(check, {})
        if not record.get("path") or not record.get("sha256"):
            raise ValueError(f"{check} artifact required")
        path = report_path.parent / record["path"]
        if sha256(path) != record["sha256"]:
            raise ValueError(f"{check} artifact hash mismatch")

    block = gemmi.cif.read_file(str(cif)).sole_block()
    labels = list(block.find_values("_atom_site_label"))
    if not labels or len(set(labels)) != len(labels):
        raise ValueError("missing or duplicated labels")
    for axis in "xyz":
        coords = [gemmi.cif.as_number(x) for x in block.find_values(f"_atom_site_fract_{axis}")]
        if len(coords) != len(labels) or any(not math.isfinite(x) or not 0 <= x < 1 for x in coords):
            raise ValueError("coordinates must be finite and wrapped into [0,1)")
    occ = [gemmi.cif.as_number(x) for x in block.find_values("_atom_site_occupancy")]
    if len(occ) != len(labels) or any(not math.isfinite(x) or abs(x-1) > 1e-8 for x in occ):
        raise ValueError("unresolved or missing occupancies")
    for tag in ("_atom_site_disorder_group", "_atom_site_disorder_assembly"):
        if any(x not in {".", "?", "0"} for x in block.find_values(tag)):
            raise ValueError("unresolved disorder fields")
    operations = []
    for tag in ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"):
        operations.extend(block.find_values(tag))
    if not operations or any(gemmi.cif.as_string(x).replace(" ", "").lower() != "x,y,z" for x in operations):
        raise ValueError("identity-only P1 symmetry required")
    names = [block.find_value(t) for t in ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt")]
    if not any(names) or any(gemmi.cif.as_string(x).replace(" ", "").upper() != "P1" for x in names if x):
        raise ValueError("P1 declaration required")
    for tag in ("_symmetry_Int_Tables_number", "_space_group_IT_number"):
        if block.find_value(tag) and gemmi.cif.as_number(block.find_value(tag)) != 1:
            raise ValueError("conflicting space-group number")
    pm = CifParser(str(cif)).parse_structures(primitive=False)[0]
    ase = read(str(cif))
    if not pm.is_ordered or len(pm) != len(labels) or len(ase) != len(pm) or dict(pm.composition.get_el_amt_dict()) != dict(Counter(ase.get_chemical_symbols())):
        raise ValueError("Pymatgen/ASE/raw atom counts or compositions disagree")
    if not np.allclose(pm.lattice.matrix @ pm.lattice.matrix.T, ase.cell.array @ ase.cell.array.T, atol=1e-5):
        raise ValueError("parser cell metrics disagree")
    _, _, _, distances = pm.get_neighbor_list(0.75)
    if len(distances):
        raise ValueError("sub-0.75 A periodic contacts require resolution")
    symmetry = spglib.get_symmetry_dataset((ase.cell.array, ase.get_scaled_positions(), ase.numbers), symprec=.01)
    if symmetry is None:
        raise ValueError("Spglib cross-check failed")
    if charges:
        values = [gemmi.cif.as_number(x) for x in block.find_values("_atom_site_charge")]
        if len(values) != len(labels) or any(not math.isfinite(x) for x in values) or abs(sum(values)) > 1e-6:
            raise ValueError("finite neutral CIF charges required")
    return {"status": "passed", "candidate_sha256": sha256(cif), "atom_count": len(pm),
            "composition": pm.composition.get_el_amt_dict(), "recovered_space_group": symmetry.international,
            "cleaning_report_sha256": sha256(report_path)}
