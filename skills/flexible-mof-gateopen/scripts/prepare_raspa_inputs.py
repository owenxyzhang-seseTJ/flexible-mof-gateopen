#!/usr/bin/env python3
"""Write audited RASPA3 fixed-host single-guest input directories.

The default is write-only.  Add ``--run`` only when the caller explicitly
wants to launch a local RASPA executable.  The framework is rigid; this helper
does not claim to run a flexible-host or kinetic simulation.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
try:
    from flexible_mof_gateopen.preflight import load_gas_profile, pressure_gate, structure_gate
except ModuleNotFoundError:
    from preflight import load_gas_profile, pressure_gate, structure_gate


CO2_JSON = {
    "CriticalTemperature": 304.1282,
    "CriticalPressure": 7377300.0,
    "AcentricFactor": 0.22394,
    "Type": "rigid",
    "pseudoAtoms": [["O_co2", [0.0, 0.0, 1.149]], ["C_co2", [0.0, 0.0, 0.0]], ["O_co2", [0.0, 0.0, -1.149]]],
}


def pressure_tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def simulation(temperature: float, pressure_bar: float, name: str) -> dict:
    return {
        "SimulationType": "MonteCarlo",
        "NumberOfCycles": 20000,
        "NumberOfInitializationCycles": 3000,
        "NumberOfEquilibrationCycles": 5000,
        "PrintEvery": 5000,
        "Systems": [{
            "Type": "Framework",
            "Name": name,
            "HeliumVoidFraction": 0.30,
            "NumberOfUnitCells": [1, 1, 1],
            "ExternalTemperature": temperature,
            "ExternalPressure": pressure_bar * 1.0e5,
            "ChargeMethod": "Ewald",
            "UseChargesFrom": "CIF_File",
        }],
        "Components": [{
            "Name": "CO2",
            "MoleculeDefinition": "ExampleDefinitions",
            "FugacityCoefficient": 1.0,
            "TranslationProbability": 0.5,
            "RotationProbability": 0.5,
            "ReinsertionProbability": 0.5,
            "SwapProbability": 1.0,
            "WidomProbability": 1.0,
            "CreateNumberOfMolecules": 0,
        }],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--phase", default="host")
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--pressures", type=float, nargs="+", required=True)
    parser.add_argument("--force-field", type=Path, required=True)
    parser.add_argument("--raspa-bin", type=Path, default=None)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--cleaning-report", type=Path, required=True)
    parser.add_argument("--gas-profile", type=Path, help="Required for any guest other than bundled CO2")
    parser.add_argument("--molecule-json", type=Path, help="Validated RASPA guest model for a custom gas")
    parser.add_argument("--unit-cells", type=int, nargs=3, required=True)
    parser.add_argument("--cutoff", type=float, required=True, help="Largest vdW/Coulomb real-space cutoff in Angstrom")
    args = parser.parse_args()
    if not args.cif.is_file() or not args.force_field.is_file():
        raise SystemExit("CIF and force-field files must exist")
    profile = load_gas_profile(args.gas_profile)
    pressure_audit = pressure_gate(profile, args.temperature, args.pressures)
    structure_audit = structure_gate(args.cif, args.cleaning_report, charges=True)
    if profile["name"] != "CO2" and args.molecule_json is None:
        raise ValueError("custom adsorbate requires --molecule-json; CO2 coordinates cannot be reused")
    molecule = json.loads(args.molecule_json.read_text()) if args.molecule_json else CO2_JSON
    ff = json.loads(args.force_field.read_text())
    import numpy as np
    from ase.io import read
    pseudo = {p["name"]: p for p in ff.get("PseudoAtoms", [])}
    guest_atoms = molecule.get("pseudoAtoms", [])
    if not guest_atoms or any(a[0] not in pseudo for a in guest_atoms):
        raise ValueError("guest atom types must be present in force_field PseudoAtoms")
    if any(len(a[1]) != 3 or not np.isfinite(np.array(a[1], dtype=float)).all() for a in guest_atoms):
        raise ValueError("guest coordinates must be finite 3D vectors")
    guest_charge = sum(float(pseudo[a[0]]["charge"]) for a in guest_atoms)
    if not np.isfinite(guest_charge) or abs(guest_charge) > 1e-6:
        raise ValueError("single-component vapor guest must have finite neutral charge")
    if any(n <= 0 for n in args.unit_cells) or not np.isfinite(args.cutoff) or args.cutoff <= 0:
        raise ValueError("positive supercell counts and cutoff required")
    lattice = read(str(args.cif)).cell.array * np.array(args.unit_cells)[:, None]
    heights = 1 / np.linalg.norm(np.linalg.inv(lattice), axis=0)
    if min(heights) < 2 * args.cutoff:
        raise ValueError("supercell perpendicular heights must be >= twice the real-space cutoff")
    # Explicitly bind the supplied real-space range rather than relying on defaults.
    ff.update({"CutOffFrameworkVDW": args.cutoff, "CutOffMoleculeVDW": args.cutoff, "CutOffCoulomb": args.cutoff})
    if not args.phase or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in args.phase):
        raise ValueError("phase must be a safe filename stem")
    if args.run and (args.raspa_bin is None or not args.raspa_bin.is_file()):
        raise ValueError("--run requires an existing RASPA executable")
    points = [args.output_dir / f"P_{pressure_tag(p)}bar" for p in args.pressures]
    if len(set(points)) != len(points) or any(p.exists() for p in points):
        raise ValueError("output pressure directories already exist or collide; use a new output directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for pressure in args.pressures:
        point = args.output_dir / f"P_{pressure_tag(pressure)}bar"
        point.mkdir(parents=True, exist_ok=True)
        cif_name = f"{args.phase}.cif"
        shutil.copy2(args.cif, point / cif_name)
        (point / "force_field.json").write_text(json.dumps(ff, indent=2) + "\n")
        (point / f"{profile['name']}.json").write_text(json.dumps(molecule, indent=2) + "\n")
        sim = simulation(args.temperature, pressure, args.phase)
        sim["Systems"][0]["NumberOfUnitCells"] = args.unit_cells
        sim["Components"][0]["Name"] = profile["name"]
        (point / "simulation.json").write_text(json.dumps(sim, indent=2) + "\n", encoding="utf-8")
        (point / "preflight.json").write_text(json.dumps({"structure": structure_audit, "reservoir": pressure_audit}, indent=2) + "\n")
        if args.run:
            if args.raspa_bin is None:
                raise SystemExit("--raspa-bin is required with --run")
            if not args.raspa_bin.is_file():
                raise SystemExit(f"RASPA executable not found: {args.raspa_bin}")
            env = os.environ.copy()
            with (point / "run.log").open("w", encoding="utf-8") as log:
                subprocess.run([str(args.raspa_bin)], cwd=point, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f"wrote {len(args.pressures)} RASPA input directories under {args.output_dir}")


if __name__ == "__main__":
    main()
