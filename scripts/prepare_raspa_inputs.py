#!/usr/bin/env python3
"""Write portable RASPA3 fixed-host CO2 input directories.

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
    args = parser.parse_args()
    if not args.cif.is_file() or not args.force_field.is_file():
        raise SystemExit("CIF and force-field files must exist")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for pressure in args.pressures:
        point = args.output_dir / f"P_{pressure_tag(pressure)}bar"
        point.mkdir(parents=True, exist_ok=True)
        cif_name = f"{args.phase}.cif"
        shutil.copy2(args.cif, point / cif_name)
        shutil.copy2(args.force_field, point / "force_field.json")
        (point / "CO2.json").write_text(json.dumps(CO2_JSON, indent=2) + "\n", encoding="utf-8")
        (point / "simulation.json").write_text(json.dumps(simulation(args.temperature, pressure, args.phase), indent=2) + "\n", encoding="utf-8")
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
