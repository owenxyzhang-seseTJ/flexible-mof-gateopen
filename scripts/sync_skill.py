#!/usr/bin/env python3
"""Synchronize portable skill helpers from canonical repository sources."""
import argparse
from pathlib import Path
import shutil
from build_agent_skill import sync_text,manifest

ROOT = Path(__file__).resolve().parents[1]
SKILL = ROOT / "skills/flexible-mof-gateopen"
HELPERS = ["audit_cif.py", "check_cleaning.py", "fetch_cif.py", "parse_raspa_output.py",
           "prepare_raspa_inputs.py", "prepare_guest_rdkit.py", "quick_free_energy_screen.py", "run_pacmof.py", "validate_example.py",
           "advanced_free_energy.py", "validate_advanced.py", "audit_published_mil53_entropy.py",
           "prepare_public_mil53.py", "relax_host.py", "run_workflow.py", "make_mil53_demo_spec.py", "refine_guest_tblite.py", "quick_host_free_energy.py",
           "run_empty_widom.py", "screen_closed_exclusion.py"]


def sync(install_dir=None):
    if (SKILL/"MANIFEST.json").exists():
        # Curated fast bundle: never reintroduce legacy trial/UI/slow-route files.
        sync_text(SKILL)
        manifest(SKILL)
        if install_dir is not None:
            destination=Path(install_dir).expanduser().resolve()
            if destination.name!="flexible-mof-gateopen" or destination==SKILL.resolve():
                raise ValueError("install-dir must be a separate flexible-mof-gateopen directory")
            shutil.copytree(SKILL,destination,dirs_exist_ok=True,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        return
    for name in HELPERS:
        shutil.copy2(ROOT / "scripts" / name, SKILL / "scripts" / name)
    for name in ("free_energy.py", "cif_audit.py", "preflight.py", "reservoir.py", "advanced_thermo.py", "host_models.py", "framework_ti.py", "adsorption_fit.py", "elastic_network.py"):
        shutil.copy2(ROOT / "src/flexible_mof_gateopen" / name, SKILL / "scripts" / name)
    for source in (ROOT / "references").glob("*.md"):
        shutil.copy2(source, SKILL / "references" / source.name)
    if install_dir is not None:
        destination = Path(install_dir).expanduser().resolve()
        if destination.name != "flexible-mof-gateopen" or destination == SKILL.resolve():
            raise ValueError("install-dir must be a separate flexible-mof-gateopen directory")
        shutil.copytree(SKILL, destination, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--install-dir", type=Path)
    args = parser.parse_args()
    sync(args.install_dir)
