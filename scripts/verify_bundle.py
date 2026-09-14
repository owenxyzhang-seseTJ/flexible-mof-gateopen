#!/usr/bin/env python3
"""Verify all bundled artifacts before loading trusted pickle/torch weights."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path


def verify(root):
    manifest=json.loads((root/"MANIFEST.json").read_text())
    failures=[]
    for name,expected in manifest["files"].items():
        p=root/name
        if not p.is_file():failures.append(name+": missing");continue
        h=hashlib.sha256()
        with p.open("rb") as f:
            for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
        if h.hexdigest()!=expected["sha256"]:failures.append(name+": hash mismatch")
    versions={}
    for name in ["numpy","scipy","ase","pymatgen","gemmi","spglib","rdkit","mace-torch","torch","scikit-learn","CoolProp"]:
        try:versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:versions[name]=None
    compatible=versions["scikit-learn"]=="1.8.0" and all(versions.values())
    return dict(integrity_passed=not failures,failures=failures,file_count=len(manifest["files"]),
                runtime_dependencies_present_and_sklearn_compatible=compatible,versions=versions,
                platform=platform.platform(),
                warning="Bundled RASPA executable is macOS arm64-specific; other systems must build supplied source. This check does not certify material accuracy or binary dynamic-library availability.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--root",type=Path,default=Path(__file__).resolve().parent.parent)
    a=p.parse_args();r=verify(a.root);print(json.dumps(r,indent=2));raise SystemExit(0 if r["integrity_passed"] and r["runtime_dependencies_present_and_sklearn_compatible"] else 2)
