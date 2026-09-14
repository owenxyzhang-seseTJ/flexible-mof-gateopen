#!/usr/bin/env python3
"""Build a curated source-and-weights agent skill; never include UI/trial data."""
import argparse
from datetime import datetime,timezone
import hashlib
from importlib.metadata import distributions,version
import json
from pathlib import Path
import platform
import shutil
import sys
import tarfile

ROOT=Path(__file__).resolve().parents[1]
HELPERS=["run_fast_screen.py","run_workflow.py","quick_host_free_energy.py",
    "prepare_guest_rdkit.py","run_pacmof.py","prepare_raspa_inputs.py","parse_raspa_output.py",
    "run_empty_widom.py","check_cleaning.py","audit_cif.py","fetch_cif.py",
    "verify_bundle.py","build_principles_pdf.py","self_test.py"]
MODULES=["preflight.py","cif_audit.py","reservoir.py","advanced_thermo.py","host_models.py",
         "elastic_network.py","adsorption_fit.py"]
REFERENCES=["agent_contract.md","principles.zh-CN.md","principles.tex","cleaning_gate.md","vapor_and_guest.md"]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(1024*1024),b""):h.update(data)
    return h.hexdigest()


def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")


def manifest(root):
    files={str(p.relative_to(root)):dict(sha256=sha(p),bytes=p.stat().st_size)
        for p in sorted(root.rglob("*")) if p.is_file() and p.name!="MANIFEST.json" and "__pycache__" not in p.parts}
    dump(root/"MANIFEST.json",dict(schema_version=1,package="flexible-mof-gateopen",version="1.0.0-fast",
        files=files,total_bytes=sum(x["bytes"] for x in files.values())))


def sync_text(root):
    def copy(source,target):
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    for name in ("SKILL.md","README.md","README.zh-CN.md"):
        copy(ROOT/"packaging"/name,root/name)
    for name in HELPERS:copy(ROOT/"scripts"/name,root/"scripts"/name)
    for name in MODULES:copy(ROOT/"src/flexible_mof_gateopen"/name,root/"scripts"/name)
    for name in REFERENCES:copy(ROOT/"references"/name,root/"references"/name)
    copy(ROOT/"LICENSE",root/"LICENSE")
    copy(ROOT/"output/pdf/principles.pdf",root/"assets/principles.pdf")
    metadata=ROOT/"skills/flexible-mof-gateopen/agents/openai.yaml"
    if metadata.is_file() and metadata.resolve()!=(root/"agents/openai.yaml").resolve():
        copy(metadata,root/"agents/openai.yaml")


def install(bundle,targets,backup):
    """Replace only explicitly named skill directories, with recoverable backups."""
    from verify_bundle import verify
    bundle=bundle.resolve();backup=backup.resolve()
    if not verify(bundle)["integrity_passed"]:raise ValueError("cannot install corrupt bundle")
    backup.mkdir(parents=True,exist_ok=False)
    for i,target in enumerate(targets):
        target=target.resolve()
        if target.name!="flexible-mof-gateopen" or target==bundle or target in bundle.parents:
            raise ValueError("exact, separate flexible-mof-gateopen target required")
        if target.exists():shutil.move(str(target),str(backup/f"skill_{i}"))
        shutil.copytree(bundle,target,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
    for name in ("README.md","README.zh-CN.md"):
        shutil.copy2(ROOT/name,backup/name)
        shutil.copy2(ROOT/"packaging"/name,ROOT/name)


def build(out,pacmof,clean,raspa,weights,sources):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    sync_text(out)
    for name in ("pacmof.py","__init__.py","_version.py","LICENSE.md","data/Model_RF_DDEC.pkl"):
        target=out/"vendor/pacmof"/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(pacmof/name,target)
    for name in ("SKILL.md","references/cleaning-decisions.md","scripts/audit_cif.py"):
        target=out/"vendor/clean-cif"/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(clean/name,target)
    shutil.copytree(raspa,out/"vendor/raspa3",ignore=shutil.ignore_patterns(".DS_Store","__pycache__"))
    target=out/"assets/models/mofs_v2.model";target.parent.mkdir(parents=True);shutil.copy2(weights,target)
    for name in ("raspa3-3.1.0.tar.gz","mace_torch-0.3.16.tar.gz"):
        target=out/"vendor/sources"/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(sources/name,target)
    for name in ("MACE_MODEL_LICENSE.txt","MACE_MODEL_LICENSE_AND_CITATION.md"):
        target=out/"licenses"/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(sources/name,target)
    with tarfile.open(sources/"raspa3-3.1.0.tar.gz") as archive:
        member=next(m for m in archive.getmembers() if m.name.endswith("/LICENSE.md"))
        (out/"licenses/RASPA3_LICENSE.md").write_bytes(archive.extractfile(member).read())
    # Runtime libraries are not redistributed; record exact reproducibility surface.
    selected=["numpy","scipy","ase","pymatgen","gemmi","spglib","rdkit","mace-torch","torch",
              "scikit-learn","CoolProp","joblib","pandas","tqdm","dask","distributed"]
    (out/"requirements.txt").write_text("# Tested fast-screen runtime; Python 3.14.4 on macOS arm64\n"+"\n".join(f"{n}=={version(n)}" for n in selected)+"\n")
    dump(out/"runtime-lock.json",dict(python=sys.version,platform=platform.platform(),
        packages={d.metadata["Name"]:d.version for d in distributions() if d.metadata["Name"]},
        native_raspa=dict(version="3.1.0",platform="macOS arm64",external_libraries=["BLIS","LAPACK","zlib","OpenMP","libc++","OpenCL"],
            portability="binary has a local conda-env RPATH; build supplied source or provide compatible native runtime on another machine")))
    dump(out/"MODEL_PROVENANCE.json",dict(
        MACE=dict(path="assets/models/mofs_v2.model",sha256=sha(weights),source="https://github.com/ddmms/data/tree/main/mace-mof-0/v2",
            doi="10.1038/s41524-025-01611-8",license="CC-BY-4.0",head="pbe_d3",extra_D3=False,role="static host energy; not entropy"),
        PACMOF=dict(path="vendor/pacmof/data/Model_RF_DDEC.pkl",sha256=sha(pacmof/"data/Model_RF_DDEC.pkl"),
            doi="10.1021/acs.jctc.0c01229",license="vendor/pacmof/LICENSE.md",sklearn="1.8.0",
            provenance="Exact locally supplied RF-DDEC serialized artifact, retained byte-for-byte; not asserted byte-identical to an original publication release",role="charges only"),
        entropy=dict(role="periodic elastic-network proxy",weights="no learned entropy weights; parameters in scripts/elastic_network.py",calibrated=False),
        RASPA=dict(version="3.1.0",commit="3b81db0e6fe7119732d1525d83a572462a1e8c73",source="https://github.com/iRASPA/RASPA3")))
    dump(out/"assets/workflow.template.json",dict(temperature_K=298.15,pressure_bar=[.001,.003,.01,.03,.1,.3,1.,3.,10.],
        phases={n:dict(cif=f"{n}.cif",cleaning_report=f"{n}.cleaning.json",formula_units=None) for n in ("closed","open")},
        guest=dict(name="CO2",smiles="O=C=O",geometry_sdf="CO2_validated.sdf",geometry_source="SUPPLY_SOURCE",source="SUPPLY_ADSORPTION_MODEL_SOURCE"),
        force_field="force_field.json",force_field_source="SUPPLY_FORCE_FIELD_SOURCE",host=dict(method="elastic-network",relax=False,proxy_mesh=[3,3,3]),
        charges=dict(method="pacmof"),henry_review={}))
    manifest(out)
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("out","pacmof-root","clean-cif-root","raspa-root","weights","source-dir"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args();build(a.out,a.pacmof_root,a.clean_cif_root,a.raspa_root,a.weights,a.source_dir)
    print(json.dumps(dict(bundle=str(a.out.resolve()),manifest=str(a.out/"MANIFEST.json")),indent=2))
