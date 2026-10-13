# flexible-mof-gateopen | Agent execution contract

[中文](README.zh-CN.md)

Purpose: fast, auditable **conditional equilibrium opening pressure P_eq**.
No web application, MD, DFT or finite-displacement phonons in the packaged route.

`closed/open CIF + clean-cif reviews + gas/T + force field + pressure grid`
→ `source charges/PACMOF + RDKit/sourced guest geometry`
→ `static MACE + elastic-network entropy proxy + two RASPA isotherms`
→ `osmotic integration + PΔV → P_eq and sensitivity`.

## Install

Install the screening library from PyPI:

```bash
python -m pip install flexible-mof-gateopen
```

Core has no third-party dependencies. Install only the extras a route needs, for
example `flexible-mof-gateopen[structure,advanced]` or
`flexible-mof-gateopen[ml]`. Every extra is declared in
[pyproject.toml](pyproject.toml).

The workflow **scripts** and the Codex skill bundle are not shipped in the
wheel. For the full workflow, clone the repository and run the commands in the
next section from the checkout:

```bash
git clone https://github.com/owenxyzhang-seseTJ/flexible-mof-gateopen.git
cd flexible-mof-gateopen
python -m pip install -r requirements.txt
```

## Execute

Run inside the skill directory with a configured Python environment:

```bash
python -m pip install -r requirements.txt
python scripts/verify_bundle.py
python scripts/self_test.py
python scripts/run_fast_screen.py --spec /absolute/job.json --output-dir /absolute/new_run
```

Fill `assets/workflow.template.json` using [agent_contract.md](references/agent_contract.md).
Paths in a job are relative to that JSON; bundled model defaults resolve from
the skill root. Override `raspa.binary` on a different platform. Never run an
unfilled template or fabricate chemistry evidence. The source repository's
curated skill is under `skills/flexible-mof-gateopen`.

| Need | Script |
|---|---|
| Full workflow | `run_fast_screen.py` |
| Host ΔE + entropy/free-energy proxy | `quick_host_free_energy.py` |
| Clean-CIF gate | `check_cleaning.py` |
| Charge assignment | `run_pacmof.py` |
| Guest geometry | `prepare_guest_rdkit.py` |
| GCMC input/run | `prepare_raspa_inputs.py` |
| Empty-host dilute uptake audit | `run_empty_widom.py` |
| Bundle/runtime checks | `verify_bundle.py` |
| Analytic regressions without material simulation | `self_test.py` |

## Mandatory decisions

- Use bundled clean-cif instructions/auditor. Chemistry decisions remain reviewed;
  do not blindly remove disorder/solvent. Phases need identical composition and
  activation per formula unit; Z refers to the CIF cell, not the simulation box.
- PACMOF RF-DDEC requires **scikit-learn 1.8.0** for the bundled serialized model.
  Source charges require neutrality, valid atom mapping and provenance.
- RDKit geometry does not provide adsorption LJ/charges. Custom rigid guests need
  explicit atom-type mapping, molecule JSON and model/geometry provenance.
- VOCs require sourced Psat(T); never exceed it. GCMC and analysis must use the
  same actual fugacity. Do not extrapolate above simulated pressure.
- A closed zero-uptake scenario is legitimate but explicit, with residual-uptake
  sensitivity. No accepted swaps is not exclusion; Widom is dilute evidence,
  not a rigorous finite-pressure upper bound.

## Outputs and status

`run/workflow_summary.json` contains status/failures/analysis; stages retain
immutable attempts, hashes and native outputs. Exit 0 means a conditional screen
was produced; exit 2 means unresolved execution/scientific gates. A completed run
with no opening roots means no predicted equilibrium opening **within the sampled
range**, not at every pressure. Dense evaluation nodes do not replace GCMC.

`ΔF_screen=ΔE_static+ΔU_vib,proxy−TΔS_proxy` combines **different models**.
The entropy proxy has no calibrated MOF free-energy accuracy. Scenario ranges
are not confidence intervals; do not infer hysteresis, barriers or experimental
accuracy. Report calibration and independent validation separately.

## Included / external

Included: required project source; MACE-MOF0-v2 and PACMOF RF-DDEC weights;
PACMOF source; clean-cif instructions/auditor; RASPA 3.1.0 and MACE runtime source
archives; macOS arm64 RASPA executable/share tree; pinned requirements/runtime
record; licenses; model provenance; manifest; JSON template; detailed theory PDF
and Markdown. No development trial records, training datasets or application UI.

Python packages and OS shared libraries remain dependencies, not an offline OS
distribution. The native binary may need the recorded native libraries; build the
supplied source elsewhere. No web service is needed after installing dependencies.
Review third-party licenses before redistribution. Only load trusted pickle/Torch
files. The [Theory PDF](assets/principles.pdf) explains derivations and limitations;
its reproducible [XeLaTeX source](references/principles.tex) is included for
text/formula inspection and rebuilding.
