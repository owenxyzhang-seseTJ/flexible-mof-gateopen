# End-to-end workflow / 完整工作流

## Entry points / 入口

```bash
python -m pip install -e '.[structure,guest,advanced,mace,xtb,ui,dev]'
python scripts/run_workflow.py --spec workflow.json --output-dir output/my_run
python -m streamlit run scripts/workflow_app.py --server.address 127.0.0.1 --server.port 8513
```

The local-only UI imports two CIFs, a specification and associated JSON files,
then executes the same CLI. Do not expose it to the public network. PACMOF,
trusted host ML weights and native RASPA3 must already be installed. The locally
tested PACMOF RF-DDEC model requires scikit-learn 1.8.0; mismatched serialized
estimator versions are rejected, not ignored. GFN2-xTB is run in a separate
process to avoid mixed torch/tblite OpenMP runtimes on macOS.

界面上传两个 CIF、配置 JSON 和配套文件，保存在新建本地目录。通用结构仍需要
clean-cif 审核记录及源/目标审计 JSON；未知无序、缺氢、溶剂和反离子不会自动放行。
内置自动化学修复目前仅适用于下述公开 MIL-53(Al) 无氢超晶胞。其他材料的结构
清理应由 clean-cif 技能完成后接入；不能用全占有率或无差别删除分子代替审核。

## Stages / 阶段

1. Preserve uploaded originals; inspect raw loops; resolve ordered activation
   model; run independent Pymatgen/ASE/Spglib and chemistry checks.
2. Optional fixed-cell FIRE relaxation, followed by chemistry/connectivity
   re-audit. Automatic re-certification currently supports only MIL-53(Al).
   For other MOFs supply relaxed, reviewed CIFs and set `host.relax=false`.
3. Actual PACMOF inference on final geometry; insert **only** the charge column,
   preserving coordinates, occupancies and cell precision; record code/weight
   hashes, version, neutrality and geometry invariance.
4. RDKit ETKDGv3 conformers → MMFF94s/UFF optimization → optional GFN2-xTB
   refinement → independent bond checks → explicit RASPA atom-type mapping.
   CO₂ MMFF may converge to an incorrect 1.405 Å C–O bond; this now fails the
   geometry gate. Do not equate optimizer convergence with chemical correctness.
5. Fixed-host RASPA at actual reservoir fugacities; automatic supercell sizing,
   cutoff checks, independent seeds, absolute loadings, move diagnostics,
   host–guest energy decomposition and Widom Henry diagnostics.
6. Same host ML model/head for both phases: static energy and finite-displacement
   Phonopy free energy/entropy. Substantial imaginary modes block harmonic ΔF.
7. Descriptive Langmuir fit of N(f) and U_host–guest=aN+bN² fit of actual energy;
   numerical integration of **raw N(f)** plus host ΔF and PΔV for all equilibrium
   crossings. Only closed→open stability changes are opening branches.

每一阶段保留 `attempt_NNN/stage.log`、失败时的 `failure.json`、成功时包含文件哈希的
`completed.json`。相同命令重跑可验证并复用已完成阶段，失败阶段创建新尝试；更改输入
或代码需新建运行目录。恢复是阶段级，不是 RASPA 二进制/优化器步级重启。
最终输出 `workflow_summary.json`；退出码 2 表示没有放行材料压力预测。

## Specification / 配置

Paths are relative to the specification; formula_units is per uploaded CIF
cell, not the whole RASPA supercell. Below is a schema example, not real input.

```json
{
  "temperature_K": 298.15,
  "pressure_bar": [0.001, 0.003, 0.01, 0.03, 0.1, 1.0],
  "reservoir": "coolprop",
  "phases": {
    "closed": {"cif": "closed.cif", "cleaning_report": "closed.cleaning.json", "formula_units": 4},
    "open": {"cif": "open.cif", "cleaning_report": "open.cleaning.json", "formula_units": 4}
  },
  "guest": {
    "name": "CO2", "smiles": "O=C=O", "refinement": "gfn2-xtb",
    "source": "Actual guest force-field reference; changed rigid geometry requires validation"
  },
  "force_field": "force_field.json",
  "force_field_source": "Actual host/guest LJ, charge and mixing references",
  "pacmof_root": "/path/to/trusted/pacmof",
  "host": {
    "relax": false,
    "calculator": {"backend": "mace", "model": "/path/to/mofs_v2.model", "head": "pbe_d3", "dispersion": false},
    "phonons": {"supercell": [2, 2, 2], "mesh": [8, 8, 8], "displacement": 0.01, "imaginary_tolerance_THz": 0.0001, "force_tolerance": 0.03}
  },
  "raspa": {
    "binary": "/path/to/raspa3", "cutoff_A": 12,
    "cycles": 20000, "initialization_cycles": 3000, "equilibration_cycles": 5000,
    "seeds": [417, 941], "point_timeout_seconds": 1800
  },
  "henry_review": {}
}
```

For VOCs add `gas_profile` with sourced Psat(T) or valid Antoine coefficients;
see [vapor_and_guest.md](vapor_and_guest.md). Generation and analysis reject
P>Psat. CoolProp also checks subcritical CO₂ vapor validity. Custom guests need
`guest.molecule_json` and `guest.atom_types` in **RDKit atom order including H**.
The geometry bridge supports rigid guests only. MMFF charges/LJ are not silently
transferred into RASPA. The built-in CO₂ LJ/charges combined with refined gas
geometry remains an unvalidated transfer, not a new validated adsorption model.

`henry_review` takes per-phase evidence strings from actual low-fugacity
refinement/repeat checks. The first two N/f values must also agree within 25%.
This loose numeric check alone does not prove convergence. Zero uptake without
accepted swaps triggers an audit, not an automatic physical failure. An explicit
`adsorption_models.closed` zero/residual scenario can replace the closed sampling
and Henry gate; the open-state and host gates remain active. Positive measured
closed loading contradicts and blocks the zero model. Widom values are diagnostic,
not silently substituted for GCMC or treated as rigorous upper bounds.
See [zero/near-zero closed states](nonadsorbing_closed_phase.zh-CN.md) for schema,
independent empty-host Widom runs, sensitivity and the actual MIL-53 result.

## Public MIL-53 demonstration / 公开示例

Sources: <https://github.com/hnbabaei/Flexible-MOFs>, commit
`d8c445037b3f4be1a60e86fd2caf2ef71dcaff87`.
`MIL-53(Al)_collapsed_4-4-4.cif` and `...expanded...` contain
Al256C2048O1280, **no H**. The adapter verifies all 64 translations, reduces to
Al4C32O20, adds 16 aromatic H and 4 bridging OH H, and requires Al4C32H20O20,
six-coordinate Al and one intended host for each H. No unrelated solvent or
ligand deletion is performed. This is an ordered activated **computational
model**, not a new crystal refinement or proof of experimental phase stability.

```bash
python scripts/make_mil53_demo_spec.py \
  --closed-source /path/to/closed_supercell.cif \
  --open-source /path/to/open_supercell.cif \
  --model /path/to/mofs_v2.model \
  --pacmof-root /path/to/pacmof \
  --raspa-bin /path/to/raspa3 \
  --clean-cif-auditor /path/to/clean-cif/scripts/audit_cif.py \
  --output-dir output/mil53_config --cycles 200
python scripts/run_workflow.py --spec output/mil53_config/workflow.json --output-dir output/mil53_run
```

The small demo uses 268 K, 0.001–20 bar, two seeds, 8 Å cutoff, 1×1×1 phonon
supercell and 2³ q mesh. It is a **real-backend smoke test**, not a converged
material prediction. CO₂ Psat(268 K)≈30.33 bar from CoolProp. Increase cutoff,
cycles, supercell and mesh and review low-pressure sampling for quantitative use.

MACE-MOF0-v2 weights: <https://github.com/ddmms/data/tree/main/mace-mof-0/v2>;
Elena et al., [10.1038/s41524-025-01611-8](https://doi.org/10.1038/s41524-025-01611-8),
CC BY 4.0. Select `pbe_d3`; additional D3 is rejected to avoid double counting.
Successful loading of weights does not prove phase ΔF accuracy.

## Energy definitions / 各种能量的区别

Use **open minus closed**, per identical framework formula unit:

\[
\Delta F_{host}^{harm}=\Delta E_{static}+\Delta U_{vib}-T\Delta S_{vib}
=\Delta E_{static}+\Delta F_{vib}.
\]

不能省略 ΔU_vib（包括零点能、热激发），也不能在声子失败时用静态 ΔE 冒充 ΔF。

\[
\Delta\Omega(P,T)=\Delta F_{host}(T)+P\Delta V
-RT\int_0^{f(P,T)}(n_{open}-n_{closed})\,d\ln f.
\]

Roots are `P_eq`, not experimental `P_open/P_close`. Missing roots within the
sampled range are a legitimate outcome; do not extrapolate or fit the host term
to manufacture one. Barriers, kinetics and hysteresis are not determined by
two rigid structures.

RASPA `molecules/cell` and energies in K refer to the whole simulation box;
`molecules/uc` is the CIF cell. The parser cross-checks them, divides extensive
values by Z×supercell_product, and converts U/kB with R. Native confidence limits
and between-seed spread are retained, not relabeled as standard deviations or a
validated pressure confidence interval. A single-T Langmuir b is an affinity;
it does not uniquely identify LJ interaction energy. The fitted ensemble energy
is not adsorption free energy or isosteric heat.

See [free_energy_methods.zh-CN.md](free_energy_methods.zh-CN.md) and
[framework_entropy_plan.zh-CN.md](framework_entropy_plan.zh-CN.md) for detailed
host thermodynamics and limitations.
