# Advanced free energy workflows / 高级自由能工作流

This is an executable screening toolkit, not a validated gate-pressure predictor for every MOF.
本文区分**算法实现、真实软件接口验证、具体材料定量验证**。前两项已经执行；本次没有获得 MIL-53(Al) 或 X-pcu-5-Zn 独立收敛的框架自由能。

## 1. 安装与统一入口

建议使用单独的 Python 3.11/3.12 环境，不改变现有 RASPA 环境：

```bash
python -m pip install -e '.[structure,advanced,ml,mace,dev]'
python scripts/advanced_free_energy.py --help
python scripts/validate_advanced.py --ml --out output/software-checks.json
```

`advanced` 安装 NumPy、SciPy、CoolProp、ASE、Phonopy；`ml` 安装 CHGNet；`mace` 安装 MACE。首次运行可能下载官方预训练权重；只加载可信来源的模型，部分模型文件使用 Python 序列化，不可加载陌生权重。权重不随本仓库分发，许可证以原作者规定为准。需要 MACE D3 时另装 `torch-dftd`，并明确模型是否已经包含色散，禁止重复添加。

所有命令采用：

```bash
python scripts/advanced_free_energy.py COMMAND --spec specification.json --out result.json
```

支持 `host-energy`、`phonons`、`qha`、`mc-ti`、`reservoir`、`phase-diagram`、`landscape`。结果不覆盖已有文件；spec 内路径相对于 spec 本身。结果附输入 SHA256。ML/声子/框架 MC 的 CIF 必须有通过检查的 clean-cif 报告，不提供跳过选项。库层的 ASE 对象接口用于集成和合成测试，不等于 CIF 化学审核。

## 2. 预训练势与 ΔE

`host-energy` 示例结构：

```json
{
  "calculator": {"backend": "chgnet", "model": "0.3.0", "device": "cpu", "dispersion": false},
  "states": [
    {"name": "closed", "cif": "closed.cif", "cleaning_report": "closed.cleaning.json", "formula_units": 4},
    {"name": "open", "cif": "open.cif", "cleaning_report": "open.cleaning.json", "formula_units": 2}
  ]
}
```

可替换为 `{"backend":"mace","model":"small","device":"cpu","dispersion":false}`，也支持经核实的 MACE 别名或本地 checkpoint 路径。默认小模型用于快速接口检验，不代表推荐其直接用于柔性 MOF 定量预测。CHGNet 的 `model` 是其官方版本名。

计算同组成、同活化状态各相的单点能量，并以第一相为参照：

`ΔE_i = 96.4853321233 × (E_i[eV/cell]/Z_i − E_ref[eV/cell]/Z_ref)`，单位 kJ/mol-f.u.。

输出原始能量、每化学式单元能量、体积、最大原子力、ASE Voigt 应力（eV/Å³）、模型版本及载入张量的哈希。组成/Z 不匹配立即停止。PACMOF 的电荷不作为 CHGNet/MACE 输入；电荷审核不能替代势函数验证。

这是 **ΔE 静态筛选**，不是 ΔF(T)。不自动松弛，因为不受约束的松弛可能让开/闭孔都落入同一势阱。应先在各相保持可解释的约束/拓扑，检查固定体积的内部坐标收敛，再把重审的结构送入声子计算。对于不稳定或强非谐的呼吸模式，不能用静态能差宣称温度相关自由能。

材料校验至少包括：色散处理是否适配、开闭孔 DFT 相对能、力/应力、势阱是否保持、软模稳定性。不同模型的绝对能量不能互相相减；可比较各自预测的相内 ΔE，作为模型敏感性而非校准置信区间。

## 3. Phonopy 谐振自由能

`phonons` 使用上述 calculator/states，并增加：

```json
{
  "phonons": {
    "temperatures": [0, 268, 298.15],
    "supercell": [2, 2, 2],
    "mesh": [8, 8, 8],
    "displacement": 0.01,
    "imaginary_tolerance_THz": 0.1,
    "force_tolerance": 0.03
  }
}
```

流程：确认输入原子力小于阈值 → Phonopy 生成对称不等价有限位移超胞 → ML 势求力 → 去除数值净力漂移 → 生成/对称化力常数 → 网格声子 → `F_vib(T)`。

`F_vib = Σ_modes [ℏω/2 + kBT ln(1 − exp(−ℏω/kBT))]`，包括零点能。

强制使用 identity primitive matrix，使 Phonopy 的 kJ/mol 指向**输入晶胞**，再除以输入的 Z；不是除以位移超胞的 Z。输出 `F_total = E_static + F_vib`、熵、Cv、最低频率、位移数和最大净力漂移。

同次命令包含两相时，以第一相为参考，另输出 `harmonic_state_differences`：静态 ΔE、振动 ΔU、ΔS、−TΔS 和谐振 ΔF。该配对汇总采用更严格的 −10⁻⁴ THz 虚频门槛，不允许把静态 ΔE−TΔS 当作已包含振动内能的自由能。具体熵方法、已找到的 MOF 专用模型及公开 MIL-53 数据核验见[框架熵方案](framework_entropy_plan.zh-CN.md)。

小于 `−imaginary_tolerance_THz` 的虚频返回 `unstable_phonons` 和 `thermal:null`，不取绝对值。接近零的频率按显式 cutoff 排除；这对低频柔性模式很敏感，必须收敛阈值与 q 网格，不能单靠“没有大虚频”证明稳定。默认阈值仅为数值初筛参数。网格中未采样到的不稳定波矢仍可能遗漏。

## 4. 准谐自由能 QHA

先对**同一势阱**的至少 5 个体积点生成、审核结构，固定体积优化内部坐标，分别运行上述声子计算；把各点 `F_total(T)` 组织为：

```json
{
  "slices": [{
    "phase": "closed", "temperature_K": 300, "pressure_bar": 1,
    "volume_A3_fu": [90, 95, 100, 105, 110],
    "F_total_kJ_mol_fu": [1, 0.25, 0, 0.25, 1],
    "all_phonons_stable": true,
    "source": "SYNTHETIC schema example; replace with computed static+phonon records"
  }]
}
```

`qha` 实现 `G_i(T,P)=min_V[F_i(T,V)+PV]`。采用分段保形 PCHIP，并在每段内和所有节点检查最小值；不是 Birch–Murnaghan EOS 拟合。极小值落在体积边界时标记 `boundary_minimum`，应扩大体积范围，不能外推。开/闭孔分别最小化，禁止跨相强行做一个单井拟合。

当前自动化的是力常数/谐振 F 与**给定 F(V,T) 表的准谐最小化**；没有自动生成、优化并审核所有体积结构。QHA 输出 G 已含 PV，若传入后续 osmotic 模型，必须先还原为所选体积的 F，避免再加一遍 PV。更严格的吸附诱导体积变化应最小化 `F(V,T)+PV+Ω_ads(V,T,f)`，不能只优化空框架后假定吸附下体积不变。

## 5. 框架 MC 热力学积分：实际采样但约束明确

`mc-ti` 使用 calculator/states，并增加：

```json
{
  "mc": {
    "temperature_K": 300, "spring_eV_A2": 2.0, "bound_A": 0.3,
    "lambdas": [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1],
    "warmup": 1000, "samples": 5000, "stride": 5,
    "step_A": 0.05, "blocks": 10, "seed": 32
  }
}
```

这是实际 Metropolis MC，不是只积分用户手填的均值。每个 λ 下调用同一 ASE 框架 ML 势，采用单 Cartesian 坐标的对称位移提议及 Boltzmann 接受率。计数单位是**单坐标提议**，不是全框架 sweep；大晶胞需要大幅增加步数。默认参数不保证收敛。

为让参考配分函数精确且不遗漏约束，本版采用：固定晶胞；固定 atom 0；其余每个坐标相对输入结构限制在 `[-a,a]`；所有原子有标签；积分测度 `∏dx/(1 Å)`。不自动识别真实开/闭孔势阱。

设 d=3(N−1)、k 为弹簧常数、β=1/(kBT)：

`U_ref(x)=k Σx_j²/2`

`Z_ref=[sqrt(2π/(βk)) erf(a sqrt(βk/2))/(1 Å)]^d`

`F_ref=−kBT ln Z_ref`

`U_λ=(1−λ)U_ref+λU_ML`

`F_conf=F_ref+∫₀¹ <U_ML−U_ref>_λ dλ`

参考能和 MC 使用同一硬壁范围，因此有限界面的 erf 修正不能省掉。结果单位 eV/输入晶胞，乘 96.4853321233/Z 得到 kJ/mol-f.u.；**这仍是受约束的构型自由能**。

不包含动量积分、解除锚定/质心约束、不可分辨性/晶体置换或有限尺寸修正，不能称为完整无约束 Frenkel–Ladd ΔF。比较不同晶胞大小、不同固定原子、不同化学组成或不同硬壁时，不能假设这些项自动抵消。即使两相组成相同，也要检查锚定及边界依赖。

输出 λ 窗口均值、接受率、分块 SE 和梯形积分。SE 不包含 λ 离散误差；连续窗口可能相关，需多种 λ 网格、不同初态/正反向扫描、独立 seeds、不同 block 长度及范围 a 收敛。框架 MC 周期不是实验时间。

目前未实现完全自动的无约束固体绝对自由能、变胞框架 MC、umbrella/metadynamics 或完整采样路径的相态审核。这些不应标成“材料已预测”。

## 6. 真实储气相与机械项

`reservoir` 的最小 spec：

```json
{"temperature_K": 268, "pressure_bar": [0.1, 1, 10, 19], "reservoir_backend": "coolprop"}
```

默认气体 CO2；VOC 必须提供 `gas_profile` 路径，沿用有来源的饱和压规则，可在 profile 内明确 `coolprop_fluid` 名称，例如 `Ethanol`。混合流体不支持。CoolProp 使用 HEOS，不任意选择立方 EOS 的最大根；液相/两相、超饱和及三相点以下范围拒绝。等于 Psat 时使用饱和**蒸汽**密度计算 fugacity，支持饱和端点。

RASPA 输入直接增加 `--reservoir coolprop`：保存温度、实际 P、φ、f=φP、软件版本；将 φ 写入 `FugacityCoefficient`。旧 `quick_free_energy_screen.py` 仍维持 ideal-only 约定，不改变归档结果。

后处理一律以真实模拟时的 `fugacity_bar` 为吸附数据横轴：

- 若旧模拟设 φ=1，数据横轴就是旧 P；在分析中视为 N(f)，按新 EOS 的 f(P_real) 查询这条曲线。
- 若新模拟已经使用 CoolProp，横轴是已记录的 φP，不能再次乘 φ。
- 禁止在原 P 处保留 N 不变、只换积分轴。目标 f 超过数据范围就拒绝，不外推。

机械项为 `P[bar] V[Å³/f.u.] × 6.02214076e−5` kJ/mol-f.u.，自动加到每相 Ω；相减即 PΔV。当前假定外压等于储气相压力、各状态体积固定、无额外应力场。体积可直接采用 host-energy 输出；晶胞体积必须先除 Z。

## 7. 多相、多分支及温度切片

`phase-diagram` spec，以下全部是格式演示数值：

```json
{
  "slices": [{
    "temperature_K": 300, "pressure_bar": [0.01, 0.1, 1, 10], "reservoir_backend": "ideal",
    "states": [
      {"name":"closed", "host_F_kJ_mol_fu":0, "host_source":"synthetic", "volume_A3_fu":100,
       "fugacity_bar":[0.001,0.01,0.1,1,10], "loading_molecules_fu":[0.0001,0.001,0.01,0.1,0.5], "henry_confirmed":true},
      {"name":"open", "host_F_kJ_mol_fu":3, "host_source":"synthetic", "volume_A3_fu":200,
       "fugacity_bar":[0.001,0.01,0.1,1,10], "loading_molecules_fu":[0.001,0.01,0.1,1,2], "henry_confirmed":true}
    ]
  }]
}
```

每个温度独立提供对应 host F 与等温线，不复用其他温度的 N。允许加入更多相。`host_F` 可统一减去同一参考常数，须为不含机械 PV 的 Helmholtz 自由能。表中 loading 已按 molecules/f.u. 归一，与旧 CSV 的 molecules/cell 不同。`henry_confirmed` 是用户依据低压收敛作出的判断，不是自动验证。

对每相计算 `Ω_i=F_host,i+PV_i−RT∫N_i dlnf`，输出节点全局最低相和每对相的全部分段插值根，标明交叉是否处于全局稳定包络上，并给出交叉前后稳定相。亚稳两相交叉不等于真实平衡相界。每相 Ω 对 ln P 用 PCHIP，同一插值同时用于根与稳定性分类；需加密压力节点确认收敛。当前输出 JSON 数值相图，不自动画图，也不在温度之间插值。

## 8. 中间构象景观、能垒和 spinodal

推荐在 `landscape` spec 中提供与上一节单个 slice 相同的 `temperature_K/pressure_bar/states`，每个 state 另外有 `q`；加上：

```json
{"landscape_source":"validated constrained F(q) and fixed-q GCMC provenance",
 "coordinate_measure":"F(q) under stated order-parameter definition; Jacobian convention documented"}
```

至少 5 个严格递增 q，至少 2 个压力。每个 q 的 host F、体积和 N(q,f) 都是独立输入，程序构建 `Ω(q,P)`。先用相同化学组成/原子映射生成并检查中间结构，逐个清理审核；本版**不自动把两个 CIF 线性插值成可信路径**。若只提供 ML 静态 E(q)，结果只能称为静态势能路径，不是有限温度自由能景观。

也可直接输入 `q`、`pressure_bar` 和形状 `[n_pressure,n_q]` 的 `omega_kJ_mol_fu`；仍需来源、坐标测度、温度和气体档案。仅解析模型测试允许 `analytic_test:true` 不经过蒸汽压审核，不适用于材料输入。

每个压力用 cubic spline 求 `∂Ω/∂q=0`，按曲率分类极小/极大；沿此 1D 路径计算相邻极小到鞍点的差值。端点不能自动当成正常内部势阱，另列 `boundary_minima_q`。

至少 4 个压力切片时，对出现极小值数目变化的区间，用二维 `(ln P,q)` spline 联立求解 `Ω_q=0, Ω_qq=0`；检查 `Ω_qqq` 与 `Ω_q,lnP` 非零，排除普通的退化临界点。输出 `spline_spinodals`，是**给定 1D 插值景观的 saddle-node**；`candidate_spinodals` 保留粗网格范围及边界风险。

局限：路径外的低能通道、多维不稳定性、晶核形成、缺陷和速率尚未预测；有限网格可能漏掉两个抵消的极小值数目变化。必须细化 P/q、移动 q 边界并检查插值敏感性。不能把数值 spinodal 直接等同实验 P_open，也不能用一个每化学式单元的均匀变形能垒直接推算晶体成核速率。

## 9. 本次真实验证及剩余项

机器可读记录在 `examples/advanced_validation/software_checks.json`，由 `validate_advanced.py --ml` 生成，**不是 MOF 预测结果**：

| 检查 | 实际结果 |
|---|---|
| CHGNet 0.3.0 / package 0.4.2 | 硅晶胞及受应变晶胞能量/力/应力、ΔE 正常；记录权重哈希 |
| MACE-MP small / package 0.3.16 | 下载官方权重并完成同样检查 |
| CHGNet→Phonopy 4.5.0 | Si 2×2×2 超胞、4³ 网格；最低频率约 −1.62×10⁻⁷ THz；300 K F_vib≈4.2734 kJ/mol-Si |
| CoolProp 8.0.0 | CO2 268 K、19 bar：φ≈0.867332，f≈16.4793 bar，Psat≈30.3336 bar |
| MC-TI，3D 谐振解析检查 | ΔF≈0.0272784 eV；解析值≈0.0268789 eV；固定 seed=32，仍有采样及积分误差 |
| MC→实际 CHGNet | Si 3 个 λ、每窗 20 样本，仅接口联通测试，绝非收敛 F |
| 双势阱 `(q²−1)²−h q` | h_spinodal≈1.5395994；解析值 8/(3√3)≈1.5396007 |
| QHA 抛物线 | 恢复解析极小值 V=100 Å³/f.u. |

测试环境 macOS arm64、Python 3.14.4；MACE/Torch 对 Python 3.14 的 JIT 给出未来兼容性警告，虽本次能运行，日常建议 3.11/3.12。测试值会随后端版本略变。

尚未完成：具体柔性 MOF 开/闭孔 ΔF 的 DFT 校准及有限温度收敛；完整无约束/变胞 MC 绝对自由能；自动生成并验证中间路径；多维自由能采样和实验时间尺度预测；混合气相平衡。没有把这些标成已集成。

## 10. 在线检索依据（2026-09-13）

| 来源 | 已核实用途 |
|---|---|
| [CHGNet official](https://github.com/CederGroupHub/chgnet), [paper DOI](https://doi.org/10.1038/s42256-023-00716-3) | 预训练权重、ASE 能量/力/应力接口；不是 MOF 专属校准 |
| [MACE official](https://github.com/ACEsuit/mace), [foundation interface](https://github.com/ACEsuit/mace/blob/main/mace/calculators/foundations_models.py) | `mace_mp`、模型别名/权重地址、可选色散接口 |
| [MOFSimBench code/data](https://github.com/AI4ChemS/mofsim-bench), [paper DOI](https://doi.org/10.1038/s41524-025-01872-3) | MOF 上多种 uMLIP 的结构/稳定性/体模量/热容评测；calculator.py 和 calculators.yaml 提供真实接口依据 |
| [Phonopy API](https://phonopy.github.io/phonopy/phonopy-module.html), [QHA docs](https://phonopy.github.io/phonopy/qha.html) | 有限位移、声子热力学、至少 5 个体积点的准谐分析 |
| [CoolProp low-level API](https://coolprop.org/coolprop/LowLevelAPI.html) | 相态、饱和态和纯流体 fugacity |
| Frenkel & Ladd, 1984, [DOI 10.1063/1.448024](https://doi.org/10.1063/1.448024) | Einstein 晶体固体自由能热力学积分原理；当前受约束实现不是完整复现 |
| Vega & Noya, 2007, [DOI 10.1063/1.2790426](https://doi.org/10.1063/1.2790426) | Einstein molecule 参考约束与修正，说明固定粒子约束不能忽略 |
| Bousquet, Coudert & Boutin, 2012, [DOI 10.1063/1.4738776](https://doi.org/10.1063/1.4738776) | 吸附诱导结构相变的自由能景观与亚稳态 |

CrossRef 宽泛的 MACE/MOF 关键词检索给出不少不相关结果，改用官方代码和 MOFSimBench 的确切 DOI 核验。初次检索未确认 `mofs_v1.model` 下载地址；随后追溯原始论文，现已找到[作者发布的 v1/v2 权重和声子数据](https://github.com/ddmms/data/tree/main/mace-mof-0)。现有 MACE 适配器接受本地权重，但尚未实际装载该专用模型。已下载并核验其中的 MIL-53(Al) 熵表，详见[框架熵估算方案](framework_entropy_plan.zh-CN.md)。

没有找到能仅靠两个未经校准的 CIF 就可靠直接输出任意柔性 MOF `ΔF(T)`、spinodal 和实验开关门压力的一体化预训练模型。现有软件提供计算环节，不省略结构、模型和采样验证。
