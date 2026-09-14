# 空框架开孔/闭孔熵差：文献依据、实际数据核验与可执行方案

日期：2026-09-13。范围仅限空框架，不涉及客体吸附熵。目标是同组成、同活化状态下的 `ΔS_host(T)=S_open(T)−S_closed(T)`，以及由此参与决定的 `ΔF_host(T)`。

**最新用户约束：快速筛选，不使用 MD/DFT。**当前已切换到
[静态 ML 能量＋弹性网络代理熵](fast_entropy_screening.zh-CN.md)，并完成数秒级实跑。
下文保留方法背景，不表示现在会启动 MD。MACE-MOF0-v2 已在后续工作中载入并用于
两相静态能量/声子测试，后续状态见新文档；此前“尚未装载”的表述属于历史记录。

## 1. 结论：可以估算，不需要先训练一个直接预测熵的模型

首选 **预训练保守 ML 势 → 两相的声子 → 量子振动熵**。若存在显著热膨胀，分别对两相做准谐处理；若出现真实软模/强非谐性，则以短程 ML-MD 的协方差熵做筛查，最后用相内热力学积分核验。

要纠正一个容易产生的简化：不能只拿静态 `ΔE` 减去 `TΔS` 就称为完整有限温度结果。谐振层面应计算：

`ΔF_host(T) ≈ ΔE_static + ΔF_vib(T)`

其中 `ΔF_vib = ΔU_vib − TΔS_vib`，包括零点能和振动热能。声子软件能同时给出这些量，无需遗漏其中任意一项。

## 2. 找到了现成模型，而且有 MIL-53(Al) 的真实公开热力学数据

Elena 等人的 [2025 年 MACE-MP-MOF0 论文](https://doi.org/10.1038/s41524-025-01611-8) 专门针对 MOF 声子及准谐性质；模型基于 127 种 MOF 的数据微调，涉及 24 种元素。作者明确报告基础 MACE-MP-0b 的虚频/声子问题，以及微调后的改善；也指出 QHA 仍不能充分描述 MOF 的非谐效应。

这轮追溯原论文 Data/Code availability 后确认：

- [官方数据与模型仓库](https://github.com/ddmms/data/tree/main/mace-mof-0)
- [v1 模型](https://github.com/ddmms/data/blob/main/mace-mof-0/v1/mofs_v1.model)
- [v2 模型](https://github.com/ddmms/data/blob/main/mace-mof-0/v2/mofs_v2.model)
- [MIL-53 声子数据包](https://github.com/ddmms/data/blob/main/mace-mof-0/phonons/MIL-53.zip)

作者仓库采用 CC BY 4.0，使用模型/数据需要引用原论文。上轮“未确认 MOF 专用权重地址”的状态现已更新为**已确认公开来源**。本轮没有新下载/装载这些权重做推断，而是下载并分析了原作者声子数据。

### 实际核验结果，不是材料预测承诺

公开数据包里有 MACE-MP-0b、MACE-MP-MOF0 v1/v2 对一个 MIL-53 条目的声子带和热力学表。声子 primitive cell 组成核实为 `Al4C32H20O20`，即 4 个 `Al(OH)(BDC)` 化学式单元，共 76 个原子。

| 公开数据中的模型 | 300 K 振动熵，J mol-f.u.⁻¹ K⁻¹ | 采样声子路径上的最低频率，THz | 解释 |
|---|---:|---:|---|
| MACE-MP-0b | 215.6712 | −0.289918 | 有显著虚频，不作为稳定谐振参考 |
| MACE-MP-MOF0 v1 | 212.6807 | −5.62×10⁻⁵ | 在该采样路径上只剩数值量级的负频 |
| MACE-MP-MOF0 v2 | 212.4884 | −2.67×10⁻⁵ | 在该采样路径上只剩数值量级的负频 |

这些是**读取原作者模型计算的单个结构熵**，不是本轮重新计算的熵，也不是开孔/闭孔的熵差。两个模型相近不能证明误差小；不同模型的松弛几何也可能不同。带路径没有大虚频不保证整个 Brillouin 区都稳定。

核验还发现热力学表头将 F 标为 J/mol，而数值满足 **F 为 kJ/mol、S 为 J/mol/K** 的约定。以 v2 的 290/310 K 中心差分为例：`−1000×dF/dT≈849.9204 J/mol-primitive/K`，与表中 300 K 的 `S=849.9538` 一致；若直接采用 F 表头的 J/mol，会错 1000 倍。再除以 Z=4，得到 212.4884 J/mol-f.u./K。保留了原表头、哈希和单位诊断，不静默改写源数据。

可复核文件：`examples/advanced_validation/published_mil53_entropy_audit.json`。

```bash
python scripts/audit_published_mil53_entropy.py \
  --archive output/entropy_research/MIL-53-api.zip \
  --out output/rechecked-entropy.json
```

脚本绑定此次源压缩包 SHA256；不修改、解压覆盖或重新计算源结构。

## 3. 推荐路线 A：两相量子声子熵，最适合作为第一版

对每相，用同一 ML 势和色散约定、对应稳定体积，求力常数及频率 ν。不能给已在含 D3 数据上训练的模型不加判断地再加 D3。

令 `x=hν/(kBT)`，单个量子谐振模式的熵为：

`s(ν,T)=kB[x/(exp(x)−1)−ln(1−exp(−x))]`

各 q 点按 Brillouin 区权重加总所有声子支，除以对应 primitive cell 的化学式单元数，即得到每 f.u. 的 S。常规 Phonopy `thermal_properties` 已完成这一统计求和。最后相减，不必手工逐一配对两相的振动模。

具体执行顺序：

1. 用 clean-cif 审查开孔/闭孔的无序、溶剂、H 和组成；不能把湿态与干态相减。
2. 用同一势在各自势阱内优化。先固定晶胞优化原子位置；若做变胞优化，需监测相态，避免两相都塌入同一个势阱。
3. 有限位移初始可比较 0.005 Å 和 0.01 Å；检查超胞和 q 网格，不能只依靠 Γ 点。原论文使用超 Γ 点网格，MIL-53 归档配置记录了 3×1×2 超胞和 0.001 Å 位移；这些是原作者设置，不是通用收敛保证。
4. 输出 `S_open(T), S_closed(T), ΔS(T), ΔE_static, ΔF_vib(T), ΔF_total(T)`，以及每相最低频率、归一化、模型权重与误差来源。
5. 在目标温度附近重复，并以 ΔS/ΔF 的变化而非单独 Cv 的变化判定数值收敛。

已有 `advanced_free_energy.py phonons` 可调用 CHGNet 或本地 MACE 模型并输出 `S_J_mol_fu_K`。提供官方 v2 权重路径时可使用现有 MACE 本地模型接口；**这不等于本次已验证 v2 运行接口或目标两相的自由能**。

### 同时计算力场 ΔE，并直接合成为 ΔF

根据后续要求，已让同一次 `phonons` 命令自动合并两相静态能和振动项。将闭孔放在 `states` 第一项、开孔放在第二项；程序使用同一个 calculator 实例计算两相，不混用不同势的能量零点。新增 `harmonic_state_differences` 输出：

- `delta_E_static_kJ_mol_fu`：力场静态开孔减闭孔能差；
- `delta_U_vib_kJ_mol_fu`：含零点能/热激发的振动内能差；
- `delta_S_vib_J_mol_fu_K` 及 `minus_T_delta_S_kJ_mol_fu`；
- `delta_F_vib_kJ_mol_fu`；
- `delta_F_host_harmonic_kJ_mol_fu = delta_E_static + delta_U_vib − T delta_S_vib`。

显著虚频会阻止生成配对 ΔF；配对门槛为 −10⁻⁴ THz，仍需专门检验低频数值收敛。此输出是同势、同组成的**谐振自由能估计**，不含非谐、构型熵或机械 PV。

这里“框架力场”可使用预训练 ML 势。若使用传统力场，则需经过适配的 MOF-FF、QuickFF 或完整 UFF4MOF 参数，包括键、键角、扭转/面外项及非键项；现有 RASPA 用于刚性框架吸附的 LJ/电荷参数不能自动提供可信的框架变形能。尚未接入通用的传统键合 MOF 力场生成器。

### 低频模式为什么关键？

300 K 下，一个模式从 20 cm⁻¹ 降到 10 cm⁻¹，按上式算出的熵增加约 `5.7608 J/mol-mode/K`，其 `−TΔS` 贡献约 `−1.7282 kJ/mol-mode`。这是计算的单模式例子，不是 MIL-53 数据；实际需 q 权重和每 f.u. 的模式归一。

这说明少量软模就可能显著改变两相自由能。不能因“高频大部分谱相似”就忽略熵，也不能把开孔相熵较高当成普遍定律。

## 4. 路线 B：准谐方法处理热膨胀，但不是所有非谐性

分别对两个势阱取至少 5 个体积点，固定体积优化内部坐标，计算：

`F_i(T,V)=E_i(V)+F_vib,i(T,V)`

在共同外压 P 下，各自最小化 `G_i(T,P)=min_V[F_i(T,V)+PV]`。平衡体积处，可由 `S_i(T,P)=−(∂G_i/∂T)_P` 求熵；不要对任意变化体积路径的 F 直接求导，再当成固定体积熵。

不要跨开闭孔两个势阱拟合成单个平滑 EOS。QHA 能描述频率随体积的变化，但不能自动包含配体大幅跳转、多井混合或强烈的温度依赖非谐振动。原 MACE-MP-MOF0 论文自身也指出这一限制。

## 5. 路线 C：ML-MD 协方差熵，作为非谐问题的快速筛查

Huang 与 Widom 的 [2022 年研究](https://doi.org/10.3390/e24050618) 从原子位移分布和协方差估算晶体振动熵，并讨论了从协方差提取有效振动态密度及量子统计处理。

可对每个相做相内 NVT ML-MD，计算去除整体平移后的位移协方差：

`C=<δr δrᵀ>`

若两相具有相同的原子集合、独立自由度、测度和约束，且用经典高斯近似，则构型熵差中相同常数抵消：

`ΔS_config,Gaussian = (R/2) ln(det C_open / det C_closed)`

这给出一个比完整 TI 更轻量的估计；但必须满足：

- 周期轨迹正确解缠，原子映射一致；使用相同组成的可比超胞，再归一为每 f.u.。
- 去除相同的刚体平移自由度；固定晶胞的振动分析不要混入未定义的盒子形变。
- 有效独立样本足够估计协方差谱；秩亏不能靠任意加小常数就宣称收敛。
- 多时间长度、独立轨迹和分块检查。可以从数十 ps 试跑，但不能预先保证该时长足以采到慢转动和呼吸。
- 高斯近似对给定协方差的单相位置熵给出上界；两个上界之差**不是**熵差的上界或下界。

经典 MD 在室温不能正确量子化所有 C–H/O–H 高频模式；有效频率后加 Bose 统计是近似，并不自动成为精确量子非谐熵。此路线目前是建议，尚未接入自动轨迹处理器。

## 6. 路线 D：相内热力学积分，处理真正的非谐项

对每相建立同一约束下可解析的谐振参考，再对参考势到 ML 势作耦合：

`U_λ=(1−λ)U_ref+λU_ML`

`F_ML−F_ref=∫₀¹ <U_ML−U_ref>_λ dλ`

参见 Frenkel–Ladd 固体参考体系方法 [1984](https://doi.org/10.1063/1.448024)、Einstein molecule 约束修正 [2007](https://doi.org/10.1063/1.2790426)，以及 Cheng–Ceriotti 的固体绝对自由能路线 [2018](https://doi.org/10.1103/PhysRevB.97.054102)。

有两种获得熵的方式：

- 在相同固定体积/统计约束下，以完整匹配的 `S=(<U>−F)/T` 计算；如果 F 是仅构型自由能，就只能配仅势能得到构型熵，不能混入动量项后仍称相同量。
- 在相同约束下计算多个邻近温度的 F，取 `S=−∂F/∂T`；对压力固定的 G 则用 `−∂G/∂T|P`。温差太小会放大 MC 噪声，须检查步长和拟合误差。

常用折中是 `F≈F_harmonic,quantum + (F_ML,classical−F_harmonic,classical)`，即量子谐振基础加经典非谐修正。它要求相同参考、约束和相域；不能称为精确量子非谐结果。高频非谐/核量子效应重要时，需更昂贵的路径积分或其他量子非谐处理。

现有项目的 MC-TI 是固定一个原子、固定晶胞和有限位移范围的构型积分。解除约束和有限尺寸修正没有自动补齐，所以尚不能把其两相结果直接称为完整框架 ΔF。这一限制不会因采样步数增加而自动消失。

## 7. 构象熵与无序：额外检查，不能重复计算

若每相内存在多个可及构型，完整配分函数需对构型求和。离散分布项为 `S_mix=−R Σp_i ln p_i`，但必须结合各构型内部熵，并明确 p_i 的采样与归一。

不能直接把晶体学 occupancy 当作平衡概率；也不能先对完整多井轨迹估计熵，再把同一构型混合熵额外加一遍。清理 CIF 得到一个有序代表结构后，应保留其他合理无序构型用于敏感性分析，不能把“结构已清理”误认为构象熵为零。

## 8. 接下来最小且有效的实际计算

**优先以公开 MACE-MP-MOF0-v2 为候选，计算经过审核的 MIL-53(Al) 开孔/闭孔两相的声子熵差。**这比直接上长时间框架 MC 更适合得到第一版结果。

第一版结果交付应是某个明确温度下的 `ΔS_vib` 与 `ΔF_host^harmonic`，附带软模/数值收敛和模型适用性说明。若一相存在无法通过合理相内优化消除的真实虚频，就转向非谐采样，不人为取绝对值或截成小正数。

本轮已经证明公开资源中有可用的 MIL-53(Al) 熵数据及 MOF 专用模型；尚未得到一对已确认相态、同条件且收敛的开闭孔结果。因此不在此填入未经计算的 ΔS 或 ΔF。

## 文献与查证范围

1. Elena et al. *Machine learned potential for high-throughput phonon calculations of metal—organic frameworks*. npj Comput. Mater. 11, 125 (2025). [10.1038/s41524-025-01611-8](https://doi.org/10.1038/s41524-025-01611-8)。已核对正文相关段落、官方仓库及实际数据。
2. Hoffman et al. *Elucidating the Vibrational Fingerprint of the Flexible Metal–Organic Framework MIL-53(Al) Using a Combined Experimental/Computational Approach*. J. Phys. Chem. C 122, 2734–2746 (2018). [10.1021/acs.jpcc.7b11031](https://doi.org/10.1021/acs.jpcc.7b11031)。已核对元数据与摘要；证明两相低频谱及呼吸相关集体模不同，不作为本轮熵数值来源。
3. Cockayne. *Thermodynamics of the Flexible Metal–Organic Framework Material MIL-53(Cr) From First-Principles*. J. Phys. Chem. C 121, 4312–4317 (2017). [10.1021/acs.jpcc.6b11692](https://doi.org/10.1021/acs.jpcc.6b11692)。本轮核实元数据，未取得可据以引用的全文熵值；Cr 结果不能直接移植到 Al。
4. Huang & Widom. *Vibrational Entropy of Crystalline Solids from Covariance of Atomic Displacements*. Entropy 24, 618 (2022). [10.3390/e24050618](https://doi.org/10.3390/e24050618)。已核对摘要；非针对 MOF 的既定精度保证。
5. [Phonopy 统计热力学公式](https://phonopy.github.io/phonopy/formulation.html)。已核对 F 与 S 的表达式。
6. Cheng & Ceriotti. *Computing the absolute Gibbs free energy in atomistic simulations: Applications to defects in solids*. Phys. Rev. B 97, 054102 (2018). [10.1103/PhysRevB.97.054102](https://doi.org/10.1103/PhysRevB.97.054102)。本轮核实元数据，方法作为固体自由能参考，未声称全文复现。

宽泛 CrossRef 查询产生了较多无关结果，已排除；arXiv 接口出现 HTTP 429，改用确切 DOI、期刊及作者 GitHub。源压缩包 raw 地址首次超时，改用 GitHub blob raw API 成功下载并核对 SHA256。
