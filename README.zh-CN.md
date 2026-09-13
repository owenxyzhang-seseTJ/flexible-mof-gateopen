# 柔性 MOF 开门吸附

中文 | [English](README.md)

用于准备柔性 MOF 吸附计算、筛选平衡相变压力的 Codex skill。

**清理 CIF → 加电荷 → 客体模型 → RASPA → 吸附自由能 + 框架自由能 → P_eq**

## 安装

```bash
git clone https://github.com/owenxyzhang-seseTJ/flexible-mof-gateopen.git
cd flexible-mof-gateopen
python -m pip install -e '.[dev,structure,guest]'
cp -R skills/flexible-mof-gateopen ~/.codex/skills/
```

还需安装 **clean-cif** 技能。RASPA3、PACMOF 单独安装。仅分析自由能时只需 Python 3.10+。

## 在 Codex 中使用

> 用 $flexible-mof-gateopen 处理我的开孔/闭孔 CIF，计算 298.15 K 乙醇吸附。先用 $clean-cif 处理无序和溶剂，再用 PACMOF 加电荷、RDKit 优化客体。压力不超过有来源的 Psat(T)，按我提供的框架自由能范围给出 P_eq。

提供两相结构、气体、温度、可靠的客体模型/力场，以及框架自由能或明确标为假设的待筛选范围。

| 步骤 | 要求与输出 |
|---|---|
| 框架 | 强制审查无序和溶剂；输出 P1 CIF、双解析器审计及绑定哈希的记录 |
| 电荷 | 使用有依据的原始电荷或 PACMOF；检查电中性及几何不变 |
| 客体 | RDKit 多构象 → MMFF94s/UFF → SDF、XYZ、收敛记录 |
| VOC 压力 | 使用实际温度下有来源的饱和蒸汽压；拒绝 P > Psat |
| GCMC | 显式设置超胞/截断、Ewald、CIF 电荷；自定义气体使用对应模型 |
| 自由能 | 按统一化学式单元积分绝对吸附量，给出所需框架自由能与条件 P_eq |

清理记录缺失、文件哈希不一致或无序/溶剂未解决时，停止生成计算输入。RDKit 优化几何后仍需独立配置 RASPA 力场参数。

## 运行已有数据示例

```bash
python scripts/quick_free_energy_screen.py \
  --csv examples/xpcu5zn_co2_268K/data/fixed_host_isotherms_268K.csv \
  --temperature 268 --open-formula-units 2 --closed-formula-units 4 \
  --target-pressure 19 --host-penalty 67.9
```

X-pcu-5-Zn/CO₂ 归档数据在 19 bar 对应的框架自由能惩罚约为 **67.93 kJ/mol（框架化学式单元）**；假设该惩罚为 67.9 kJ/mol，得到 **P_eq ≈ 18.96 bar**。这是旧数据分析示例，尚未按新增清理门槛认证。[数值修正说明](references/method_boundaries.md)。

预测实验开门/关门压力还需要经验证的框架自由能景观及能垒或失稳判据。仅有两个刚性结构不能确定滞后和切换速率。

## 详细说明

- [自由能预测方法：原理、ML 路线、计算与误差](references/free_energy_methods.zh-CN.md)
- [结构、无序与溶剂处理](references/cleaning_gate.md)
- [VOC 饱和蒸汽压与 RDKit](references/vapor_and_guest.md)
- [RASPA 命令](references/raspa_input.md)
- [Skill 执行说明](skills/flexible-mof-gateopen/SKILL.md)

```bash
python -m pytest -q
python scripts/validate_example.py
python skills/flexible-mof-gateopen/scripts/quick_validate.py skills/flexible-mof-gateopen
```
