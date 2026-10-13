# 柔性 MOF 快速筛选 | Agent 执行约定

[English](README.md)

目标：从闭孔/开孔 CIF 快速估算条件性平衡开孔压力 P_eq。纯 Agent/CLI 技能，
不含网页服务；快速入口不执行 MD、DFT、框架优化或有限位移声子。

## 安装

从 PyPI 安装筛选库：

```bash
python -m pip install flexible-mof-gateopen
```

核心包不依赖第三方库。按路线需要加装 extra，例如
`flexible-mof-gateopen[structure,advanced]` 或 `flexible-mof-gateopen[ml]`，
完整列表见 [pyproject.toml](pyproject.toml)。

工作流**脚本**和 Codex skill 包不在 wheel 内。需要完整流程时，克隆仓库，并在
下一节中从仓库根目录执行命令：

```bash
git clone https://github.com/owenxyzhang-seseTJ/flexible-mof-gateopen.git
cd flexible-mof-gateopen
python -m pip install -r requirements.txt
```

## 输入、执行与输出

输入：两相已审核 CIF、clean-cif 记录、每个 CIF 的 Z、气体及温度、吸附力场及
来源、压力网格。两相必须具有相同的每化学式单元组成和活化状态。

流程：clean-cif → 已有电荷或 PACMOF → RDKit/有来源的气体几何 → 静态 MACE
能量＋周期弹性网络熵 → 两相 RASPA 等温线 → 真实逸度积分＋PΔV → P_eq/敏感性。

在 skill 目录中使用已配置的 Python：

```bash
python -m pip install -r requirements.txt
python scripts/verify_bundle.py
python scripts/self_test.py
python scripts/run_fast_screen.py --spec /absolute/job.json --output-dir /absolute/new_run
```

配置从 `assets/workflow.template.json` 填写；字段见
[Agent 合同](references/agent_contract.md)。JSON 内路径相对于该 JSON；内置模型
默认从 skill 根目录解析。不同平台需要指定 `raspa.binary`。不要运行未填写的模板。

## AI 必须保留的判断

- 读取内置 `vendor/clean-cif/SKILL.md` 和清理决策，先处理无序/溶剂，再赋电荷。
  不能误删框架配体、结构水、配位物种或反离子，不能伪造审核通过记录。
- PACMOF RF-DDEC 使用随包权重及 scikit-learn 1.8.0。赋电荷保持结构精度。
- RDKit 只优化几何；吸附 LJ、电荷、原子类型和刚性分子定义需要明确来源。
- VOC 的压力必须不超过该温度的饱和蒸汽压；禁止改标签伪造真实逸度数据。
- 闭孔可不吸附，但零吸附分支必须显式声明并计算残余吸附敏感性。Widom 不能
  排除高压协同填充，也不能把采样零当作严格物理零。
- 检查 Henry、独立种子、循环数、交点附近压力分辨率。程序默认参数不是收敛证明。

结果为 `run/workflow_summary.json` 及分阶段结构、等温线、能量/代理熵和相变分析。
退出 0 表示完成条件性筛选；退出 2 表示仍有阻断项。没有根只表示采样范围内无根。
源码仓库使用 `skills/flexible-mof-gateopen` 中的打包入口；开发历史不进入 skill。

## 精度与内容边界

ΔF_screen=ΔE_static+ΔU_vib,proxy−TΔS_proxy。静态能和熵来自不同模型，尚非
经校准的材料自由能。参数敏感性不是置信区间；P_eq 不是实测开门/关门压力，
两端点不能确定迟滞、能垒或动力学。详细推导见 [原理 PDF](assets/principles.pdf)。

包含所需项目源码、MACE/PACMOF 权重、PACMOF 和 RASPA/MACE 源码、clean-cif、
清单/哈希/许可证、依赖记录与原理文档。不包含试验过程和网页。原理文档同时提供
可检索 PDF 和可复现的 [XeLaTeX 源码](references/principles.tex)，便于检查中文字体
与公式。Python 依赖和操作系统动态库仍需安装；附带 RASPA 二进制仅针对 macOS arm64，
不宣称跨平台免装。
