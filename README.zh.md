# cGAS-TRIM41 分子动力学研究

> 通过计算模拟探究裸鼹鼠 cGAS 四个氨基酸变异（D431S/K479E/L495Y/K498T）如何通过长程变构效应调控 TRIM41 介导的泛素化。
>
> 基于：Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. **Science**.

---

## 最终产出

| 文档 | 说明 |
|------|------|
| 📄 **[`paper/paper_final.md`](paper/paper_final.md)** | **论文终稿**——核心结论与支撑逻辑 |
| 📋 **[`docs/00-project/lab_report.md`](docs/00-project/lab_report.md)** | **实验报告**——完整实验记录、决策因果与经验教训 |
| 📝 **[`docs/00-project/retrospective.md`](docs/00-project/retrospective.md)** | **项目回顾**——得失总结与关键转折点 |

---

## 核心发现

1. **四个突变位点不在 TRIM41 物理界面上。** 它们距结合界面 30–39 Å。真正的界面在 cGAS N 端区域。

2. **4mut 通过长程变构效应引起 N 端构象位移。** N 端界面区域位移高达 12.3 Å，而突变位点自身移动 <0.7 Å（AF3 单体对比）。

3. **结合亲和力不变。** 人源 cGAS 的 MD 几何指标（COM、contacts、Rg）无显著变化。MM-GBSA ΔΔG = 4.5 ± 10.1 kcal/mol（p = 0.50）。

4. **4mut 重塑了构象动态网络。** DCCM 显示长程耦合改变，PCA 显示构象空间分布偏移。ΔRMSF 无全局显著变化。

5. **K315 是唯一几何上可达的泛素化靶点。** 对 cGAS 全部 36 个赖氨酸的扫描显示，K315 是唯一在 Ub-G76 催化中心 25 Å 以内的残基。

6. **Umbrella sampling PMF：4mut 使 K315 自由能最小值偏移 −2.5 Å**（WT: 21.5 Å → 4mut: 19.0 Å），到达 15 Å 的自由能代价降低 1.1 kcal/mol。催化准备度（K315<19 Å + 正确攻击角）：**WT 12.6% → 4mut 53.7%（4.3 倍）。**

---

## 项目结构

```
├── paper/
│   └── paper_final.md       # 论文终稿
├── sequences/               # 所有体系的 FASTA 文件
├── structures/              # AF3 预测、对接结果、PDB 模板
├── scripts/                 # 构建、MD、分析、对接、验证脚本
│   ├── lib/                 # 共享库（路径、统计、绑图、MDA 工具）
│   └── archive/             # 已归档的废弃和实验脚本
├── data/
│   ├── md_runs/             # MD 轨迹（二元体系 + 四元复合物）
│   ├── analysis/            # 分析输出（NPZ、JSON、PNG）
│   └── structures/          # 构建的体系拓扑文件
├── figures/                 # 发表级图表
└── docs/                    # 完整文档（见 docs/README.md）
```

---

## 模拟体系

| 体系 | cGAS | TRIM41 | 副本 | 状态 |
|------|------|--------|------|------|
| Hsap_WT | 人源 WT | SPRY | 3 × 200 ns | ✅ |
| Hsap_4mut | 人源 4mut | SPRY | 3 × 200 ns | ✅ |
| Hgal_WT | 裸鼹鼠 WT | SPRY | 3 × 200 ns | ✅ |
| Hgal_4mut_rev | 裸鼹鼠 rev | SPRY | 3 × 200 ns | ✅ |
| 四元 WT | 人源 WT | RING+SPRY+E2~Ub | 1 × 50 ns | ✅ |
| 四元 4mut | 人源 4mut | RING+SPRY+E2~Ub | 1 × 50 ns | ✅ |
| US K315 (WT) | — | RING+SPRY+E2~Ub | 11 窗口 × 10 ns | ✅ |
| US K315 (4mut) | — | RING+SPRY+E2~Ub | 10 窗口 × 10 ns | ✅ |

---

## 快速开始

```bash
conda activate cgas-md

# 运行生产 MD
CUDA_VISIBLE_DEVICES=0 python scripts/02_md/run_production.py \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
    --name Hsap_WT_rep1 --outdir data/md_runs/Hsap_WT/rep1 \
    --prod-ns 200 --platform CUDA

# 四系统对比
python scripts/03_analysis/compare_four_systems_fast.py

# WHAM PMF 分析
python scripts/03_analysis/wham_pmf.py --system WT --windows 12,...,22
```

---

## 环境

| 组件 | 规格 |
|------|------|
| CPU | AMD EPYC 7702, 64 核 |
| GPU | 4× NVIDIA RTX 3090 (24 GB) |
| MD 引擎 | OpenMM 8.5.1 + GROMACS 2026.0 |
| 力场 | Amber ff19SB + OPC 水模型 |
| 分析 | MDAnalysis 2.10.0, AmberTools 24.8 |

## 引用

Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. **Science**.
