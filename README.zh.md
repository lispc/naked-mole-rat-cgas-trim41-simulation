# cGAS-TRIM41 分子动力学研究

> 计算探究裸鼹鼠（*Heterocephalus glaber*）cGAS 中四个氨基酸变异如何影响 TRIM41 介导的泛素化。
>
> 基于：*A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging* (Chen et al., Science 2025)

## 项目结构

```
├── AGENTS.md              # AI agent 指令
├── GEMINI.md              # 指向 AGENTS.md 的符号链接
├── README.md              # 本文件（英文版）
├── README.zh.md           # 本文件（中文版）
├── main.pdf               # 项目手稿
│
├── sequences/             # 用于结构预测的 FASTA 文件
│   ├── Hsap_cGAS_WT.fasta
│   ├── Hsap_cGAS_4mut.fasta
│   ├── Hgal_cGAS_WT.fasta
│   ├── Hgal_cGAS_4mut_rev.fasta
│   └── TRIM41_WT.fasta
│
├── structures/            # 预测结构和对接结果
│   ├── af3_raw/           # 原始 AlphaFold3 输出
│   └── docking/           # 蛋白-蛋白对接结果
│
├── scripts/               # 分析和模拟脚本
│   ├── README.md          # 完整脚本索引
│   └── AGENTS.md          # 脚本特定规范
│
├── data/
│   ├── md_runs/           # OpenMM MD 轨迹
│   │   ├── Hsap_WT/
│   │   ├── Hsap_4mut/
│   │   ├── Hgal_WT/
│   │   ├── Hgal_4mut_rev/
│   │   └── Hsap_WT_S305phos/
│   ├── md_runs_gmx/       # GROMACS MD 轨迹（旧版转换）
│   ├── md_runs_gmx2026/   # GROMACS 2026 原生力场验证
│   └── analysis/          # 分析输出（JSON、PNG、NPZ）
│
├── figures/               # 可直接用于发表的图表
│
└── docs/                  # 文档（见下方索引）
    ├── README.md
    ├── 00-project/        # 项目日志和元信息
    ├── 10-reports/        # 已完成的实验报告
    ├── 20-protocols/      # 方案和执行计划
    ├── 30-diagnostics/    # 故障排查和差异诊断
    ├── 40-reviews/        # 外部同行评审反馈
    └── 50-infra/          # 硬件、软件、环境
```

## 文档索引

| 文档 | 内容 |
|------|------|
| [`docs/00-project/project_log.md`](docs/00-project/project_log.md) | **主日志**：当前状态、时间线和决策 |
| [`docs/10-reports/docking_report.md`](docs/10-reports/docking_report.md) | 蛋白-蛋白对接：ClusPro、SDOCK2.0、LightDock、Rosetta |
| [`docs/10-reports/af3_report.md`](docs/10-reports/af3_report.md) | AlphaFold3 结构预测：ipTM/pTM、突变映射 |
| [`docs/10-reports/interface_analysis_report.md`](docs/10-reports/interface_analysis_report.md) | 界面分析和变构机制证据 |
| [`docs/20-protocols/phosphorylation_md_plan.md`](docs/20-protocols/phosphorylation_md_plan.md) | 磷酸化 MD 方案和运行记录 |
| [`docs/30-diagnostics/gromacs_openmm_divergence_diagnosis.md`](docs/30-diagnostics/gromacs_openmm_divergence_diagnosis.md) | GROMACS 与 OpenMM 差异：CMAP 修复、PBC 陷阱 |
| [`docs/50-infra/hardware_benchmark.md`](docs/50-infra/hardware_benchmark.md) | 硬件基准和 MD 性能（~152 ns/day） |
| [`docs/50-infra/software_versions.md`](docs/50-infra/software_versions.md) | 软件版本锁定文件 |

完整索引见 [`docs/README.md`](docs/README.md)。

## 模拟体系

所有体系均使用 cGAS C-末端结构域构建体（残基 200-554，355 个氨基酸）与 TRIM41 SPRY 的复合物，采用 ff19SB + OPC 水模型制备。

| 体系 | cGAS | TRIM41 | 突变 | 副本数 | 状态 |
|------|------|--------|------|--------|------|
| **Hsap_WT** | 人源 WT | WT | — | 3 | 已完成（3 × 200 ns） |
| **Hsap_4mut** | 人源 | WT | D431S, K479E, L495Y, K498T | 3 | 已完成（3 × 200 ns） |
| **Hgal_WT** | 裸鼹鼠 WT | WT | — | 3 | 已完成（3 × 200 ns） |
| **Hgal_4mut_rev** | 裸鼹鼠 | WT | S463D, E511K, Y527L, T530K | 2 | 已完成（2 × 200 ns） |
| **Hsap_WT_S305phos** | 人源 WT, SEP@305 | WT | — | 3 | ✅ 完成（3 × 200 ns） |
| **Hsap_WT_S305E** | 人源 WT, S305E | WT | — | 3 | 进行中（0 ns / 200 ns） |

*NMR = 裸鼹鼠（*Heterocephalus glaber*）。论文编号使用 NMR 坐标（554 aa）；人源对应 522 aa。*

**跨引擎验证：**
- GROMACS 2026 原生 `amber19sb.ff` + OPC：Hsap_WT rep1，~83 ns / 200 ns（进行中）。CMAP 修复已验证：COM/Rg 与 OpenMM 偏差 <1%，RMSD 比率 1.34×。

**Boltz-2 验证：**
- 本地 Boltz-2（v2.2.1）安装在 `boltz` conda 环境中。全长 cGAS-TRIM41：ipTM 0.17（与 AF3 一致）。截断构建体：ipTM 0.33，但 TRIM41 放置位置与 AF3 定性不同（RMSD 21.34 Å，Jaccard 0.00）。证实了对此瞬态复合物采用对接+MD 策略的正确性。

## 关键发现

### 1. 突变位点不在物理界面上

四个突变位点（NMR 编号中的残基 463、511、527、530）在序列上距离界面区域（残基 228-266）**超过 200 个残基**，在空间上距离 **30-39 Å**。它们无法直接与 TRIM41 接触。

| 突变 | NMR 残基编号 | 距最近界面残基的距离 | 距 TRIM41 的距离 |
|------|-------------|---------------------|-----------------|
| S463 | 482 | 32 个残基 | ~29 Å |
| E511 | 530 | 80 个残基 | ~31 Å |
| Y527 | 546 | 96 个残基 | ~31 Å |
| T530 | 549 | 99 个残基 | ~39 Å |

### 2. 界面位于 cGAS N-末端

物理结合界面涉及 cGAS 残基 **228-266** 和 TRIM41 残基 **82-190**：
- TRIM41-157/158 ↔ cGAS-258/259（占有率 >60%，最稳定）
- TRIM41-187 ↔ cGAS-236（占有率 ~39%）

### 3. 人源体系比裸鼹鼠更具动态性

| 指标 | Hgal_WT | Hsap_WT | Hsap_4mut |
|------|---------|---------|-----------|
| RMSD | 5.16 ± 0.72 Å | **8.94 ± 1.58 Å** | **9.76 ± 2.21 Å** |
| COM | 37.4 ± 0.9 Å | **46.6 ± 2.4 Å** | **49.0 ± 2.8 Å** |

4mut 使人源复合物**更不稳定**（RMSD +9%，COM +5%，p < 1e-30）。

### 4. 活性位点几何结构在物种间存在差异

| 位点 | Hgal_WT | Hgal_4mut_rev | Hsap_WT | Hsap_4mut |
|------|---------|---------------|---------|-----------|
| S463/D431 | 29.4 Å | **18.6 Å** (−10.8) | 29.4 Å | 32.0 Å (+2.6) |
| E511/K479 | 31.5 Å | **21.9 Å** (−9.7) | 36.6 Å | 40.0 Å (+3.4) |

**Hgal_4mut_rev** 使活性位点**更靠近** TRIM41；**Hsap_4mut** 使其**远离**。

### 5. Lys-334 是人源 cGAS 中首要的泛素化候选位点

赖氨酸可及性分析显示 **Lys-334** 是距离 TRIM41 RING 结构域最近的 cGAS 赖氨酸：
- Hsap WT: 10.4 Å
- Hsap 4mut: **6.4 Å** (−4.0 Å)

4mut 可能将泛素化**集中**到这一单个位点。

### 6. S305 磷酸化导致复合物解离（出乎意料）

S305-phos（~138 ns）的初步分析显示所有三个副本均发生**解离**：

| 体系 | COM (Å) | Rg (Å) | 行为 |
|------|---------|--------|------|
| WT（参考） | 45.1 ± 3.2 | 31.1 ± 1.3 | 稳定结合 |
| S305-phos rep1 | **67.6 ± 1.4** | 38.9 ± 0.6 | 正在解离 |
| S305-phos rep2 | **89.8 ± 10.0** | 48.4 ± 4.3 | 完全解离 |
| S305-phos rep3 | **71.3 ± 3.4** | 40.2 ± 1.4 | 正在解离 |

> 注：截至 2026-04-23，3 个 replica 均在约 140 ns / 200 ns。

这与 Zhen et al. (2023) 的结论相矛盾，后者报告 CHK2 在 S305 处的磷酸化**促进** cGAS-TRIM41 结合。可能的解释：(1) 力场对溶液中 −2 电荷排斥的高估，(2) 需要细胞核/DNA 环境才能稳定，(3) 解离是 DNA 介导招募之前的中间态。

### 变构机制假说

由于突变位点远离界面，其作用必须是**变构的**：

```
C-末端 4mut  →  改变全局动态  →  改变 RING-赖氨酸几何结构  →  差异泛素化
```

## 快速开始

### 环境配置

```bash
# 激活 MD 环境
conda activate cgas-md
python --version  # Python 3.11

# 对于 GROMACS
conda activate gmx
gmx --version  # GROMACS 2026.0
```

### 构建 MD 体系

```bash
conda activate cgas-md
python scripts/build_system.py \
    --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
    --name Hgal_domain
```

### 运行生产 MD

```bash
conda activate cgas-md
CUDA_VISIBLE_DEVICES=0 python scripts/run_md.py \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
    --name Hsap_WT_rep1 \
    --outdir data/md_runs/Hsap_WT/rep1 \
    --prod-ns 200 \
    --platform CUDA
```

**多 GPU 副本：** 使用 `CUDA_VISIBLE_DEVICES` 隔离 GPU。脚本硬编码了 `CudaDeviceIndex='0'`。

### 分析单个体系

```bash
conda activate cgas-md
python scripts/analyze_system.py \
    --system Hsap_WT \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --trajectories data/md_runs/Hsap_WT/rep1/*.dcd \
    --replica-names rep1 \
    --cgas-range 219 573 \
    --active-sites '{"D431": 450, "K479": 498, "L495": 514, "K498": 517}' \
    --dt-ns 0.1 \
    --outdir data/analysis/my_run
```

### 比较两个体系

```bash
conda activate cgas-md
python scripts/compare_systems.py \
    --a data/analysis/run_A/Hsap_WT_summary.json \
    --b data/analysis/run_B/Hsap_4mut_summary.json \
    --name-a Hsap_WT --name-b Hsap_4mut \
    --outdir data/analysis/comparison
```

## 环境

### Conda 环境

| 环境 | 用途 | 关键包 |
|------|------|--------|
| `cgas-md` | OpenMM MD、AmberTools、分析 | OpenMM 8.5.1, AmberTools 24.8, MDAnalysis 2.10.0, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.8 |
| `gmx` | GROMACS 2026.0 CUDA | GROMACS 2026.0 (nompi_cuda) |
| `rosetta` | PyRosetta 对接 | PyRosetta 2025.06 |
| `boltz` | **Boltz-2 结构预测** | **Boltz-2.2.1, PyTorch 2.11.0, CUDA 13** |

**始终使用 conda 环境。不要使用系统 Python。**

### 硬件

| 组件 | 规格 |
|------|------|
| CPU | AMD EPYC 7702, 64 核 |
| GPU | 4× NVIDIA RTX 3090 (各 24 GB 显存) |
| 操作系统 | Linux x86_64, CUDA 13.0 |
| 性能 | ~152 ns/day @ 116k 原子（OpenMM，截断结构域） |

## 引用

Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. **Science**.
