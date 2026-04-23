# cGAS-TRIM41 Molecular Dynamics Study

> Computational investigation of how 4 amino acid variants in naked mole-rat cGAS affect TRIM41-mediated ubiquitination.
> 
> Based on: *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging* (Chen et al., Science 2025)

## Project Structure

```
├── sequences/          # FASTA files for structure prediction
│   ├── Hsap_cGAS_WT.fasta
│   ├── Hsap_cGAS_4mut.fasta
│   ├── Hgal_cGAS_WT.fasta
│   ├── Hgal_cGAS_4mut_rev.fasta
│   └── TRIM41_WT.fasta
├── structures/         # Predicted structures and docking results
│   ├── af3_raw/        # Raw AF3 outputs (4 jobs)
│   └── docking/        # Protein-protein docking results
│       ├── cluspro/    # ClusPro web server results
│       ├── sdock2/     # SDOCK2.0 attempt (abandoned)
│       └── lightdock/  # LightDock results (successful)
├── scripts/            # Analysis and simulation scripts
│   ├── process_af3_results.py      # Automated AF3 result processing
│   ├── analyze_lightdock.py        # LightDock pose analysis
│   ├── build_system.py             # MD system preparation
│   ├── run_md.py                   # Production MD with OpenMM
│   ├── analyze_trajectory.py       # Trajectory analysis
│   └── run_mmpbsa.py               # MM-GBSA binding energy
├── data/
│   ├── md_runs/        # MD system files and trajectories
│   ├── analysis/       # Analysis outputs
│   └── rosetta/        # Rosetta mutation scanning
├── tools/              # Third-party tools
│   └── SDOCK2.0/       # Git clone (compiled but unused)
└── docs/               # Documentation
```

## Documentation Index

| Document | Content |
|----------|---------|
| [`docs/project_log.md`](docs/project_log.md) | **主索引**：项目概况、当前状态总览、时间线、文件速查 |
| [`docs/docking_report.md`](docs/docking_report.md) | **蛋白对接完整记录**：ClusPro（失败）、SDOCK2.0（放弃）、LightDock（成功）；3 套体系的运行参数和结果分析 |
| [`docs/af3_report.md`](docs/af3_report.md) | **AF3 结构预测**：4 个 job 的 ipTM/pTM 数据、突变映射表、AF3 序列问题发现（Job 2/4 实际提交 WT 序列） |
| [`docs/hardware_benchmark.md`](docs/hardware_benchmark.md) | **硬件与性能**：M3 Pro OpenCL 基准测试、MD 速度估算（~95 ns/day @ 80k atoms）、内存/存储估算 |
| [`docs/software_versions.md`](docs/software_versions.md) | **软件版本清单**：Conda env、pip、Homebrew 全部包的版本号 |
| [`docs/execution_plan_v1.md`](docs/execution_plan_v1.md) | **执行方案 v1.0**：6 个阶段的详细计划、时间线、人力分配 |
| [`docs/computational_workflow.md`](docs/computational_workflow.md) | **原始方案设计**：方法学选择理由、工具对比、技术路线 |
| [`docs/paper_notes_cgas_trim41.md`](docs/paper_notes_cgas_trim41.md) | **论文理解笔记**：Chen et al. Science 2025 的关键发现、实验因果链、突变信息 |
| [`docs/cluspro_submission_guide.md`](docs/cluspro_submission_guide.md) | **ClusPro 提交指南**：6 个 job 的详细参数、活性残基列表、截断后编号说明 |

## Systems Simulated

| System | cGAS | TRIM41 | Mutations |
|--------|------|--------|-----------|
| Hsap_WT | Human WT | WT | — |
| Hsap_4mut | Human → NMR | WT | C463S, K479E, L495Y, K498T |
| Hgal_WT | NMR WT | WT | — |
| Hgal_4mut_rev | NMR → Human | WT | S463C, E511K, Y527L, T530K |

*Note: NMR = naked mole-rat (Heterocephalus glaber). Paper numbering uses NMR coordinates (554 aa); human is 522 aa.*

## Key Finding (So Far)

**NMR cGAS mutations create a compact TRIM41-binding patch (~18 Å) vs human's dispersed sites (~28 Å)**

- NMR: residues 463/511/527/530 cluster into a single ~18Å patch → TRIM41 can contact all 4 simultaneously
- Human: residues 463/479 and 495/498 are ~28Å apart on opposite faces → TRIM41 cannot contact all 4 at once
- This geometric change may explain the reduced ubiquitination efficiency (not affinity, but E2-Ub transfer geometry)

## Quick Start

### Environment Setup

```bash
# Conda base
export CONDA_PATH="/Users/zhangzhuo/miniforge3"
source $CONDA_PATH/etc/profile.d/conda.sh

# MD environment (OpenMM, AmberTools, analysis)
conda activate cgas-md
python --version  # Python 3.13.13

# Docking environment (LightDock)
conda activate py311
python --version  # Python 3.11.15
```

### Process AF3 results
```bash
conda activate cgas-md
python scripts/process_af3_results.py --job-dir structures/af3_raw/job1_Hsap_WT
```

### Analyze LightDock poses
```bash
conda activate py311
python structures/docking/lightdock/analyze_lightdock.py \
    structures/docking/lightdock/Hgal_domain \
    structures/af3_raw/job3_Hgal_WT/trim41_SPRY_413-630.pdb \
    hgal
```

### Build MD system
```bash
conda activate cgas-md
python scripts/build_system.py \
    --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
    --name Hgal_domain
```

## Environment Details

### Conda Installation

| Item | Path |
|------|------|
| Conda | `/Users/zhangzhuo/miniforge3/bin/conda` |
| Conda version | 26.1.1 |
| Platform | osx-arm64 (Apple Silicon) |

### Conda Environments

#### `cgas-md` — MD Simulation & Analysis

| Package | Version |
|---------|---------|
| Python | 3.13.13 |
| OpenMM | 8.5.1 (conda-forge, apple build) |
| AmberTools | 24.8 |
| MDAnalysis | 2.10.0 |
| MDTraj | 1.11.1 |
| NumPy | 2.4.3 |
| SciPy | 1.17.1 |
| Matplotlib | 3.10.8 |
| Seaborn | 0.13.2 |
| Pandas | 2.3.3 |
| Biopython | 1.87 |

Path: `/Users/zhangzhuo/miniforge3/envs/cgas-md/bin/python`

#### `py311` — Protein Docking

| Package | Version |
|---------|---------|
| Python | 3.11.15 |
| LightDock | 3.0 |
| Biopython | 1.87 |
| ProDy | 2.6.1 |
| FreeSASA | 2.2.1 |
| Cython | 3.2.4 |
| NumPy | 2.4.4 |
| SciPy | 1.17.1 |

Path: `/Users/zhangzhuo/miniforge3/envs/py311/bin/python`

### System Software

| Software | Version | Path |
|----------|---------|------|
| macOS | 15.3 (Darwin 25.3.0) | — |
| Apple Clang | 21.0.0 (clang-2100.0.123.102) | `/usr/bin/gcc` |
| System Python | 3.14.3 | `/opt/homebrew/bin/python3` |
| Homebrew | — | `/opt/homebrew` |
| FFTW | 3.3.11 | `/opt/homebrew/Cellar/fftw/3.3.11` (for SDOCK2.0, unused) |

## Hardware

| Component | Specification |
|-----------|--------------|
| CPU | Apple M3 Pro (12-core) |
| GPU | 18-core Apple Silicon GPU (via OpenCL) |
| Memory | 36 GB unified memory |
| OS | macOS 15.3 |
| Backend | OpenCL (Metal unavailable on conda-forge osx-arm64) |
| Performance | ~95 ns/day @ 80k atoms (domain-truncated) |

## Citation

Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. Science.
