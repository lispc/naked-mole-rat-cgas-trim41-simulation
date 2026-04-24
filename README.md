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
│   ├── analyze_system.py           # Generic single-system MD analysis (NEW)
│   ├── compare_systems.py          # Cross-system comparison (NEW)
│   ├── generate_mock_trajectory.py # Mock trajectory generator (NEW)
│   ├── analyze_trajectory.py       # Legacy trajectory analysis
│   └── run_mmpbsa.py               # MM-GBSA binding energy
├── data/
│   ├── md_runs/        # MD system files and trajectories
│   │   ├── Hgal_domain/rep{1,2,3}/  # Hgal WT 3×200ns (COMPLETED)
│   │   ├── Hsap_WT/rep1/            # Hsap WT 1×200ns (RUNNING)
│   │   ├── Hsap_4mut/rep1/          # Hsap 4mut 1×200ns (RUNNING)
│   │   └── Hgal_4mut_rev/rep1/      # Hgal 4mut_rev 1×200ns (RUNNING)
│   ├── analysis/       # Analysis outputs
│   │   ├── final_200ns/             # Hgal WT final analysis
│   │   └── mock_run/                # Pipeline dry-run outputs
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
| [`docs/af3_report.md`](docs/af3_report.md) | **AF3 结构预测**：4 个 job 的 ipTM/pTM 数据、突变映射表 |
| [`docs/hardware_benchmark.md`](docs/hardware_benchmark.md) | **硬件与性能**：4×RTX 3090 基准测试、MD 速度估算（~152 ns/day @ 116k atoms） |
| [`docs/execution_plan_v1.md`](docs/execution_plan_v1.md) | **执行方案 v1.0**：6 个阶段的详细计划、时间线 |
| [`docs/computational_workflow.md`](docs/computational_workflow.md) | **原始方案设计**：方法学选择理由、工具对比、技术路线 |
| [`docs/paper_notes_cgas_trim41.md`](docs/paper_notes_cgas_trim41.md) | **论文理解笔记**：Chen et al. Science 2025 的关键发现、实验因果链、突变信息 |
| [`docs/cluspro_submission_guide.md`](docs/cluspro_submission_guide.md) | **ClusPro 提交指南**：6 个 job 的详细参数、活性残基列表 |

## Systems Simulated

| System | cGAS | TRIM41 | Mutations | Status |
|--------|------|--------|-----------|--------|
| Hgal_WT | NMR WT | WT | — | ✅ 3×200ns 完成 + 分析 |
| Hsap_WT | Human WT | WT | — | 🔄 1×200ns 运行中 |
| Hsap_4mut | Human → NMR | WT | C463S, K479E, L495Y, K498T | 🔄 1×200ns 运行中 |
| Hgal_4mut_rev | NMR → Human | WT | S463C, E511K, Y527L, T530K | 🔄 1×200ns 运行中 |

*NMR = naked mole-rat (Heterocephalus glaber). Paper numbering uses NMR coordinates.*

## Key Findings

### Finding 1: Mutations are NOT at the physical interface

**MD evidence** (Hgal WT 3×200ns):

| Mutation site | PDB resid | Distance to nearest interface residue | Distance to TRIM41 |
|--------------|-----------|--------------------------------------|-------------------|
| S463 | 482 | **32 residues** | ~29 Å |
| E511 | 530 | **80 residues** | ~31 Å |
| Y527 | 546 | **96 residues** | ~31 Å |
| T530 | 549 | **99 residues** | ~39 Å |

The 4 mutation sites (resid 482–549) are **>200 residues away** in sequence from the interface region (resid 228–266) and **30–39 Å away in space**. They cannot directly contact TRIM41.

### Finding 2: Interface is at the cGAS N-terminus

Physical interface involves cGAS resid **228–266** ↔ TRIM41 82–190:
- TRIM41-157/158 ↔ cGAS-258/259 (occupancy >60%, most stable)
- TRIM41-187 ↔ cGAS-236 (occupancy ~39%)

The C-terminal mutation sites are spatially separated from the interface by the entire cGAS domain.

### Finding 3: Allosteric mechanism hypothesis

Since mutations are far from the interface, their effect on TRIM41 recognition must be **allosteric**:
- C-terminal 4mut → N-terminal ~211–219 rearrangement (up to 12 Å displacement, per AF3 monomer comparison)
- → Interface micro-change → altered TRIM41 recognition/ubiquitination efficiency

**Direct test pending**: Hgal WT vs Hgal 4mut_rev MD comparison (expected completion: ~April 26).

## Quick Start

### Environment Setup

```bash
# Conda base
export CONDA_PATH="$HOME/miniforge3"
source $CONDA_PATH/etc/profile.d/conda.sh

# MD environment (OpenMM, AmberTools, analysis)
conda activate cgas-md
python --version  # Python 3.11.15
```

### Build MD System

```bash
conda activate cgas-md
python scripts/build_system.py \
    --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
    --name Hgal_domain
```

### Run Production MD

```bash
conda activate cgas-md
CUDA_VISIBLE_DEVICES=0 python scripts/run_md.py \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
    --name Hsap_WT_rep1 --outdir data/md_runs/Hsap_WT/rep1 \
    --prod-ns 200 --platform CUDA
```

### Analyze Single System

```bash
conda activate cgas-md
python scripts/analyze_system.py \
    --system Hgal_WT \
    --prmtop data/md_runs/Hgal_domain/Hgal_domain.prmtop \
    --trajectories data/md_runs/Hgal_domain/rep1/*.dcd \
    --replica-names rep1 \
    --cgas-range 219 573 \
    --active-sites '{"S463": 482, "E511": 530, "Y527": 546, "T530": 549}' \
    --dt-ns 0.1 \
    --outdir data/analysis/my_run
```

### Compare Two Systems

```bash
conda activate cgas-md
python scripts/compare_systems.py \
    --a data/analysis/run_A/Hgal_WT_summary.json \
    --b data/analysis/run_B/Hgal_4mut_rev_summary.json \
    --name-a Hgal_WT --name-b Hgal_4mut_rev \
    --outdir data/analysis/comparison
```

### Generate Mock Trajectory

```bash
conda activate cgas-md
python scripts/generate_mock_trajectory.py \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --coord data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
    --out-dcd data/md_runs/Hsap_WT/mock/mock_prod.dcd \
    --out-log data/md_runs/Hsap_WT/mock/mock_prod.log \
    --duration-ps 10
```

## Analysis Pipeline

```
MD trajectories → analyze_system.py → per-system JSON + PNG plots
                                    ↓
                         compare_systems.py → cross-system comparison
                                              (Welch t-test, ΔRMSF,
                                               contact differences)
```

**Metrics computed**:
- RMSD time series (per replica)
- RMSF per residue (cross-replica overlay)
- Interface contact occupancy (heatmap + top list)
- COM distance time series
- Active site distances (optional)

## Environment Details

### Conda Environments

#### `cgas-md` — MD Simulation & Analysis

| Package | Version |
|---------|---------|
| Python | 3.11.15 |
| OpenMM | 8.5.1 (conda-forge, CUDA build) |
| AmberTools | 24.8 |
| MDAnalysis | 2.10.0 |
| SciPy | 1.17.1 |
| Matplotlib | 3.10.8 |
| NumPy | 2.4.3 |

Path: `~/miniforge3/envs/cgas-md/bin/python`

## Hardware

| Component | Specification |
|-----------|--------------|
| CPU | x86_64 server CPU |
| GPU | 4× NVIDIA GeForce RTX 3090 (24GB VRAM each) |
| Memory | Server RAM (充足) |
| OS | Linux (CUDA 13.0) |
| Backend | CUDA (mixed precision) |
| Performance | ~152 ns/day @ 116k atoms (domain-truncated, measured) |

## Citation

Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. Science.
