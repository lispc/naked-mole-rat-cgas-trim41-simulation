# cGAS-TRIM41 Molecular Dynamics Study

> Computational investigation of how four amino acid variants in naked mole-rat (*Heterocephalus glaber*) cGAS affect TRIM41-mediated ubiquitination.
>
> Based on: *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging* (Chen et al., Science 2025)

## Project Structure

```
├── AGENTS.md              # AI agent instructions
├── GEMINI.md              # Symlink to AGENTS.md
├── README.md              # This file
├── main.pdf               # Project manuscript
│
├── sequences/             # FASTA files for structure prediction
│   ├── Hsap_cGAS_WT.fasta
│   ├── Hsap_cGAS_4mut.fasta
│   ├── Hgal_cGAS_WT.fasta
│   ├── Hgal_cGAS_4mut_rev.fasta
│   └── TRIM41_WT.fasta
│
├── structures/            # Predicted structures and docking results
│   ├── af3_raw/           # Raw AlphaFold3 outputs
│   └── docking/           # Protein-protein docking results
│
├── scripts/               # Analysis and simulation scripts
│   ├── README.md          # Full script index
│   └── AGENTS.md          # Script-specific conventions
│
├── data/
│   ├── md_runs/           # OpenMM MD trajectories
│   │   ├── Hsap_WT/
│   │   ├── Hsap_4mut/
│   │   ├── Hgal_WT/
│   │   ├── Hgal_4mut_rev/
│   │   └── Hsap_WT_S305phos/
│   ├── md_runs_gmx/       # GROMACS MD trajectories (legacy conversion)
│   ├── md_runs_gmx2026/   # GROMACS 2026 native force field validation
│   └── analysis/          # Analysis outputs (JSON, PNG, NPZ)
│
├── figures/               # Publication-ready figures
│
└── docs/                  # Documentation (see index below)
    ├── README.md
    ├── 00-project/        # Project logs and meta-information
    ├── 10-reports/        # Completed experimental reports
    ├── 20-protocols/      # Protocols and execution plans
    ├── 30-diagnostics/    # Troubleshooting and divergence diagnosis
    ├── 40-reviews/        # External peer review feedback
    └── 50-infra/          # Hardware, software, environment
```

## Documentation Index

| Document | Content |
|----------|---------|
| [`docs/00-project/project_log.md`](docs/00-project/project_log.md) | **Primary log**: Current status, timeline, and decisions |
| [`docs/10-reports/docking_report.md`](docs/10-reports/docking_report.md) | Protein-protein docking: ClusPro, SDOCK2.0, LightDock, Rosetta |
| [`docs/10-reports/af3_report.md`](docs/10-reports/af3_report.md) | AlphaFold3 structure prediction: ipTM/pTM, mutation mapping |
| [`docs/10-reports/interface_analysis_report.md`](docs/10-reports/interface_analysis_report.md) | Interface analysis and allosteric mechanism evidence |
| [`docs/20-protocols/phosphorylation_md_plan.md`](docs/20-protocols/phosphorylation_md_plan.md) | Phosphorylation MD protocol and running records |
| [`docs/30-diagnostics/gromacs_openmm_divergence_diagnosis.md`](docs/30-diagnostics/gromacs_openmm_divergence_diagnosis.md) | GROMACS vs OpenMM divergence: CMAP fix, PBC pitfalls |
| [`docs/50-infra/hardware_benchmark.md`](docs/50-infra/hardware_benchmark.md) | Hardware benchmarks and MD performance (~152 ns/day) |
| [`docs/50-infra/software_versions.md`](docs/50-infra/software_versions.md) | Software version lockfile |

See [`docs/README.md`](docs/README.md) for the full index.

## Simulated Systems

All systems use the cGAS C-terminal domain construct (residues 200-554, 355 aa) in complex with TRIM41 SPRY, prepared with ff19SB + OPC water model.

| System | cGAS | TRIM41 | Mutations | Replicas | Status |
|--------|------|--------|-----------|----------|--------|
| **Hsap_WT** | Human WT | WT | — | 3 | Completed (3 × 200 ns) |
| **Hsap_4mut** | Human | WT | D431S, K479E, L495Y, K498T | 3 | Completed (3 × 200 ns) |
| **Hgal_WT** | NMR WT | WT | — | 3 | In progress (~183 / 200 ns rep1) |
| **Hgal_4mut_rev** | NMR | WT | S463D, E511K, Y527L, T530K | 3 | In progress (~84 / 200 ns rep1) |
| **Hsap_WT_S305phos** | Human WT, SEP@305 | WT | — | 3 | Completed (3 × 200 ns) |
| **Hsap_WT_S305E** | Human WT, S305E | WT | — | 3 | Completed (3 × 200 ns) |
| **Quaternary MVP (WT)** | Human WT + TRIM25 RING proxy | UBE2D1~Ub | — | 1 | In progress (0.8 / 50 ns) |
| **Quaternary MVP (4mut)** | Human 4mut + TRIM25 RING proxy | UBE2D1~Ub | D431S, K479E, L495Y, K498T | 1 | In progress (0.3 / 50 ns) |

*NMR = naked mole-rat (*Heterocephalus glaber*). Paper numbering uses NMR coordinates (554 aa); human corresponds to 522 aa.*

**Cross-engine validation:**
- GROMACS 2026 native `amber19sb.ff` + OPC: Hsap_WT rep1, ~83 ns / 200 ns (ongoing). CMAP fix validated: COM/Rg match OpenMM within 1%, RMSD ratio 1.34×.

**Boltz-2 validation:**
- Local Boltz-2 (v2.2.1) installed in `boltz` conda environment. Full-length cGAS-TRIM41: ipTM 0.17 (identical to AF3). Domain-truncated: ipTM 0.33, but TRIM41 placement qualitatively differs from AF3 (RMSD 21.34 Å, Jaccard 0.00). Confirms docking+MD approach is correct for this transient complex.

## Key Findings

### 1. Mutations are NOT at the physical interface

The four mutation sites (residues 463, 511, 527, 530 in NMR numbering) are **>200 residues away** in sequence from the interface region (residues 228-266) and **30-39 Å away** in 3D space. They cannot directly contact TRIM41.

| Mutation | NMR Resid | Distance to nearest interface residue | Distance to TRIM41 |
|----------|-----------|--------------------------------------|-------------------|
| S463 | 482 | 32 residues | ~29 Å |
| E511 | 530 | 80 residues | ~31 Å |
| Y527 | 546 | 96 residues | ~31 Å |
| T530 | 549 | 99 residues | ~39 Å |

### 2. Interface is at the cGAS N-terminus

The physical binding interface involves cGAS residues **228-266** and TRIM41 residues **82-190**:
- TRIM41-157/158 ↔ cGAS-258/259 (occupancy >60%, most stable)
- TRIM41-187 ↔ cGAS-236 (occupancy ~39%)

### 3. Human systems are more dynamic than naked mole-rat

| Metric | Hgal_WT | Hsap_WT | Hsap_4mut |
|--------|---------|---------|-----------|
| RMSD | 5.16 ± 0.72 Å | **8.94 ± 1.58 Å** | **9.76 ± 2.21 Å** |
| COM | 37.4 ± 0.9 Å | **46.6 ± 2.4 Å** | **49.0 ± 2.8 Å** |

The 4mut makes the human complex **even less stable** (RMSD +9%, COM +5%, p < 1e-30).

### 4. Active site geometry diverges between species

| Site | Hgal_WT | Hgal_4mut_rev | Hsap_WT | Hsap_4mut |
|------|---------|---------------|---------|-----------|
| S463/D431 | 29.4 Å | **18.6 Å** (−10.8) | 29.4 Å | 32.0 Å (+2.6) |
| E511/K479 | 31.5 Å | **21.9 Å** (−9.7) | 36.6 Å | 40.0 Å (+3.4) |

**Hgal_4mut_rev** brings active sites **closer** to TRIM41; **Hsap_4mut** pushes them **farther**.

### 5. Lys-334 is the top ubiquitination candidate in human cGAS

Lysine accessibility analysis reveals **Lys-334** as the closest cGAS lysine to the TRIM41 RING domain:
- Hsap WT: 10.4 Å
- Hsap 4mut: **6.4 Å** (−4.0 Å)

The 4mut may **concentrate** ubiquitination to this single site.

### 6. S305 phosphorylation causes complex dissociation (unexpected)

Preliminary analysis of S305-phos (~138 ns) shows **dissociation** in all three replicas:

| System | COM (Å) | Rg (Å) | Behavior |
|--------|---------|--------|----------|
| WT (reference) | 45.1 ± 3.2 | 31.1 ± 1.3 | Stable bound |
| S305-phos rep1 | **67.6 ± 1.4** | 38.9 ± 0.6 | Dissociating |
| S305-phos rep2 | **89.8 ± 10.0** | 48.4 ± 4.3 | Fully dissociated |
| S305-phos rep3 | **71.3 ± 3.4** | 40.2 ± 1.4 | Dissociating |

> Note: The 3 replicas are at ~140 ns / 200 ns as of 2026-04-23.

This contradicts Zhen et al. (2023), who reported that CHK2 phosphorylation at S305 **promotes** cGAS-TRIM41 binding. Possible explanations: (1) force field overestimation of −2 charge repulsion in solution, (2) nuclear/DNA environment required for stabilization, (3) dissociation is an intermediate state prior to DNA-mediated recruitment.

### Allosteric mechanism hypothesis

Since the mutations are far from the interface, their effect must be **allosteric**:

```
C-terminal 4mut  →  altered global dynamics  →  changed RING-to-Lys geometry  →  differential ubiquitination
```

## Quick Start

### Environment Setup

```bash
# Activate the MD environment
conda activate cgas-md
python --version  # Python 3.11

# For GROMACS
conda activate gmx
gmx --version  # GROMACS 2026.0
```

### Build an MD System

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
    --name Hsap_WT_rep1 \
    --outdir data/md_runs/Hsap_WT/rep1 \
    --prod-ns 200 \
    --platform CUDA
```

**Multi-GPU replicas:** Use `CUDA_VISIBLE_DEVICES` to isolate GPUs. The script hardcodes `CudaDeviceIndex='0'`.

### Analyze a Single System

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

### Compare Two Systems

```bash
conda activate cgas-md
python scripts/compare_systems.py \
    --a data/analysis/run_A/Hsap_WT_summary.json \
    --b data/analysis/run_B/Hsap_4mut_summary.json \
    --name-a Hsap_WT --name-b Hsap_4mut \
    --outdir data/analysis/comparison
```

## Environment

### Conda Environments

| Environment | Purpose | Key Packages |
|-------------|---------|--------------|
| `cgas-md` | OpenMM MD, AmberTools, analysis | OpenMM 8.5.1, AmberTools 24.8, MDAnalysis 2.10.0, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.8 |
| `gmx` | GROMACS 2026.0 CUDA | GROMACS 2026.0 (nompi_cuda) |
| `rosetta` | PyRosetta docking | PyRosetta 2025.06 |
| `boltz` | **Boltz-2 structure prediction** | **Boltz-2.2.1, PyTorch 2.11.0, CUDA 13** |

**Always use conda environments. Do not use the system Python.**

### Hardware

| Component | Specification |
|-----------|--------------|
| CPU | AMD EPYC 7702, 64 cores |
| GPU | 4× NVIDIA RTX 3090 (24 GB VRAM each) |
| OS | Linux x86_64, CUDA 13.0 |
| Performance | ~152 ns/day @ 116k atoms (OpenMM, domain-truncated) |

## Citation

Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. **Science**.
