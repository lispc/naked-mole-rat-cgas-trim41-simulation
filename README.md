# cGAS-TRIM41 Molecular Dynamics Study

> Computational investigation of how four amino acid variants affect cGAS-TRIM41 interaction and TRIM41-mediated ubiquitination.
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
│   ├── Hsap_cGAS_4mut.fasta           # ✅ Correct: D431S, K479E, L495Y, K498T
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
│   │   └── Hsap_WT_S305E/
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
| [`docs/10-reports/allosteric_mechanism_analysis.md`](docs/10-reports/allosteric_mechanism_analysis.md) | **Latest**: Chai-1 vs AF3 allostery + three-mechanism analysis |
| [`docs/10-reports/af3_mutation_analysis.md`](docs/10-reports/af3_mutation_analysis.md) | AF3 structure prediction: 4mut does NOT change active-site geometry |
| [`docs/10-reports/interface_analysis_report.md`](docs/10-reports/interface_analysis_report.md) | Interface analysis and allosteric mechanism evidence |
| [`docs/10-reports/rosetta_mutational_scan_report.md`](docs/10-reports/rosetta_mutational_scan_report.md) | Rosetta scan: I_sc unchanged, allosteric effect confirmed |
| [`docs/10-reports/docking_report.md`](docs/10-reports/docking_report.md) | Protein-protein docking: ClusPro, SDOCK2.0, LightDock, Rosetta |
| [`docs/20-protocols/phosphorylation_md_plan.md`](docs/20-protocols/phosphorylation_md_plan.md) | Phosphorylation MD protocol and running records |
| [`docs/30-diagnostics/md_realism_gap.md`](docs/30-diagnostics/md_realism_gap.md) | Gap analysis: simulation vs experimental realism |
| [`docs/50-infra/hardware_benchmark.md`](docs/50-infra/hardware_benchmark.md) | Hardware benchmarks and MD performance |

See [`docs/README.md`](docs/README.md) for the full index.

## Simulated Systems

All binary systems use the cGAS C-terminal domain construct in complex with TRIM41 SPRY, prepared with ff19SB + OPC water model.

| System | cGAS | TRIM41 | Mutations | Replicas | Status |
|--------|------|--------|-----------|----------|--------|
| **Hsap_WT** | Human WT | WT | — | 3 | Completed (3 × 200 ns) |
| **Hsap_4mut** | Human | WT | D431S, K479E, L495Y, K498T | 3 | Completed (3 × 200 ns) |
| **Hgal_WT** | NMR WT | WT | — | 3 | Completed (rep1 ~104 ns, rep2/rep3 ~96 ns) |
| **Hgal_4mut_rev** | NMR | WT | S463D, E511K, Y527L, T530K | 3 | Completed (rep1 ~104 ns, rep2/rep3 ~61 ns) |
| **Hsap_WT_S305E** | Human WT, S305E | WT | — | 3 | Completed (3 × 200 ns) |
| **Hsap_WT_S305phos** | Human WT, SEP@305 | WT | — | 3 | Completed (3 × 200 ns) |

*NMR = naked mole-rat (*Heterocephalus glaber*).*

**⚠️ Known Issues:**
- ~~`Hgal_WT_rep2_prod.dcd` is corrupted (zero magic bytes)~~ ✅ **FIXED** (2026-05-03): Repaired via binary frame extraction + mdtraj rewrite; simulation restarted from 111 ns checkpoint
- Hsap systems originally used incorrect residue C463S; corrected to D431S in later analyses. See [`docs/00-project/4mut_correction_log.md`](docs/00-project/4mut_correction_log.md)
- ~~All "Lys-334" references in US code/docs should read **"Lys-315"** (prmtop renumbering artifact)~~ ✅ **FIXED**: Active codebase (build/analysis scripts) all use K315; only archived deprecated scripts retain old label. See [`docs/reviews/0427-response.md`](docs/reviews/0427-response.md)

---

## Key Findings (Updated 2026-04-23)

### 1. Mutations are NOT at the physical interface

The four mutation sites are **>20 Å away** from the cGAS-SPRY interface in all prediction methods (Rosetta, Chai-1, AF3). They cannot directly contact TRIM41.

| Evidence | Distance to SPRY | Interpretation |
|----------|-----------------|----------------|
| Rosetta I_sc | I_sc unchanged (1e-12 level) | No direct physical contact |
| Chai-1 truncated | 24–39 Å | Far from interface |
| AF3 full-length | ~7 Å (but see caveat below) | May be full-length conformational artifact |

**Caveat**: AF3 full-length places D431/L495 ~7 Å from SPRY, but Chai-1 truncated (same sequence, no NTD) places them at 24–32 Å. The AF3 "short-range" positioning may be an NTD-constraint artifact. Chai-1 is more consistent with Rosetta energy evidence.

### 2. The allosteric mechanism is DYNAMIC, not structural

**Static structure analyses** (AF3 monomer, Chai-1, Rosetta scan) show:
- 4mut site RMSD < 0.6 Å (local structure unchanged)
- DNA-binding groove SASA unchanged (Δ = 0%)
- Catalytic site Rg unchanged (Δ < 0.1 Å)

**Dynamic analyses** from MD trajectories show:
- **DCCM**: WT shows anti-correlated motion between N-terminal (305-307) and C-terminal 4mut region (474-478, corr ≈ −0.44); 4mut_rev **eliminates** this anti-correlation
- **ΔRMSF**: Hgal 4mut_rev makes the N-terminal region significantly more flexible (80% of residues with |ΔRMSF| > 0.5 Å) and the C-terminal region more rigid (mean ΔRMSF = −6.4 Å)
- **PCA**: PC1 captures 72% of variance; WT and 4mut_rev occupy distinct conformational subspaces

**Conclusion**: 4mut acts as a **dynamic allosteric switch**—it does not change what cGAS looks like (static structure), but it changes how cGAS moves (dynamic coupling network).

### 3. Species differences are LARGE

| Metric | Hgal_WT | Hsap_WT | Hsap_4mut |
|--------|---------|---------|-----------|
| Global RMSD (species) | — | 20.9 Å vs Hgal | — |
| Interface location | N-terminal (203-249) | Mid-domain (288-315) | — |

**Hgal and Hsap cGAS have fundamentally different folds.** Hgal MD results cannot be directly extrapolated to Hsap.

### 4. S305 phosphorylation causes complex dissociation (unexpected)

S305-phos (SEP, −2 charge) causes **complete dissociation** in <100 ns. S305E (−1 charge, Glu) causes **weakened binding** (COM −3.7 Å, H-bonds −2.0) but **not dissociation**.

This contradicts Zhen et al. (2023), who reported phosphorylation **promotes** binding. Possible explanation: charge density matters (−2 phosphate vs −1 Glu).

### 5. Revised mechanism hypothesis

```
C-terminal 4mut
    ↓  Local structure unchanged (<0.6 Å)
Dynamic coupling network重塑
    ↓  DCCM anti-correlation eliminated; RMSF pattern changed
N-terminal conformational sampling altered
    ↓  Interface dynamics changed
TRIM41 recognition / ubiquitination efficiency altered
```

**Evidence strength**:
- ⭐⭐⭐⭐⭐ Rosetta I_sc unchanged (no direct contact)
- ⭐⭐⭐⭐⭐ AF3 N-terminal 12 Å displacement (long-range structural effect)
- ⭐⭐⭐⭐ Chai-1 WT vs 4mut interface difference marginal (binding mode unchanged)
- ⭐⭐⭐ DCCM anti-correlation change (dynamic coupling altered)
- ⭐⭐⭐ ΔRMSF / PCA (conformational subspace separation)

---

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

**Multi-GPU replicas:** Use `CUDA_VISIBLE_DEVICES` to isolate GPUs.

### Analyze a Single System

```bash
conda activate cgas-md
python scripts/analyze_system.py \
    --system Hsap_WT \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --trajectories data/md_runs/Hsap_WT/rep1/*.dcd \
    --replica-names rep1 \
    --cgas-range 219 541 \
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

---

## Environment

### Conda Environments

| Environment | Purpose | Key Packages |
|-------------|---------|--------------|
| `cgas-md` | OpenMM MD, AmberTools, analysis | OpenMM 8.5.1, AmberTools 24.8, MDAnalysis 2.10.0, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.8 |
| `gmx` | GROMACS 2026.0 CUDA | GROMACS 2026.0 (nompi_cuda) |
| `rosetta` | PyRosetta docking | PyRosetta 2025.06 |
| `boltz` | Boltz-2 structure prediction | Boltz-2.2.1, PyTorch 2.11.0, CUDA 13 |

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
