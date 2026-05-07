# cGAS-TRIM41 Molecular Dynamics Study

> Computational investigation of how four naked mole-rat cGAS variants modulate TRIM41-mediated ubiquitination via long-range allostery.
>
> Based on: Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. **Science**.

---

## Primary Outputs

| Document | Description |
|----------|-------------|
| 📄 **[`paper/paper_final.md`](paper/paper_final.md)** | **Final paper** — core conclusions and supporting logic |
| 📋 **[`docs/00-project/lab_report.md`](docs/00-project/lab_report.md)** | **Lab report** — full experimental record, decisions, and lessons learned |
| 📝 **[`docs/00-project/retrospective.md`](docs/00-project/retrospective.md)** | **Retrospective** — what went well, what didn't, and key takeaways |

---

## Key Findings

1. **The four variant sites (D431S/K479E/L495Y/K498T) do NOT contact TRIM41.** They are 30–39 Å from the binding interface. The physical interface is in the N-terminal region.

2. **4mut induces long-range allosteric conformational shifts.** The N-terminal interface region moves up to 12.3 Å while the variant sites move <0.7 Å (AF3 monomer comparison).

3. **Binding affinity is unchanged.** MD geometric metrics (COM, contacts, Rg) are unaltered in human cGAS. MM-GBSA ΔΔG = 4.5 ± 10.1 kcal/mol (p = 0.50).

4. **4mut reshapes conformational dynamics.** DCCM shows altered long-range coupling. PCA reveals distinct but overlapping conformational spaces. ΔRMSF shows no globally significant flexibility changes.

5. **K315 is the only geometrically accessible ubiquitination target.** A full scan of all 36 cGAS lysines shows K315 is the only residue within 25 Å of the Ub-G76 catalytic center.

6. **Umbrella sampling PMF: 4mut shifts the K315 free energy minimum by −2.5 Å** (WT: 21.5 Å → 4mut: 19.0 Å) and reduces the free energy cost to reach 15 Å by 1.1 kcal/mol. Catalytic readiness (K315<19 Å + correct angle): **WT 12.6% → 4mut 53.7% (4.3×).**

---

## Project Structure

```
├── paper/
│   └── paper_final.md       # Final manuscript
├── sequences/               # FASTA files for all systems
├── structures/              # AF3 predictions, docking results, PDB templates
├── scripts/                 # Build, MD, analysis, docking, validation
│   ├── lib/                 # Shared library (paths, stats, plotting, MDA tools)
│   └── archive/             # Deprecated and experimental scripts
├── data/
│   ├── md_runs/             # MD trajectories (binary systems + quaternary)
│   ├── md_runs_gmx2026/     # GROMACS cross-validation
│   ├── analysis/            # Analysis outputs (NPZ, JSON, PNG)
│   └── structures/          # Built system topologies
├── figures/                 # Publication-ready figures
└── docs/                    # Full documentation (see docs/README.md)
```

---

## Simulated Systems

| System | cGAS | TRIM41 | Replicas | Status |
|--------|------|--------|----------|--------|
| Hsap_WT | Human WT | SPRY | 3 × 200 ns | ✅ |
| Hsap_4mut | Human 4mut | SPRY | 3 × 200 ns | ✅ |
| Hgal_WT | NMR WT | SPRY | 3 × 200 ns | ✅ |
| Hgal_4mut_rev | NMR rev | SPRY | 3 × 200 ns | ✅ |
| Quaternary WT | Human WT | RING+SPRY+E2~Ub | 1 × 50 ns | ✅ |
| Quaternary 4mut | Human 4mut | RING+SPRY+E2~Ub | 1 × 50 ns | ✅ |
| US K315 (WT) | — | RING+SPRY+E2~Ub | 11 windows × 10 ns | ✅ |
| US K315 (4mut) | — | RING+SPRY+E2~Ub | 10 windows × 10 ns | ✅ |

---

## Quick Start

```bash
conda activate cgas-md

# Run production MD
CUDA_VISIBLE_DEVICES=0 python scripts/02_md/run_production.py \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
    --name Hsap_WT_rep1 --outdir data/md_runs/Hsap_WT/rep1 \
    --prod-ns 200 --platform CUDA

# Four-system comparison
python scripts/03_analysis/compare_four_systems_fast.py

# WHAM PMF analysis
python scripts/03_analysis/wham_pmf.py --system WT --windows 12,...,22
```

---

## Environment

| Component | Specification |
|-----------|--------------|
| CPU | AMD EPYC 7702, 64 cores |
| GPU | 4× NVIDIA RTX 3090 (24 GB) |
| MD Engine | OpenMM 8.5.1 + GROMACS 2026.0 |
| Force Field | Amber ff19SB + OPC water |
| Analysis | MDAnalysis 2.10.0, AmberTools 24.8 |

## Citation

Chen et al. (2025). *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*. **Science**.
