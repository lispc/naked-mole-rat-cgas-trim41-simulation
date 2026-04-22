# cGAS-TRIM41 Molecular Dynamics Study

> Computational investigation of how 4 amino acid variants in naked mole-rat cGAS affect TRIM41-mediated ubiquitination.
> 
> Based on: *A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging* (Chen et al., Science 2025)

## Project Structure

```
├── sequences/          # FASTA files for structure prediction
├── structures/         # Predicted structures (AF3)
│   └── af3_raw/        # Raw AF3 outputs
├── scripts/            # Analysis and simulation scripts
│   ├── analyze_af3.py       # AF3 structure quality analysis
│   ├── build_system.py      # MD system preparation (Amber/OpenMM)
│   ├── run_md.py            # Production MD with OpenMM
│   ├── analyze_trajectory.py # Trajectory analysis (RMSD, RMSF, etc.)
│   └── run_mmpbsa.py        # MM-GBSA/PBSA binding free energy
├── data/
│   ├── md_runs/        # MD system files and trajectories
│   ├── analysis/       # Analysis outputs
│   └── rosetta/        # Rosetta mutation scanning
├── docs/               # Documentation and planning
└── README.md           # This file
```

## Systems Simulated

| System | cGAS | TRIM41 | Mutations |
|--------|------|--------|-----------|
| Hsap_WT | Human WT | WT | — |
| Hsap_4mut | Human → NMR | WT | C463S, K479E, L495Y, K498T |
| Hgal_WT | NMR WT | WT | — |
| Hgal_4mut_rev | NMR → Human | WT | S463C, E511K, Y527L, T530K |

*Note: NMR = naked mole-rat (Heterocephalus glaber)*

## Quick Start

### Activate environment
```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md
```

### Analyze AF3 structures (after downloading results)
```bash
python scripts/analyze_af3.py
```

### Build MD system
```bash
python scripts/build_system.py --pdb structures/af3_raw/job1_Hsap_WT/ranked_0.pdb --name Hsap_WT
```

### Run production MD
```bash
python scripts/run_md.py --system-dir data/md_runs/Hsap_WT --name Hsap_WT_rep1 --seed 42
```

### Analyze trajectory
```bash
python scripts/analyze_trajectory.py \
    --topology data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --trajectory data/md_runs/Hsap_WT/production/Hsap_WT_rep1.dcd \
    --name Hsap_WT_rep1 \
    --highlight 463 479 495 498
```

### MM-GBSA binding energy
```bash
python scripts/run_mmpbsa.py \
    --topology data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --trajectory data/md_runs/Hsap_WT/production/Hsap_WT_rep1.dcd \
    --receptor-mask ":1-522" \
    --ligand-mask ":523-1152" \
    --name Hsap_WT_rep1
```

## Software Versions

- OpenMM 8.5.1
- AmberTools 24
- MDAnalysis 2.10.0
- MDTraj 1.11.1
- Python 3.13

## Hardware

- Apple M3 Pro (18-core GPU)
- 36 GB unified memory
- OpenCL backend (OpenMM)
