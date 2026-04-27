# GEMINI.md - cGAS-TRIM41 Molecular Dynamics Project

## Project Overview
This project is a computational investigation of the molecular mechanism by which four amino acid variants in naked mole-rat (NMR, *Heterocephalus glaber*) cGAS affect its interaction with the E3 ubiquitin ligase TRIM41. The study is based on the findings of Chen et al. (Science 2025). The core hypothesis is that these mutations, though distant from the physical interface, allosterically influence the geometric accessibility of cGAS lysine residues for TRIM41-mediated ubiquitination.

### Main Technologies
- **Simulation**: OpenMM 8.5.1 (CUDA/OpenCL)
- **Topology/System Prep**: AmberTools 24.8 (tleap), PDBFixer
- **Analysis**: MDAnalysis, NumPy, SciPy, Matplotlib
- **Protein Docking**: LightDock, Rosetta (PyRosetta 2025.06)
- **Structure Prediction**: AlphaFold3 (AF3)

## Project Structure
- `sequences/`: FASTA files for wild-type (WT) and mutant variants.
- `structures/`:
    - `af3_raw/`: Raw structure predictions from AlphaFold3.
    - `docking/`: Results from LightDock and Rosetta protein-protein docking.
- `scripts/`:
    - `build_system.py`: Prepares the molecular mechanics system using Amber ff19SB and OPC water.
    - `run_md.py`: Executes production molecular dynamics simulations.
    - `analyze_system.py`: Performs single-system analysis (RMSD, RMSF, contacts, distances).
    - `compare_systems.py`: Conducts cross-system statistical comparisons (Welch t-test, ΔRMSF).
    - `run_umbrella_sampling.py`: Launches Umbrella Sampling to calculate binding free energies (PMF).
- `data/`:
    - `md_runs/`: Topology (`.prmtop`), coordinates (`.rst7`), and trajectories (`.dcd`).
    - `analysis/`: JSON summaries and PNG plots from analysis scripts.
- `docs/`: Comprehensive project logs, docking reports, and theoretical analysis.

## Operational Workflows

### 1. Environment Setup
The project uses two primary Conda environments:
- `cgas-md`: For MD simulation and analysis.
- `rosetta`: For mutation scanning and docking verification.

### 2. Building a System
```bash
conda activate cgas-md
python scripts/build_system.py \
    --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
    --name Hgal_domain
```

### 3. Running Production MD
Simulations are typically run for 200ns per replica using CUDA.
```bash
CUDA_VISIBLE_DEVICES=0 python scripts/run_md.py \
    --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
    --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
    --name Hsap_WT_rep1 --outdir data/md_runs/Hsap_WT/rep1 \
    --prod-ns 200 --platform CUDA
```

### 4. Analysis Pipeline
1.  **Single System**: `python scripts/analyze_system.py --system Hgal_WT ...`
2.  **Comparison**: `python scripts/compare_systems.py --a <json_a> --b <json_b> ...`
3.  **Visualization**: Use PyMOL scripts in `scripts/` (e.g., `visualize_active_residues.pml`).

## Research Conventions
- **Coordinate Systems**: Paper numbering typically uses NMR cGAS coordinates (554 aa); human corresponds to 522 aa.
- **Interface Definition**: The physical interface is primarily at the cGAS N-terminus (resid 228–266).
- **Mutations (Corrected)**:
    - Hsap 4mut: D431S, K479E, L495Y, K498T.
    - Hgal 4mut_rev: S463D, E511K, Y527L, T530K.
- **Naming**: Systems are suffixed with `_WT`, `_4mut`, or `_4mut_rev`.

## Key Files & Documentation
- `README.md`: High-level summary and quick start.
- `docs/project_log.md`: The primary index for all research decisions and status updates.
- `docs/interface_analysis_report.md`: Detailed evidence for the allosteric mechanism.
- `docs/4mut_correction_log.md`: Critical record of mutation site mapping corrections.

## Ongoing Work (Status as of 2026-04-26)
- **Umbrella Sampling**: Calculating PMF for the TRIM41 RING to cGAS Lys-334 distance.
- **Comparison**: Finalizing the four-system comparison (Hgal WT/rev vs Hsap WT/mut).
- **Verification**: Validating the species-specific differences in domain dynamics.
