# Scripts Directory

Organized by workflow stage. All scripts are runnable from the project root.

## Quick Start

```bash
# System build
python scripts/01_build/build_system.py --help

# Run production MD
python scripts/02_md/run_production.py --help

# Batch analysis
python scripts/03_analysis/batch_analyze.py --help

# MM-GBSA
python scripts/03_analysis/mmpbsa/run_batch.py --help
```

## Directory Layout

### `lib/` — Shared Library

| Module | Purpose |
|--------|---------|
| `lib/paths.py` | Project-wide path constants (`BASE`, `MD_DIR`, `ANALYSIS_DIR`, system definitions) |
| `lib/mda_tools.py` | MDAnalysis helpers (`load_universe`, `align_trajectory`, `get_com_distance`, `count_hbonds_simple`) |
| `lib/plot_style.py` | Matplotlib/Seaborn setup (`setup_mpl`, `savefig`, color constants) |
| `lib/stats.py` | Statistical utilities (`correlated_ttest`, `effective_sample_size`) |

### `01_build/` — System Building

| Script | Purpose |
|--------|---------|
| `build_system.py` | Build solvated MD system from PDB via tleap + OpenMM (ff19SB + OPC) |
| `build_from_docking.py` | Build system from Rosetta docking pose |
| `build_phosphorylated.py` | Build S305-phos (phosphoserine) system |
| `build_s305e.py` | Build S305E (glutamate mimic) system |
| `minimize.py` | Energy minimization with OpenMM |
| `mutate.py` | Apply point mutations to PDB |

### `02_md/` — MD Execution

| Script | Purpose |
|--------|---------|
| `run_production.py` | Production MD with OpenMM (heating → NPT → production) |
| `run_replicas.py` | Launch multiple replicas with different random seeds |
| `watch_and_launch.py` | Monitor and auto-restart crashed runs |

### `03_analysis/` — Trajectory Analysis

| Script | Purpose |
|--------|---------|
| `batch_analyze.py` | Fast batch analysis of multiple replicas (RMSD, COM, Rg, RMSF, contacts) |
| `deep_analysis.py` | **200ns deep analysis** (H-bonds timeline, PCA, DCCM, state classification) |
| `analyze_system.py` | Per-system analysis (RMSF, active site distances) |
| `analyze_trajectory.py` | Comprehensive single-trajectory analysis (RMSD, RMSF, H-bond, SASA, PCA) |
| `compare_systems.py` | Compare two systems (WT vs mutant) |
| `quick_compare_s305phos_wt.py` | Quick S305-phos vs WT comparison |
| `cluster.py` | GMM clustering of trajectory frames |
| `pca.py` | Principal component analysis |
| `dccm.py` | Dynamic cross-correlation matrix |

#### `03_analysis/mmpbsa/` — MM-GBSA

| Script | Purpose |
|--------|---------|
| `run_single.py` | Single-replica MM-GBSA |
| `run_batch.py` | Batch MM-GBSA with Delta binding energy (uses `-rp` + `-lp`) |

### `04_validation/` — Engine Validation

| Script | Purpose |
|--------|---------|
| `compare_engines.py` | **Final corrected** GROMACS vs OpenMM comparison (200ns, proper alignment) |
| `compare_engines_gmx2026.py` | GROMACS 2026 native vs OpenMM comparison |
| `convert_amber_to_gromacs.py` | Convert AMBER topology to GROMACS |
| `fix_gromacs_cmap.py` | Fix CMAP terms for AMBER force field in GROMACS |
| `verify_setup.py` | Verify OpenMM installation and GPU backend |

### `05_docking/` — Protein Docking

| Script | Purpose |
|--------|---------|
| `rosetta_dock.py` | Rosetta docking (4mut variant) |
| `rosetta_fastrelax_mutscan.py` | FastRelax + mutational scanning |
| `rosetta_mutational_scan.py` | Full mutational scan |
| `prepare_rosetta_input.py` | Prepare receptor + ligand PDBs |
| `prepare_4mut_input.py` | Prepare 4mut docking input |
| `analyze_poses.py` | Analyze LightDock top poses |
| `analyze_results.py` | Analyze docking results |

### `06_structure/` — Structure Analysis

| Script | Purpose |
|--------|---------|
| `analyze_af3.py` | Analyze AlphaFold3 predictions (ipTM, pLDDT, interface residues) |
| `process_af3.py` | Process raw AF3 output |
| `compare_predictions.py` | Compare Boltz-2 vs AF3 structures |
| `calc_distances.py` | Calculate CA-CA distances between active residues |

### Quaternary MVP Scripts

Located in `scripts/02_md/` and `scripts/03_analysis/`:

| Script | Purpose |
|--------|---------|
| `build_quaternary_mvp.py` | Assemble E2~Ub-TRIM41-cGAS four-component complex from PDB templates |
| `minimize_quaternary_mvp.py` | OpenMM energy minimization (5,000 steps, CUDA) |
| `run_quaternary_mvp.py` | Production MD runner (50 ns test, NVT→NPT, DCD output) |
| `analyze_quaternary_mvp.py` | Analyze K315→Ub distance, E2~Ub closed fraction, interface contacts |

### `archive/` — Archived Scripts

Not for routine use. Kept for reproducibility and reference.

| Subdir | Contents |
|--------|----------|
| `_versions/` | Old versions superseded by newer scripts (e.g. `compare_gmx_openmm_200ns_v1-v3`, `run_mmpbsa_batch_v1-v3`) |
| `_experiments/` | One-off exploratory scripts (umbrella sampling, mock trajectories, quick tests) |
| `_deprecated/` | Scripts no longer relevant (`analyze_c1_c4`, `analyze_close_state`) |

## Refactoring Notes

- **Old versions**: 15 versioned scripts archived in `_versions/` (saved ~4,000 lines of duplication in the active tree)
- **Library extraction**: `lib/` contains shared utilities extracted from 40+ scripts. Future scripts should `from lib.paths import ...` instead of redefining constants.
- **Remaining merge opportunities**:
  - `compare_engines.py` + `compare_engines_gmx2026.py` → unified engine comparison
  - `prepare_rosetta_input.py` + `prepare_4mut_input.py` → unified input preparation
  - `analyze_system.py` + `analyze_trajectory.py` → unified analysis core (pending careful merge to avoid breaking existing workflows)
