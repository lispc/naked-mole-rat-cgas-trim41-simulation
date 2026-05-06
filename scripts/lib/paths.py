"""Project-wide path constants. Auto-detects project root from this file's location."""
from pathlib import Path

# Auto-detect: this file is at <project>/scripts/lib/paths.py
BASE = Path(__file__).resolve().parent.parent.parent

# Top-level directories
DATA_DIR = BASE / "data"
SCRIPTS_DIR = BASE / "scripts"
STRUCTURES_DIR = BASE / "structures"
SEQUENCES_DIR = BASE / "sequences"
FIGURES_DIR = BASE / "figures"
DOCS_DIR = BASE / "docs"

# Data subdirectories
MD_DIR = DATA_DIR / "md_runs"
MD_DIR_GMX = DATA_DIR / "md_runs_gmx2026"
ANALYSIS_DIR = DATA_DIR / "analysis"
STRUCTURES_DATA_DIR = DATA_DIR / "structures"

# Key analysis subdirectories
MMPBSA_DIR = ANALYSIS_DIR / "mmpbsa"
FOUR_SYSTEM_DIR = ANALYSIS_DIR / "four_system"
GMX_OPENMM_DIR = ANALYSIS_DIR / "gmx_openmm_comparison"
HSAP_BATCH_DIR = ANALYSIS_DIR / "hsap_batch"
ALLOSTERIC_DIR = ANALYSIS_DIR / "allosteric_network"
DELTA_RMSF_DIR = ANALYSIS_DIR / "delta_rmsf"
DELTA_DCCM_DIR = ANALYSIS_DIR / "delta_dccm"
QUATERNARY_MINIMAL_DIR = MD_DIR / "quaternary_minimal"
QUATERNARY_MVP_DIR = MD_DIR / "quaternary_mvp"

# System definitions with full paths
SYSTEMS = {
    "Hsap_WT": {
        "prmtop": MD_DIR / "Hsap_WT" / "Hsap_WT.prmtop",
        "cgas_last": 198,
        "reps": [MD_DIR / "Hsap_WT" / f"rep{i}" / f"Hsap_WT_rep{i}_prod.dcd" for i in (1, 2, 3)],
    },
    "Hsap_4mut": {
        "prmtop": MD_DIR / "Hsap_4mut" / "Hsap_4mut.prmtop",
        "cgas_last": 198,
        "reps": [MD_DIR / "Hsap_4mut" / f"rep{i}" / f"Hsap_4mut_rep{i}_prod.dcd" for i in (1, 2, 3)],
    },
    "Hgal_WT": {
        "prmtop": MD_DIR / "Hgal_WT" / "rep1" / "Hgal_WT.prmtop",
        "cgas_last": 218,
        "reps": [
            MD_DIR / "Hgal_WT" / "rep1" / "Hgal_WT_rep1_prod.dcd",
            MD_DIR / "Hgal_WT" / "rep2" / "Hgal_WT_rep2_prod_fixed.dcd",
            MD_DIR / "Hgal_WT" / "rep2" / "Hgal_WT_rep2_restart.dcd",
            MD_DIR / "Hgal_WT" / "rep3" / "Hgal_WT_rep3_prod.dcd",
        ],
    },
    "Hgal_4mut_rev": {
        "prmtop": MD_DIR / "Hgal_4mut_rev" / "rep1" / "Hgal_4mut_rev.prmtop",
        "cgas_last": 218,
        "reps": [
            MD_DIR / "Hgal_4mut_rev" / "rep1" / "Hgal_4mut_rev_rep1_prod.dcd",
            MD_DIR / "Hgal_4mut_rev" / "rep2" / "Hgal_4mut_rev_rep2_prod.dcd",
            MD_DIR / "Hgal_4mut_rev" / "rep2" / "Hgal_4mut_rev_rep2_restart.dcd",
            MD_DIR / "Hgal_4mut_rev" / "rep3" / "Hgal_4mut_rev_rep3_prod.dcd",
        ],
    },
}
