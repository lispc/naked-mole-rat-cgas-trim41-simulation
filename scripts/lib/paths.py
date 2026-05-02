"""Project-wide path constants."""
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")

# Data directories
DATA_DIR = BASE / "data"
MD_DIR = DATA_DIR / "md_runs"
MD_DIR_GMX = DATA_DIR / "md_runs_gmx2026"
ANALYSIS_DIR = DATA_DIR / "analysis"
STRUCTURES_DIR = BASE / "structures"
SEQUENCES_DIR = BASE / "sequences"
FIGURES_DIR = BASE / "figures"
DOCS_DIR = BASE / "docs"

# Specific analysis subdirectories
MMPBSA_DIR = ANALYSIS_DIR / "mmpbsa"
FINAL_200NS_DIR = ANALYSIS_DIR / "final_200ns"
GMX_OPENMM_DIR = ANALYSIS_DIR / "gmx_openmm_comparison"
S305PHOS_VS_WT_DIR = ANALYSIS_DIR / "s305phos_vs_wt"
DEEP_200NS_DIR = ANALYSIS_DIR / "deep_200ns"

# System definitions
SYSTEMS = {
    "Hsap_WT": {
        "prmtop": MD_DIR / "Hsap_WT" / "Hsap_WT.prmtop",
        "protein_prmtop": MD_DIR / "Hsap_WT" / "Hsap_WT_protein.prmtop",
        "receptor_prmtop": MD_DIR / "Hsap_WT" / "Hsap_WT_receptor.prmtop",
        "ligand_prmtop": MD_DIR / "Hsap_WT" / "Hsap_WT_ligand.prmtop",
        "cgas_range": "1-218",
        "trim_range": "219-541",
        "reps": [MD_DIR / "Hsap_WT" / f"rep{i}" / f"Hsap_WT_rep{i}_prod.dcd" for i in (1, 2, 3)],
    },
    "Hsap_WT_S305phos": {
        "prmtop": MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos.prmtop",
        "protein_prmtop": MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos_protein.prmtop",
        "receptor_prmtop": MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos_receptor.prmtop",
        "ligand_prmtop": MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos_ligand.prmtop",
        "cgas_range": "1-218",
        "trim_range": "219-541",
        "reps": [MD_DIR / "Hsap_WT_S305phos" / f"rep{i}" / f"Hsap_WT_S305phos_rep{i}_prod.dcd" for i in (1, 2, 3)],
    },
    "Hsap_WT_S305E": {
        "prmtop": MD_DIR / "Hsap_WT_S305E" / "Hsap_WT_S305E.prmtop",
        "cgas_range": "1-218",
        "trim_range": "219-541",
        "reps": [MD_DIR / "Hsap_WT_S305E" / f"rep{i}" / f"Hsap_WT_S305E_rep{i}_prod.dcd" for i in (1, 2, 3)],
    },
}
