#!/usr/bin/env python3
"""
Automated processing of AlphaFold3 Server results.
Unzips, extracts monomers, fixes with PDBFixer, generates analysis plots.

Usage:
  python scripts/process_af3_results.py --job-dir structures/af3_raw/job1_Hsap_WT
"""
import argparse
import json
import subprocess
import zipfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def unzip_results(job_dir: Path):
    """Unzip AF3 result archive."""
    zips = list(job_dir.glob("*.zip"))
    if not zips:
        print(f"No zip file found in {job_dir}")
        return False
    
    zip_file = zips[0]
    print(f"Unzipping {zip_file.name}...")
    with zipfile.ZipFile(zip_file, 'r') as z:
        z.extractall(job_dir)
    print("  Done.")
    return True


def parse_confidences(job_dir: Path):
    """Parse all model confidence JSONs."""
    results = []
    for i in range(5):
        conf_file = job_dir / f"fold_*_summary_confidences_{i}.json"
        matches = list(job_dir.glob(conf_file.name))
        if not matches:
            continue
        with open(matches[0]) as f:
            conf = json.load(f)
        results.append({
            "model": i,
            "iptm": conf.get("iptm"),
            "ptm": conf.get("ptm"),
            "ranking_score": conf.get("ranking_score"),
            "fraction_disordered": conf.get("fraction_disordered"),
        })
    return results


def extract_monomers(job_dir: Path):
    """Extract chain A and B from AF3 CIF to PDB."""
    from Bio.PDB import MMCIFParser, PDBIO
    
    cif_files = list(job_dir.glob("fold_*_model_0.cif"))
    if not cif_files:
        print("No model_0.cif found")
        return None, None
    
    cif_file = cif_files[0]
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", str(cif_file))
    
    io = PDBIO()
    chain_files = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            io.set_structure(chain)
            out_file = job_dir / f"ranked_0_chain_{chain_id}.pdb"
            io.save(str(out_file))
            n_res = len(list(chain.get_residues()))
            print(f"  Chain {chain_id}: {n_res} residues -> {out_file.name}")
            chain_files[chain_id] = out_file
    
    return chain_files


def fix_pdb(input_pdb: Path, output_pdb: Path):
    """Fix PDB with PDBFixer."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    
    fixer = PDBFixer(str(input_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    n_atoms = len(list(fixer.topology.atoms()))
    n_res = len(list(fixer.topology.residues()))
    print(f"  Fixed: {n_res} residues, {n_atoms} atoms -> {output_pdb.name}")


def plot_confidences(job_dir: Path, confidences):
    """Plot confidence scores for all models."""
    models = [c["model"] for c in confidences]
    iptms = [c["iptm"] for c in confidences]
    ptms = [c["ptm"] for c in confidences]
    
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(models))
    width = 0.35
    ax.bar(x - width/2, iptms, width, label="ipTM", color="coral")
    ax.bar(x + width/2, ptms, width, label="pTM", color="skyblue")
    ax.axhline(0.70, color="green", ls="--", alpha=0.5, label="Accept (0.70)")
    ax.axhline(0.60, color="orange", ls="--", alpha=0.5, label="Marginal (0.60)")
    ax.set_xlabel("AF3 Model")
    ax.set_ylabel("Score")
    ax.set_title(f"AF3 Confidence Scores ({job_dir.name})")
    ax.set_xticks(x)
    ax.set_xticklabels(models)
    ax.legend()
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(job_dir / "confidence_scores.png", dpi=200)
    plt.close(fig)
    print("  Saved confidence_scores.png")


def plot_plddt(job_dir: Path):
    """Plot pLDDT profile from full_data_0.json."""
    data_files = list(job_dir.glob("fold_*_full_data_0.json"))
    if not data_files:
        return
    
    with open(data_files[0]) as f:
        d = json.load(f)
    
    atom_plddts = np.array(d["atom_plddts"])
    chain_ids = np.array(d["atom_chain_ids"])
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    
    chains = sorted(set(chain_ids))
    colors = ["blue", "red"]
    for i, chain in enumerate(chains):
        mask = chain_ids == chain
        idx = np.where(mask)[0]
        ax = ax1 if i == 0 else ax2
        ax.plot(idx, atom_plddts[mask], ".", ms=0.5, alpha=0.5, color=colors[i], label=f"Chain {chain}")
        ax.axhline(70, color="green", ls="--", alpha=0.3)
        ax.axhline(50, color="orange", ls="--", alpha=0.3)
        ax.set_ylabel("pLDDT")
        ax.set_ylim(0, 100)
        ax.legend(loc="upper right")
    
    ax1.set_title(f"AF3 Model 0 pLDDT Profile ({job_dir.name})")
    ax2.set_xlabel("Atom index")
    fig.tight_layout()
    fig.savefig(job_dir / "plddt_profile.png", dpi=200)
    plt.close(fig)
    print("  Saved plddt_profile.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", required=True, help="AF3 job directory")
    args = parser.parse_args()
    
    job_dir = Path(args.job_dir)
    if not job_dir.exists():
        raise FileNotFoundError(f"Job directory not found: {job_dir}")
    
    print(f"\n{'='*60}")
    print(f"Processing AF3 results: {job_dir.name}")
    print(f"{'='*60}")
    
    # Step 1: Unzip
    unzip_results(job_dir)
    
    # Step 2: Parse confidences
    print("\nParsing confidence scores...")
    confidences = parse_confidences(job_dir)
    for c in confidences:
        status = "ACCEPT" if c["iptm"] >= 0.70 else ("MARGINAL" if c["iptm"] >= 0.60 else "REJECT")
        print(f"  Model {c['model']}: ipTM={c['iptm']:.2f}, pTM={c['ptm']:.2f} [{status}]")
    
    # Step 3: Plot confidences
    if confidences:
        plot_confidences(job_dir, confidences)
    
    # Step 4: Extract monomers
    print("\nExtracting monomers from Model 0...")
    chain_files = extract_monomers(job_dir)
    
    # Step 5: Fix with PDBFixer
    if chain_files:
        print("\nFixing monomers with PDBFixer...")
        for chain_id, chain_file in chain_files.items():
            name = "cgas" if chain_id == "A" else "trim41"
            fix_pdb(chain_file, job_dir / f"{name}_fixed.pdb")
    
    # Step 6: Plot pLDDT
    print("\nPlotting pLDDT profile...")
    plot_plddt(job_dir)
    
    print(f"\n✅ Processing complete for {job_dir.name}")
    print(f"   Check {job_dir} for outputs")


if __name__ == "__main__":
    main()
