#!/usr/bin/env python3
"""
Rosetta FastRelax mutational scanning.

Uses FastRelax (backbone + sidechain) instead of Pack+Min to capture
allosteric effects on interface energy.

For Hgal: S463D + E511K + Y527L + T530K (4mut_rev)
For Hsap: D431S + K479E + L495Y + K498T (4mut)

Usage:
    conda activate rosetta
    python scripts/rosetta_fastrelax_mutscan.py
"""

import os
import sys
import json

os.environ["PYROSETTA_SILENT"] = "1"
import pyrosetta
pyrosetta.init("-mute all -ignore_unrecognized_res")

from pyrosetta import get_fa_scorefxn, pose_from_pdb
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.relax import FastRelax


def compute_interface_energy(pose, sfxn, chain_A=1, chain_B=2):
    """Compute interface score: E_total - E_A - E_B."""
    pose.energies().clear_energies()
    pose.update_residue_neighbors()
    E_total = sfxn(pose)
    
    pose_A = pose.split_by_chain(chain_A)
    pose_A.energies().clear_energies()
    pose_A.update_residue_neighbors()
    E_A = sfxn(pose_A)
    
    pose_B = pose.split_by_chain(chain_B)
    pose_B.energies().clear_energies()
    pose_B.update_residue_neighbors()
    E_B = sfxn(pose_B)
    
    I_sc = E_total - (E_A + E_B)
    return {
        "total_score": E_total,
        "E_A": E_A,
        "E_B": E_B,
        "I_sc": I_sc
    }


def apply_mutations(pose, mutations, chain="B"):
    """Apply mutations to pose."""
    aa_map = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    
    for pdb_resi, new_aa in mutations:
        pose_resi = pose.pdb_info().pdb2pose(chain[0], pdb_resi)
        if pose_resi == 0:
            raise ValueError(f"Residue {chain}{pdb_resi} not found in pose")
        
        new_aa3 = aa_map[new_aa.upper()]
        mutator = MutateResidue(pose_resi, new_aa3)
        mutator.apply(pose)
        print(f"  Mutated {chain}{pdb_resi} -> {new_aa3}")


def run_fastrelax_scan(name, pdb_file, mutations, output_dir="structures/docking/rosetta"):
    """Run FastRelax scan for one system (WT + combined 4mut only)."""
    print(f"\n{'='*70}")
    print(f"FastRelax Scan: {name}")
    print(f"PDB: {pdb_file}")
    print(f"Mutations: {' + '.join(f'{r}{a}' for r, a in mutations)}")
    print(f"{'='*70}")
    
    # Load pose
    pose = pose_from_pdb(pdb_file)
    print(f"Loaded pose: {pose.total_residue()} residues")
    
    # Score function
    sfxn = get_fa_scorefxn()
    
    # Setup FastRelax
    relax = FastRelax()
    relax.set_scorefxn(sfxn)
    relax.max_iter(200)  # Limit for speed
    
    results = {
        "system": name,
        "pdb": pdb_file,
        "mutations": [f"{r}{a}" for r, a in mutations],
        "wt": None,
        "combined": None
    }
    
    # WT
    print("\n--- WT ---")
    wt_pose = pose.clone()
    print("Running FastRelax on WT...")
    relax.apply(wt_pose)
    wt_energy = compute_interface_energy(wt_pose, sfxn)
    for k, v in wt_energy.items():
        print(f"  {k}: {v:.3f}")
    results["wt"] = wt_energy
    
    # Save relaxed WT
    wt_out = os.path.join(output_dir, f"{name.lower()}_wt_relaxed.pdb")
    wt_pose.dump_pdb(wt_out)
    print(f"  Saved: {wt_out}")
    
    # Combined 4mut
    print(f"\n--- Combined 4mut ---")
    mut_pose = pose.clone()
    apply_mutations(mut_pose, mutations)
    print("Running FastRelax on 4mut...")
    relax.apply(mut_pose)
    mut_energy = compute_interface_energy(mut_pose, sfxn)
    for k, v in mut_energy.items():
        print(f"  {k}: {v:.3f}")
    results["combined"] = mut_energy
    
    # Save relaxed mutant
    mut_out = os.path.join(output_dir, f"{name.lower()}_4mut_relaxed.pdb")
    mut_pose.dump_pdb(mut_out)
    print(f"  Saved: {mut_out}")
    
    # Comparison
    ddG = mut_energy["I_sc"] - wt_energy["I_sc"]
    ddG_total = mut_energy["total_score"] - wt_energy["total_score"]
    ddG_B = mut_energy["E_B"] - wt_energy["E_B"]
    
    print(f"\n--- Comparison ---")
    print(f"  ΔΔG (I_sc):    {ddG:+.3f} REU")
    print(f"  ΔΔG (total):   {ddG_total:+.3f} REU")
    print(f"  ΔΔG (E_B):     {ddG_B:+.3f} REU")
    
    results["ddG_I_sc"] = ddG
    results["ddG_total"] = ddG_total
    results["ddG_E_B"] = ddG_B
    
    # Save results
    out_json = os.path.join(output_dir, f"fastrelax_{name.lower()}.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {out_json}")
    
    return results


def main():
    base = "/Users/zhangzhuo/repos/personal/naked-mole-rat-cgas-trim41-simulation"
    
    # Hgal
    hgal_results = run_fastrelax_scan(
        "Hgal",
        os.path.join(base, "structures/docking/rosetta/input.pdb"),
        [(463, 'D'), (511, 'K'), (527, 'L'), (530, 'K')],
        os.path.join(base, "structures/docking/rosetta")
    )
    
    # Hsap
    hsap_results = run_fastrelax_scan(
        "Hsap",
        os.path.join(base, "structures/docking/rosetta/hsap_input.pdb"),
        [(431, 'S'), (479, 'E'), (495, 'Y'), (498, 'T')],
        os.path.join(base, "structures/docking/rosetta")
    )
    
    # Summary
    print(f"\n{'='*70}")
    print("FASTRELAX SUMMARY")
    print(f"{'='*70}")
    print(f"\nHgal WT -> 4mut_rev:")
    print(f"  WT I_sc:     {hgal_results['wt']['I_sc']:.3f} REU")
    print(f"  4mut I_sc:   {hgal_results['combined']['I_sc']:.3f} REU")
    print(f"  ΔΔG (I_sc):  {hgal_results['ddG_I_sc']:+.3f} REU")
    print(f"  ΔΔG (total): {hgal_results['ddG_total']:+.3f} REU")
    print(f"\nHsap WT -> 4mut:")
    print(f"  WT I_sc:     {hsap_results['wt']['I_sc']:.3f} REU")
    print(f"  4mut I_sc:   {hsap_results['combined']['I_sc']:.3f} REU")
    print(f"  ΔΔG (I_sc):  {hsap_results['ddG_I_sc']:+.3f} REU")
    print(f"  ΔΔG (total): {hsap_results['ddG_total']:+.3f} REU")
    
    print(f"\nComparison with Pack+Min (previous):")
    print(f"  Hgal ΔΔG: Pack+Min = 0.000, FastRelax = {hgal_results['ddG_I_sc']:+.3f}")
    print(f"  Hsap ΔΔG: Pack+Min = 0.000, FastRelax = {hsap_results['ddG_I_sc']:+.3f}")


if __name__ == "__main__":
    main()
