#!/usr/bin/env python3
"""
Rosetta mutational scanning: compute interface energy change upon 4mut.

Uses direct score function with energy cache clearing for reliable results.

For Hgal: S463D + E511K + Y527L + T530K (4mut_rev)
For Hsap: D431S + K479E + L495Y + K498T (4mut)

Usage:
    conda activate rosetta
    python scripts/rosetta_mutational_scan.py
"""

import os
import sys
import json

# Initialize PyRosetta
os.environ["PYROSETTA_SILENT"] = "1"
import pyrosetta
pyrosetta.init("-mute all -ignore_unrecognized_res")

from pyrosetta import get_fa_scorefxn, pose_from_pdb
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking, IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT
)
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector, NotResidueSelector, NeighborhoodResidueSelector
)
from pyrosetta.rosetta.core.kinematics import MoveMap


def compute_interface_energy(pose, sfxn, chain_A=1, chain_B=2):
    """Compute interface score: E_total - E_A - E_B.
    
    Clears energy caches to ensure fresh calculation.
    """
    # Force re-evaluation by clearing energy cache
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


def repack_pose(pose, sfxn, focus_residues=None, neighbor_distance=10.0):
    """Repack side chains, optionally focusing on specific residues and neighbors."""
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(IncludeCurrent())
    
    if focus_residues:
        focus_selector = ResidueIndexSelector(",".join(str(r) for r in focus_residues))
        neighbor_selector = NeighborhoodResidueSelector(focus_selector, neighbor_distance, False)
        not_neighbor = NotResidueSelector(neighbor_selector)
        prevent = PreventRepackingRLT()
        tf.push_back(OperateOnResidueSubset(prevent, not_neighbor))
    
    task = tf.create_task_and_apply_taskoperations(pose)
    packer = PackRotamersMover(sfxn, task)
    packer.apply(pose)


def minimize_pose(pose, sfxn, focus_residues=None, neighbor_distance=10.0):
    """Minimize side chains and backbone, optionally focusing on specific residues."""
    mm = MoveMap()
    mm.set_bb(False)  # Keep backbone fixed for speed
    
    if focus_residues:
        focus_selector = ResidueIndexSelector(",".join(str(r) for r in focus_residues))
        neighbor_selector = NeighborhoodResidueSelector(focus_selector, neighbor_distance, False)
        subset = neighbor_selector.apply(pose)
        
        for i in range(1, pose.total_residue() + 1):
            mm.set_chi(i, subset[i])
    else:
        mm.set_chi(True)
    
    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sfxn)
    min_mover.min_type("lbfgs_armijo_nonmonotone")
    min_mover.tolerance(0.01)
    min_mover.max_iter(200)
    min_mover.apply(pose)


def mutate_and_score(pose, sfxn, mutations, chain="B"):
    """
    Apply mutations and compute interface energy.
    mutations: list of (pdb_resi, new_aa_1letter)
    Returns: (mutant_pose_clone, energy_dict)
    """
    mutant = pose.clone()
    
    pose_residues = []
    for pdb_resi, new_aa in mutations:
        pose_resi = mutant.pdb_info().pdb2pose(chain[0], pdb_resi)
        if pose_resi == 0:
            raise ValueError(f"Residue {chain}{pdb_resi} not found in pose")
        
        aa_map = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
            'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
            'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
            'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
        }
        new_aa3 = aa_map[new_aa.upper()]
        
        mutator = MutateResidue(pose_resi, new_aa3)
        mutator.apply(mutant)
        pose_residues.append(pose_resi)
        print(f"  Mutated {chain}{pdb_resi} (pose {pose_resi}) -> {new_aa3}")
    
    # Debug: score before optimization
    pre = compute_interface_energy(mutant, sfxn)
    print(f"  Pre-opt:  total={pre['total_score']:.1f}  I_sc={pre['I_sc']:.3f}")
    
    # Repack mutated residues and neighbors
    print(f"  Repacking...")
    repack_pose(mutant, sfxn, focus_residues=pose_residues, neighbor_distance=10.0)
    post_pack = compute_interface_energy(mutant, sfxn)
    print(f"  Post-pack: total={post_pack['total_score']:.1f}  I_sc={post_pack['I_sc']:.3f}")
    
    # Minimize
    print(f"  Minimizing...")
    minimize_pose(mutant, sfxn, focus_residues=pose_residues, neighbor_distance=10.0)
    post_min = compute_interface_energy(mutant, sfxn)
    print(f"  Post-min:  total={post_min['total_score']:.1f}  I_sc={post_min['I_sc']:.3f}")
    
    return mutant, post_min


def run_scan(name, pdb_file, mutations, output_dir="structures/docking/rosetta"):
    """Run mutational scan for one system."""
    print(f"\n{'='*70}")
    print(f"System: {name}")
    print(f"PDB: {pdb_file}")
    print(f"Mutations: {' + '.join(f'{r}{a}' for r, a in mutations)}")
    print(f"{'='*70}")
    
    # Load pose
    pose = pose_from_pdb(pdb_file)
    print(f"Loaded pose: {pose.total_residue()} residues")
    
    # Score function
    sfxn = get_fa_scorefxn()
    
    # WT energy
    print("\n--- WT Energy ---")
    wt_energy = compute_interface_energy(pose, sfxn)
    for k, v in wt_energy.items():
        print(f"  {k}: {v:.3f}")
    
    results = {
        "system": name,
        "pdb": pdb_file,
        "wt": wt_energy,
        "single_mutations": {},
        "combined": None
    }
    
    # Single mutations
    print("\n--- Single Mutations ---")
    for pdb_resi, new_aa in mutations:
        mut_name = f"{pdb_resi}{new_aa}"
        print(f"\n> {mut_name}")
        _, mut_energy = mutate_and_score(pose, sfxn, [(pdb_resi, new_aa)])
        ddG = mut_energy["I_sc"] - wt_energy["I_sc"]
        print(f"  Final ΔΔG = {ddG:+.3f} REU")
        results["single_mutations"][mut_name] = {
            "energy": mut_energy,
            "ddG_I_sc": ddG
        }
    
    # Combined 4mut
    print(f"\n--- Combined 4mut ---")
    _, comb_energy = mutate_and_score(pose, sfxn, mutations)
    ddG_comb = comb_energy["I_sc"] - wt_energy["I_sc"]
    print(f"  Final ΔΔG = {ddG_comb:+.3f} REU")
    results["combined"] = {
        "energy": comb_energy,
        "ddG_I_sc": ddG_comb
    }
    
    # Save results
    out_json = os.path.join(output_dir, f"mutscan_{name.lower()}.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {out_json}")
    
    return results


def main():
    base = "/Users/zhangzhuo/repos/personal/naked-mole-rat-cgas-trim41-simulation"
    
    # Hgal: 4mut_rev = S463D + E511K + Y527L + T530K
    hgal_mutations = [(463, 'D'), (511, 'K'), (527, 'L'), (530, 'K')]
    
    hgal_results = run_scan(
        "Hgal_WT_to_4mut_rev",
        os.path.join(base, "structures/docking/rosetta/input.pdb"),
        hgal_mutations,
        os.path.join(base, "structures/docking/rosetta")
    )
    
    # Hsap: 4mut = D431S + K479E + L495Y + K498T
    hsap_mutations = [(431, 'S'), (479, 'E'), (495, 'Y'), (498, 'T')]
    
    hsap_results = run_scan(
        "Hsap_WT_to_4mut",
        os.path.join(base, "structures/docking/rosetta/hsap_input.pdb"),
        hsap_mutations,
        os.path.join(base, "structures/docking/rosetta")
    )
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"\nHgal WT -> 4mut_rev:")
    print(f"  WT I_sc:     {hgal_results['wt']['I_sc']:.3f} REU")
    print(f"  4mut I_sc:   {hgal_results['combined']['energy']['I_sc']:.3f} REU")
    print(f"  ΔΔG:         {hgal_results['combined']['ddG_I_sc']:+.3f} REU")
    print(f"\nHsap WT -> 4mut:")
    print(f"  WT I_sc:     {hsap_results['wt']['I_sc']:.3f} REU")
    print(f"  4mut I_sc:   {hsap_results['combined']['energy']['I_sc']:.3f} REU")
    print(f"  ΔΔG:         {hsap_results['combined']['ddG_I_sc']:+.3f} REU")


if __name__ == "__main__":
    main()
