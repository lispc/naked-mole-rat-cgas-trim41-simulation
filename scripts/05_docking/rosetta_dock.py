#!/usr/bin/env python3
"""
Rosetta docking for 4mut structures using PyRosetta DockingProtocol.

Usage:
    conda activate rosetta
    python scripts/rosetta_dock_4mut.py
"""

import os
import sys
import json
import numpy as np

os.environ["PYROSETTA_SILENT"] = "1"
import pyrosetta
pyrosetta.init("-mute all -ignore_unrecognized_res")

from pyrosetta import get_fa_scorefxn, pose_from_pdb
from pyrosetta.rosetta.protocols.docking import DockingProtocol


def read_pdb_atoms(filepath):
    atoms = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM  "):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    chain = line[21].strip()
                    resi = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms[(chain, resi)] = np.array([x, y, z])
    return atoms


def calc_rmsd(atoms1, atoms2):
    common = set(atoms1.keys()) & set(atoms2.keys())
    if len(common) == 0:
        return None
    coords1 = np.array([atoms1[k] for k in sorted(common)])
    coords2 = np.array([atoms2[k] for k in sorted(common)])
    c1 = coords1.mean(axis=0)
    c2 = coords2.mean(axis=0)
    coords1 -= c1
    coords2 -= c2
    H = coords1.T @ coords2
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    coords2_rot = coords2 @ R
    rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2_rot)**2, axis=1)))
    return rmsd


def calc_interface_contacts(pdb_file, chainA="A", chainB="B", cutoff=5.0):
    coords_A = []
    coords_B = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if chain == chainA:
                    coords_A.append(np.array([x, y, z]))
                elif chain == chainB:
                    coords_B.append(np.array([x, y, z]))
    
    if not coords_A or not coords_B:
        return 0
    
    coords_A = np.array(coords_A)
    coords_B = np.array(coords_B)
    
    contacts = 0
    for a in coords_A:
        dists = np.linalg.norm(coords_B - a, axis=1)
        contacts += np.sum(dists < cutoff)
    
    return contacts


def get_residue_distance(pdb_file, chainA, resA, chainB, resB):
    """Get minimum distance between two residues."""
    coords_A = []
    coords_B = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21].strip()
                resi = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if chain == chainA and resi == resA:
                    coords_A.append(np.array([x, y, z]))
                elif chain == chainB and resi == resB:
                    coords_B.append(np.array([x, y, z]))
    
    if not coords_A or not coords_B:
        return None
    
    coords_A = np.array(coords_A)
    coords_B = np.array(coords_B)
    
    min_dist = float('inf')
    for a in coords_A:
        dists = np.linalg.norm(coords_B - a, axis=1)
        min_dist = min(min_dist, dists.min())
    
    return min_dist


def run_docking(name, input_pdb, output_dir, nstruct=10, partners="A_B"):
    """Run Rosetta docking protocol."""
    print(f"\n{'='*70}")
    print(f"Docking: {name}")
    print(f"Input: {input_pdb}")
    print(f"Decoys: {nstruct}")
    print(f"{'='*70}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load starting pose
    start_pose = pose_from_pdb(input_pdb)
    print(f"Loaded pose: {start_pose.total_residue()} residues")
    
    # Setup docking protocol
    sfxn = get_fa_scorefxn()
    dock = DockingProtocol()
    dock.set_partners(partners)
    dock.set_dock_min(False)
    dock.set_low_res_protocol_only(False)
    
    # Reference for RMSD
    ref_atoms = read_pdb_atoms(input_pdb)
    
    results = []
    best_pose = None
    best_score = float('inf')
    
    for i in range(1, nstruct + 1):
        print(f"\n  Decoy {i}/{nstruct}...", end=" ")
        sys.stdout.flush()
        
        pose = start_pose.clone()
        dock.apply(pose)
        
        # Convert back to full-atom if needed (docking uses centroid mode)
        if not pose.is_fullatom():
            from pyrosetta.rosetta.core.util import switch_to_residue_type_set
            switch_to_residue_type_set(pose, 'fa_standard')
        
        score = sfxn(pose)
        
        # Save decoy
        decoy_file = os.path.join(output_dir, f"{name}_{i:04d}.pdb")
        pose.dump_pdb(decoy_file)
        
        # Calculate metrics
        decoy_atoms = read_pdb_atoms(decoy_file)
        rmsd = calc_rmsd(ref_atoms, decoy_atoms)
        contacts = calc_interface_contacts(decoy_file)
        
        # 4mut distances
        if "hgal" in name.lower():
            mut_positions = [463, 511, 527, 530]
        else:
            mut_positions = [431, 479, 495, 498]
        
        mut_dists = {}
        for p in mut_positions:
            d = get_residue_distance(decoy_file, "B", p, "A", 500)  # Approximate interface center
            if d:
                mut_dists[p] = d
        
        result = {
            "name": f"{name}_{i:04d}",
            "score": float(score),
            "rmsd": float(rmsd) if rmsd is not None else None,
            "contacts": int(contacts),
            "mut_dists": {int(k): float(v) for k, v in mut_dists.items()}
        }
        results.append(result)
        
        if score < best_score:
            best_score = score
            best_pose = decoy_file
        
        print(f"score={score:.1f} RMSD={rmsd:.2f}Å contacts={contacts}")
    
    # Save results
    out_json = os.path.join(output_dir, f"{name}_results.json")
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Summary
    print(f"\n{'='*70}")
    print(f"Results summary ({nstruct} decoys):")
    print(f"  Best score: {best_score:.1f} ({os.path.basename(best_pose)})")
    print(f"  Score range: {min(r['score'] for r in results):.1f} - {max(r['score'] for r in results):.1f}")
    print(f"  RMSD range: {min(r['rmsd'] for r in results):.2f} - {max(r['rmsd'] for r in results):.2f} Å")
    print(f"  Contacts range: {min(r['contacts'] for r in results)} - {max(r['contacts'] for r in results)}")
    
    # Best decoy 4mut distances
    best_result = min(results, key=lambda x: x['score'])
    print(f"\n  Best decoy 4mut distances to interface:")
    for p, d in sorted(best_result['mut_dists'].items()):
        print(f"    {p}: {d:.1f} Å")
    
    print(f"\n  Output: {output_dir}")
    print(f"  JSON: {out_json}")
    
    return results


def main():
    base = "/Users/zhangzhuo/repos/personal/naked-mole-rat-cgas-trim41-simulation"
    
    # Hgal 4mut_rev docking
    hgal_results = run_docking(
        "hgal_4mut_rev",
        os.path.join(base, "structures/docking/rosetta/input_4mut_rev.pdb"),
        os.path.join(base, "structures/docking/rosetta/output_4mut_rev"),
        nstruct=10,
        partners="A_B"
    )
    
    # Hsap 4mut docking
    hsap_results = run_docking(
        "hsap_4mut",
        os.path.join(base, "structures/docking/rosetta/hsap_input_4mut.pdb"),
        os.path.join(base, "structures/docking/rosetta/output_hsap_4mut"),
        nstruct=10,
        partners="A_B"
    )


if __name__ == "__main__":
    main()
