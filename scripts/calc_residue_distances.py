#!/usr/bin/env python3
"""Calculate exact CA-CA distances between active residues in cGAS CTD structures."""

from Bio.PDB import PDBParser
import numpy as np
import os

parser = PDBParser(QUIET=True)

def get_ca_coord(structure, resnum, chain_id='A'):
    """Get CA atom coordinates for a given residue number."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id or chain_id is None:
                for residue in chain:
                    if residue.id[1] == resnum:
                        if 'CA' in residue:
                            return residue['CA'].get_coord()
    return None

def calc_dist(structure, res1, res2, name=""):
    """Calculate distance between two CA atoms."""
    c1 = get_ca_coord(structure, res1)
    c2 = get_ca_coord(structure, res2)
    if c1 is None or c2 is None:
        print(f"  {name}: MISSING atom!")
        return None
    dist = np.linalg.norm(c1 - c2)
    print(f"  {name}: {dist:.2f} Å")
    return dist

def analyze_structure(pdb_path, name, residues):
    """Analyze all pairwise distances between active residues."""
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    
    structure = parser.get_structure(name, pdb_path)
    
    res_nums = list(residues.keys())
    res_names = list(residues.values())
    
    print(f"\n  Residues: {', '.join([f'{n}({r})' for n, r in zip(res_names, res_nums)])}")
    print(f"\n  Pairwise CA-CA distances:")
    
    distances = {}
    for i in range(len(res_nums)):
        for j in range(i+1, len(res_nums)):
            label = f"{res_names[i]}{res_nums[i]}-{res_names[j]}{res_nums[j]}"
            d = calc_dist(structure, res_nums[i], res_nums[j], label)
            distances[label] = d
    
    # Find max distance
    max_pair = max(distances, key=distances.get)
    max_dist = distances[max_pair]
    print(f"\n  Maximum distance: {max_pair} = {max_dist:.2f} Å")
    
    # Find min distance
    min_pair = min(distances, key=distances.get)
    min_dist = distances[min_pair]
    print(f"  Minimum distance: {min_pair} = {min_dist:.2f} Å")
    
    return distances, max_dist, min_dist

# Hgal cGAS CTD (WT)
hgal_residues = {463: 'C', 511: 'K', 527: 'L', 530: 'K'}
hgal_dists, hgal_max, hgal_min = analyze_structure(
    'structures/af3_raw/job3_Hgal_WT/cgas_CT_200-554.pdb',
    'H. glaber cGAS CTD (WT)',
    hgal_residues
)

# Hsap cGAS CTD (WT) — CORRECTED mutations: D431S, K479E, L495Y, K498T
# NOTE: Old incorrect mapping was C463S/K479E/L495Y/K498T (deleted)
hsap_residues = {431: 'D', 479: 'K', 495: 'L', 498: 'K'}
hsap_dists, hsap_max, hsap_min = analyze_structure(
    'structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb',
    'H. sapiens cGAS CTD (WT) — correct sites: D431, K479, L495, K498',
    hsap_residues
)

# Print comparison summary
print(f"\n{'='*60}")
print(f"  COMPARISON SUMMARY")
print(f"{'='*60}")
print(f"  Hgal max distance: {hgal_max:.2f} Å")
print(f"  Hsap max distance: {hsap_max:.2f} Å")
print(f"  Difference: {hsap_max - hgal_max:.2f} Å")
print(f"\n  Hgal compact patch: 4 residues within ~{hgal_max:.1f} Å")
print(f"  Hsap dispersed sites: 4 residues span ~{hsap_max:.1f} Å")

# Also save to file for figure generation
os.makedirs('figures', exist_ok=True)
with open('figures/distance_data.txt', 'w') as f:
    f.write("# Active residue CA-CA distances\n\n")
    f.write("Hgal (H. glaber cGAS CTD):\n")
    for k, v in hgal_dists.items():
        f.write(f"  {k}: {v:.2f}\n")
    f.write(f"  max: {hgal_max:.2f}\n\n")
    f.write("Hsap (H. sapiens cGAS CTD):\n")
    for k, v in hsap_dists.items():
        f.write(f"  {k}: {v:.2f}\n")
    f.write(f"  max: {hsap_max:.2f}\n")

print(f"\n  Data saved to: figures/distance_data.txt")
