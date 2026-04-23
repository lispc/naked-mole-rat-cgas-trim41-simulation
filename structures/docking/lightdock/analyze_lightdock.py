#!/usr/bin/env python3
"""Analyze LightDock poses for cGAS-TRIM41 docking."""
import sys
import glob
import numpy as np
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)

# Active residues for different species
ACTIVE_HSAP = [463, 479, 495, 498]
ACTIVE_HGAL = [463, 511, 527, 530]

def calc_ca_coord(residue):
    for atom in residue:
        if atom.id == 'CA':
            return atom.coord
    return None

def analyze_docking(dock_dir, receptor_pdb, active_residues, species_name):
    """Analyze all LightDock poses in a directory."""
    print(f"\n{'='*90}")
    print(f"LightDock Analysis: {species_name}")
    print(f"Directory: {dock_dir}")
    print(f"Active residues: {active_residues}")
    print(f"{'='*90}")
    
    # Read receptor
    trim41 = parser.get_structure("trim41", receptor_pdb)[0]['A']
    trim41_ca = {}
    for res in trim41:
        ca = calc_ca_coord(res)
        if ca is not None:
            trim41_ca[res.id[1]] = ca
    
    results = []
    pose_files = []
    for swarm_dir in sorted(glob.glob(f"{dock_dir}/swarm_*")):
        for pdb_file in sorted(glob.glob(f"{swarm_dir}/lightdock_*.pdb")):
            pose_files.append(pdb_file)
    
    print(f"Found {len(pose_files)} poses\n")
    
    for pose_file in pose_files:
        s = parser.get_structure("cgas", pose_file)[0]['A']
        cgas_ca = {}
        for res in s:
            ca = calc_ca_coord(res)
            if ca is not None:
                cgas_ca[res.id[1]] = ca
        
        # Interface
        interface_pairs = []
        for cresnum, c_ca in cgas_ca.items():
            for tresnum, t_ca in trim41_ca.items():
                dist = np.linalg.norm(c_ca - t_ca)
                if dist < 10.0:
                    interface_pairs.append((cresnum, tresnum, dist))
        
        # Active residue distances
        active_dists = {}
        for resnum in active_residues:
            if resnum in cgas_ca:
                min_dist = min(np.linalg.norm(cgas_ca[resnum] - t_ca) for t_ca in trim41_ca.values())
                active_dists[resnum] = min_dist
            else:
                active_dists[resnum] = float('inf')
        
        active_interface = [p for p in interface_pairs if p[0] in active_residues]
        max_active_dist = max(active_dists.values())
        
        results.append({
            'file': pose_file,
            'interface_pairs': len(interface_pairs),
            'active_at_interface': len(active_interface),
            'active_dists': active_dists,
            'max_active_dist': max_active_dist,
        })
    
    results.sort(key=lambda x: x['max_active_dist'])
    
    # Print top 15
    cols = " ".join([f"{r}(Å)" for r in active_residues])
    print(f"{'Pose':<40} {'Pairs':<8} {'Active@Int':<10} {cols:<40} {'MaxDist':<8} {'Status'}")
    print("-"*120)
    
    for r in results[:15]:
        status = "✅ INTERFACE" if r['max_active_dist'] < 10 else ("⚠️ CLOSE" if r['max_active_dist'] < 15 else "❌ FAR")
        dists_str = " ".join([f"{r['active_dists'][res]:<8.1f}" for res in active_residues])
        fname = r['file'].replace(f"{dock_dir}/", "")
        print(f"{fname:<40} {r['interface_pairs']:<8} {r['active_at_interface']:<10} {dists_str:<40} {r['max_active_dist']:<8.1f} {status}")
    
    good = sum(1 for r in results if r['max_active_dist'] < 10)
    print(f"\n{'='*90}")
    print(f"Total poses: {len(results)}")
    print(f"Poses with ALL active residues at interface (<10Å): {good}")
    if good > 0:
        print("✅ Found valid docking poses!")
        best = results[0]
        print(f"  Best: {best['file']} (max distance: {best['max_active_dist']:.1f}Å)")
    else:
        closest = results[0]
        print(f"⚠️ No pose has all active residues at interface")
        print(f"  Closest: {closest['file']} (max distance: {closest['max_active_dist']:.1f}Å)")
    
    return results

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: analyze_lightdock.py <dock_dir> <receptor_pdb> <species: hsap|hgal>")
        sys.exit(1)
    
    dock_dir = sys.argv[1]
    receptor_pdb = sys.argv[2]
    species = sys.argv[3].lower()
    active = ACTIVE_HSAP if species == "hsap" else ACTIVE_HGAL
    species_name = "Human cGAS + TRIM41" if species == "hsap" else "NMR cGAS + TRIM41"
    
    analyze_docking(dock_dir, receptor_pdb, active, species_name)
