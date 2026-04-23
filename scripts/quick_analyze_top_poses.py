#!/usr/bin/env python3
"""Quick analysis of LightDock top poses for Hsap restrained docking."""

import sys
import glob
import numpy as np
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)

# Hsap active residues
ACTIVE_HSAP = [463, 479, 495, 498]

def calc_ca_coord(residue):
    for atom in residue:
        if atom.id == 'CA':
            return atom.coord
    return None

dock_dir = sys.argv[1] if len(sys.argv) > 1 else "structures/docking/lightdock/Hsap_restrained"
receptor_pdb = sys.argv[2] if len(sys.argv) > 2 else "structures/af3_raw/job1_Hsap_WT/trim41_SPRY_413-630.pdb"

# Read receptor
trim41 = parser.get_structure("trim41", receptor_pdb)[0]['A']
trim41_ca = {}
for res in trim41:
    ca = calc_ca_coord(res)
    if ca is not None:
        trim41_ca[res.id[1]] = ca

# Analyze top poses
top_files = sorted(glob.glob(f"{dock_dir}/top_*.pdb"))
print(f"\n{'='*100}")
print(f"  Hsap Restrained Docking Analysis ({len(top_files)} top poses)")
print(f"  Active residues: {ACTIVE_HSAP}")
print(f"{'='*100}")

cols = " ".join([f"{r}(Å)" for r in ACTIVE_HSAP])
print(f"\n{'Pose':<12} {'Score':<10} {'Pairs':<8} {'Active@Int':<10} {cols:<45} {'MaxDist':<8} {'Status'}")
print("-"*130)

results = []
for pose_file in top_files:
    s = parser.get_structure("cgas", pose_file)[0]['A']
    cgas_ca = {}
    for res in s:
        ca = calc_ca_coord(res)
        if ca is not None:
            cgas_ca[res.id[1]] = ca
    
    # Interface pairs (<10A between any cGAS and TRIM41 CA)
    interface_pairs = []
    for cresnum, c_ca in cgas_ca.items():
        for tresnum, t_ca in trim41_ca.items():
            dist = np.linalg.norm(c_ca - t_ca)
            if dist < 10.0:
                interface_pairs.append((cresnum, tresnum, dist))
    
    # Active residue distances to nearest TRIM41 CA
    active_dists = {}
    for resnum in ACTIVE_HSAP:
        if resnum in cgas_ca:
            min_dist = min(np.linalg.norm(cgas_ca[resnum] - t_ca) for t_ca in trim41_ca.values())
            active_dists[resnum] = min_dist
        else:
            active_dists[resnum] = float('inf')
    
    active_interface = [p for p in interface_pairs if p[0] in ACTIVE_HSAP]
    max_active_dist = max(active_dists.values())
    
    fname = pose_file.split('/')[-1]
    dists_str = " ".join([f"{active_dists[res]:<8.1f}" for res in ACTIVE_HSAP])
    status = "✅ ALL<10" if max_active_dist < 10 else ("⚠️ CLOSE" if max_active_dist < 15 else "❌ FAR")
    
    print(f"{fname:<12} {'N/A':<10} {len(interface_pairs):<8} {len(active_interface):<10} {dists_str:<45} {max_active_dist:<8.1f} {status}")
    
    results.append({
        'file': fname,
        'interface_pairs': len(interface_pairs),
        'active_at_interface': len(active_interface),
        'active_dists': active_dists,
        'max_active_dist': max_active_dist,
    })

# Summary
good = sum(1 for r in results if r['max_active_dist'] < 10)
partial = sum(1 for r in results if 10 <= r['max_active_dist'] < 15)
print(f"\n{'='*100}")
print(f"Summary:")
print(f"  Total poses analyzed: {len(results)}")
print(f"  ALL 4 active residues at interface (<10Å): {good}")
print(f"  Partial (at least 1 >10Å but max <15Å): {partial}")

if good > 0:
    best = min(results, key=lambda x: x['max_active_dist'])
    print(f"\n  ✅ SUCCESS! Best pose: {best['file']}")
    print(f"     Max active residue distance: {best['max_active_dist']:.1f}Å")
    for res in ACTIVE_HSAP:
        print(f"     Residue {res}: {best['active_dists'][res]:.1f}Å")
else:
    closest = min(results, key=lambda x: x['max_active_dist'])
    print(f"\n  ⚠️ No pose has all 4 active residues at interface")
    print(f"  Closest pose: {closest['file']}")
    print(f"  Max active residue distance: {closest['max_active_dist']:.1f}Å")
    for res in ACTIVE_HSAP:
        d = closest['active_dists'][res]
        status = "✅" if d < 10 else "❌"
        print(f"     {status} Residue {res}: {d:.1f}Å")
