#!/usr/bin/env python3
"""Analyze Rosetta docking results and compare with LightDock reference."""

import numpy as np
import glob
import re

def read_pdb_atoms(filepath):
    """Read CA atoms from PDB."""
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
    """Kabsch RMSD."""
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

def read_scorefile(filepath):
    """Read Rosetta scorefile."""
    scores = {}
    with open(filepath) as f:
        header = None
        for line in f:
            if line.startswith("SCORE: total_score"):
                header = line.split()
            elif line.startswith("SCORE:") and header:
                parts = line.split()
                name = parts[-1]
                scores[name] = {}
                for i, key in enumerate(header[1:-1], start=1):
                    try:
                        scores[name][key] = float(parts[i])
                    except:
                        scores[name][key] = parts[i]
    return scores

def calc_interface_contacts(pdb_file, chainA="A", chainB="B", cutoff=5.0):
    """Count inter-chain contacts within cutoff (Å)."""
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
        return 0, 0
    
    coords_A = np.array(coords_A)
    coords_B = np.array(coords_B)
    
    # Count contacts
    contacts = 0
    for a in coords_A:
        dists = np.linalg.norm(coords_B - a, axis=1)
        contacts += np.sum(dists < cutoff)
    
    # Buried surface area approximation
    return contacts, len(coords_A), len(coords_B)

def main():
    ref_file = "input.pdb"
    ref_atoms = read_pdb_atoms(ref_file)
    print(f"Reference (LightDock): {ref_file} ({len(ref_atoms)} CA atoms)")
    
    scores = read_scorefile("output_global/global.sc")
    
    poses = sorted(glob.glob("output_global/*.pdb"))
    
    print(f"\n{'='*100}")
    print(f"{'Pose':<18} {'I_sc':>10} {'Total':>10} {'RMSD(Å)':>10} {'Irms':>10} {'Fnat':>10} {'Contacts':>10}")
    print(f"{'='*100}")
    
    results = []
    for pose_file in poses:
        basename = pose_file.split("/")[-1].replace(".pdb", "")
        pose_atoms = read_pdb_atoms(pose_file)
        rmsd = calc_rmsd(ref_atoms, pose_atoms)
        
        sc = scores.get(basename, {})
        isc = sc.get('I_sc', 0)
        total = sc.get('total_score', 0)
        irms = sc.get('Irms', 0)
        fnat = sc.get('Fnat', 0)
        
        contacts, nA, nB = calc_interface_contacts(pose_file)
        
        results.append({
            'name': basename,
            'I_sc': isc,
            'total': total,
            'rmsd': rmsd,
            'Irms': irms,
            'Fnat': fnat,
            'contacts': contacts
        })
        
        print(f"{basename:<18} {isc:>10.2f} {total:>10.2f} {rmsd:>10.3f} {irms:>10.3f} {fnat:>10.3f} {contacts:>10}")
    
    # Sort by I_sc
    results_sorted = sorted(results, key=lambda x: x['I_sc'])
    best = results_sorted[0]
    
    print(f"\n{'='*100}")
    print(f"BEST DECOY (by I_sc): {best['name']}")
    print(f"  Interface score:    {best['I_sc']:.2f} REU")
    print(f"  Total score:        {best['total']:.2f} REU")
    print(f"  CA RMSD vs LightDock: {best['rmsd']:.3f} Å")
    print(f"  Interface RMSD:     {best['Irms']:.3f} Å")
    print(f"  Native contacts:    {best['Fnat']:.3f}")
    print(f"  Inter-chain contacts: {best['contacts']}")
    print(f"{'='*100}")
    
    # Statistics
    rmsds = [r['rmsd'] for r in results]
    iscs = [r['I_sc'] for r in results]
    print(f"\nSTATISTICS (n={len(results)}):")
    print(f"  RMSD range:     {min(rmsds):.2f} - {max(rmsds):.2f} Å  (mean: {np.mean(rmsds):.2f})")
    print(f"  I_sc range:     {min(iscs):.2f} - {max(iscs):.2f} REU  (mean: {np.mean(iscs):.2f})")
    
    # Check consistency
    close_poses = [r for r in results if r['rmsd'] < 5.0]
    print(f"\nCONSISTENCY WITH LIGHTDOCK:")
    print(f"  Decoys within 5Å RMSD: {len(close_poses)}/{len(results)}")
    if close_poses:
        best_close = min(close_poses, key=lambda x: x['I_sc'])
        print(f"  Best I_sc among close decoys: {best_close['I_sc']:.2f} REU ({best_close['name']})")

if __name__ == "__main__":
    main()
