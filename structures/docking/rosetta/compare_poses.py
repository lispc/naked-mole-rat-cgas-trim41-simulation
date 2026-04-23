#!/usr/bin/env python3
"""Compare Rosetta docking poses with LightDock reference."""

import sys
import numpy as np

def read_pdb_atoms(filepath):
    """Read CA atoms from PDB, return dict: (chain, resi) -> (x, y, z)"""
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
    """Calculate RMSD between two sets of atoms."""
    common = set(atoms1.keys()) & set(atoms2.keys())
    if len(common) == 0:
        return None
    coords1 = np.array([atoms1[k] for k in sorted(common)])
    coords2 = np.array([atoms2[k] for k in sorted(common)])
    
    # Center
    c1 = coords1.mean(axis=0)
    c2 = coords2.mean(axis=0)
    coords1 -= c1
    coords2 -= c2
    
    # Kabsch algorithm
    H = coords1.T @ coords2
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    coords2_rot = coords2 @ R
    rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2_rot)**2, axis=1)))
    return rmsd

def main():
    ref_file = sys.argv[1] if len(sys.argv) > 1 else "input.pdb"
    
    ref_atoms = read_pdb_atoms(ref_file)
    print(f"Reference: {ref_file} ({len(ref_atoms)} CA atoms)")
    
    import glob
    poses = sorted(glob.glob("output_global/*.pdb"))
    
    print(f"\n{'Pose':<20} {'RMSD(Å)':>10} {'CA matched':>10}")
    print("-" * 45)
    
    for pose_file in poses:
        pose_atoms = read_pdb_atoms(pose_file)
        rmsd = calc_rmsd(ref_atoms, pose_atoms)
        common = set(ref_atoms.keys()) & set(pose_atoms.keys())
        print(f"{pose_file:<20} {rmsd:>10.3f} {len(common):>10}")

if __name__ == "__main__":
    main()
