#!/usr/bin/env python3
"""
Prepare Rosetta docking input for 4mut structures.

Steps:
1. Extract CTD region (200-554 for Hgal, 200-522 for Hsap) from 4mut monomer PDB
2. Align to WT cgas_CT coordinates using Kabsch algorithm
3. Merge with TRIM41 SPRY to create docking input
"""

import numpy as np
import argparse


def read_pdb_atoms(filepath):
    """Read PDB ATOM/HETATM records."""
    atoms = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(line)
    return atoms


def parse_ca_coords(atoms):
    """Parse CA atom coordinates from PDB lines."""
    coords = {}
    for line in atoms:
        if line[12:16].strip() == "CA":
            resi = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords[resi] = np.array([x, y, z])
    return coords


def kabsch_align(mobile_coords, target_coords, common_residues):
    """Align mobile_coords to target_coords using Kabsch algorithm.
    Returns rotation matrix R and translation vector t.
    """
    m = np.array([mobile_coords[r] for r in common_residues])
    t = np.array([target_coords[r] for r in common_residues])
    
    # Center
    m_centroid = m.mean(axis=0)
    t_centroid = t.mean(axis=0)
    m_c = m - m_centroid
    t_c = t - t_centroid
    
    # SVD
    H = m_c.T @ t_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    return R, t_centroid - m_centroid @ R


def extract_residue_range(atoms, start_resi, end_resi):
    """Extract atoms within residue range."""
    extracted = []
    for line in atoms:
        resi = int(line[22:26].strip())
        if start_resi <= resi <= end_resi:
            extracted.append(line)
    return extracted


def apply_transform(atoms, R, translation):
    """Apply rotation + translation to all atom coordinates."""
    transformed = []
    for line in atoms:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coord = np.array([x, y, z])
        new_coord = coord @ R + translation
        new_line = line[:30] + f"{new_coord[0]:8.3f}" + f"{new_coord[1]:8.3f}" + f"{new_coord[2]:8.3f}" + line[54:]
        transformed.append(new_line)
    return transformed


def change_chain(lines, new_chain):
    """Change chain ID in PDB lines."""
    new_lines = []
    for line in lines:
        new_line = line[:21] + new_chain + line[22:]
        new_lines.append(new_line)
    return new_lines


def merge_pdbs(rec_atoms, lig_atoms, out_file, lig_chain="B"):
    """Merge receptor and ligand into single PDB for Rosetta docking."""
    lig_atoms = change_chain(lig_atoms, lig_chain)
    
    with open(out_file, 'w') as f:
        f.write("REMARK   1 RECEPTOR\n")
        for atom in rec_atoms:
            f.write(atom)
        f.write("TER\n")
        f.write("REMARK   1 LIGAND\n")
        for atom in lig_atoms:
            f.write(atom)
        f.write("TER\n")
        f.write("END\n")
    
    print(f"Merged {len(rec_atoms)} receptor atoms + {len(lig_atoms)} ligand atoms -> {out_file}")


def prepare_4mut_docking(wt_cgas_pdb, mut_cgas_pdb, trim41_pdb, output_pdb,
                         ctd_start=200, ctd_end=554):
    """Prepare docking input for 4mut structure."""
    # Read structures
    wt_atoms = read_pdb_atoms(wt_cgas_pdb)
    mut_atoms = read_pdb_atoms(mut_cgas_pdb)
    trim41_atoms = read_pdb_atoms(trim41_pdb)
    
    # Extract CTD regions
    wt_ctd = extract_residue_range(wt_atoms, ctd_start, ctd_end)
    mut_ctd = extract_residue_range(mut_atoms, ctd_start, ctd_end)
    
    # Parse CA for alignment
    wt_ca = parse_ca_coords(wt_ctd)
    mut_ca = parse_ca_coords(mut_ctd)
    
    common = sorted(set(wt_ca.keys()) & set(mut_ca.keys()))
    print(f"Common CA residues for alignment: {len(common)} (range {common[0]}-{common[-1]})")
    
    # Align mutant CTD to WT CTD
    R, t = kabsch_align(mut_ca, wt_ca, common)
    
    # Apply transform to all mutant CTD atoms
    mut_ctd_aligned = apply_transform(mut_ctd, R, t)
    
    # Merge with TRIM41 (TRIM41 = receptor A, cGAS = ligand B)
    merge_pdbs(trim41_atoms, mut_ctd_aligned, output_pdb, lig_chain="B")
    
    # Verify alignment
    aligned_ca = parse_ca_coords(mut_ctd_aligned)
    rmsd = np.sqrt(np.mean([np.linalg.norm(wt_ca[r] - aligned_ca[r])**2 for r in common]))
    print(f"Alignment RMSD: {rmsd:.3f} Å")


def main():
    parser = argparse.ArgumentParser(description="Prepare 4mut docking input")
    parser.add_argument("--wt-cgas", required=True, help="WT cgas CTD PDB")
    parser.add_argument("--mut-cgas", required=True, help="Mutant cgas monomer PDB")
    parser.add_argument("--trim41", required=True, help="TRIM41 SPRY PDB")
    parser.add_argument("--output", required=True, help="Output PDB")
    parser.add_argument("--start", type=int, default=200, help="CTD start residue")
    parser.add_argument("--end", type=int, default=554, help="CTD end residue")
    args = parser.parse_args()
    
    prepare_4mut_docking(args.wt_cgas, args.mut_cgas, args.trim41, args.output,
                         args.start, args.end)


if __name__ == "__main__":
    main()
