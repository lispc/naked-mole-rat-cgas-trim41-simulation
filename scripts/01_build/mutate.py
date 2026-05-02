#!/usr/bin/env python3
"""
Manually mutate SER B 324 (S305) to GLU in Hsap_WT_amber.pdb.

Uses standard bond lengths and angles for GLU side chain:
  CB-CG = 1.52 A, CG-CD = 1.52 A, CD-OE1 = 1.25 A, CD-OE2 = 1.25 A
  CA-CB-CG ~ 109.5 deg, CB-CG-CD ~ 109.5 deg, CG-CD-OE ~ 120 deg

Usage:
    python mutate_s305e.py \
        --input data/md_runs/Hsap_WT/Hsap_WT_amber.pdb \
        --output data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_mutated.pdb
"""

import argparse
import numpy as np
from pathlib import Path


def normalize(v):
    """Normalize a vector."""
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def place_glu_sidechain(ca_pos, cb_pos):
    """
    Place GLU side chain atoms (CG, CD, OE1, OE2) based on CA and CB positions.
    
    GLU side chain: -CH2-CH2-COO-
    
    Strategy:
      1. CG: extend from CB along the CA->CB direction (trans to CA)
      2. CD: extend from CG, roughly continuing the CB->CG direction with a slight bend
      3. OE1, OE2: planar around CD, ~120 deg apart
    """
    # Standard bond lengths (Angstroms)
    BOND_CB_CG = 1.52
    BOND_CG_CD = 1.52
    BOND_CD_OE = 1.25
    
    # Direction from CA to CB
    dir_ca_cb = normalize(cb_pos - ca_pos)
    
    # CG: extend from CB away from CA
    cg_pos = cb_pos + dir_ca_cb * BOND_CB_CG
    
    # CD: extend from CG. Use a direction slightly different from CB->CG
    # to create the ~109.5 deg angle. We'll introduce a small perpendicular component.
    dir_cb_cg = normalize(cg_pos - cb_pos)
    
    # Create a perpendicular direction for the bend
    # Use cross product with an arbitrary vector to get a perpendicular direction
    arb = np.array([1.0, 0.0, 0.0])
    if np.allclose(np.abs(dir_cb_cg), np.abs(arb)):
        arb = np.array([0.0, 1.0, 0.0])
    perp = normalize(np.cross(dir_cb_cg, arb))
    if np.linalg.norm(perp) < 0.1:
        perp = normalize(np.cross(dir_cb_cg, np.array([0.0, 0.0, 1.0])))
    
    # CD direction: mostly along CB->CG but with a small perpendicular component
    # to achieve ~109.5 deg CB-CG-CD angle
    bend_angle = np.deg2rad(109.5)
    # The component along CB->CG direction
    cos_bend = np.cos(np.pi - bend_angle)  # This gives the angle relative to reverse direction
    # Actually, simpler: CD is placed at ~109.5 deg from CB-CG bond
    # Use tetrahedral geometry: CG is at center, CB and CD are two vertices
    dir_cg_cd = normalize(dir_cb_cg * np.cos(np.deg2rad(70.5)) + perp * np.sin(np.deg2rad(70.5)))
    
    cd_pos = cg_pos + dir_cg_cd * BOND_CG_CD
    
    # OE1, OE2: planar around CD, 120 deg apart
    # Direction from CG to CD
    dir_cg_cd = normalize(cd_pos - cg_pos)
    
    # Perpendicular to CG-CD for the OE plane
    arb2 = np.array([0.0, 1.0, 0.0])
    if np.allclose(np.abs(dir_cg_cd), np.abs(arb2)):
        arb2 = np.array([1.0, 0.0, 0.0])
    perp2 = normalize(np.cross(dir_cg_cd, arb2))
    if np.linalg.norm(perp2) < 0.1:
        perp2 = normalize(np.cross(dir_cg_cd, np.array([0.0, 0.0, 1.0])))
    
    # Another perpendicular in the OE plane
    perp3 = normalize(np.cross(dir_cg_cd, perp2))
    
    # OE1 and OE2 at 120 deg from CG-CD direction
    angle_oe = np.deg2rad(120.0)
    dir_cd_oe1 = normalize(dir_cg_cd * np.cos(np.pi - np.deg2rad(120)) + perp2 * np.sin(np.deg2rad(60)))
    dir_cd_oe2 = normalize(dir_cg_cd * np.cos(np.pi - np.deg2rad(120)) - perp2 * np.sin(np.deg2rad(60)))
    
    oe1_pos = cd_pos + dir_cd_oe1 * BOND_CD_OE
    oe2_pos = cd_pos + dir_cd_oe2 * BOND_CD_OE
    
    return cg_pos, cd_pos, oe1_pos, oe2_pos


def mutate_pdb(input_pdb, output_pdb, chain_id, res_seq, old_resname, new_resname):
    """Mutate a single residue in a PDB file."""
    
    with open(input_pdb) as f:
        lines = f.readlines()
    
    # Extract coordinates of target residue
    ca_pos = None
    cb_pos = None
    target_atoms = []
    other_lines = []
    
    for line in lines:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            other_lines.append(line)
            continue
            
        record = line[0:6].strip()
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21]
        res_num = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        
        if chain == chain_id and res_num == res_seq and res_name == old_resname:
            pos = np.array([x, y, z])
            target_atoms.append((atom_name, pos, line))
            if atom_name == 'CA':
                ca_pos = pos
            elif atom_name == 'CB':
                cb_pos = pos
        else:
            other_lines.append(line)
    
    if ca_pos is None or cb_pos is None:
        raise ValueError(f"Could not find CA/CB for {old_resname} {chain_id} {res_seq}")
    
    print(f"Found target: {old_resname} {chain_id} {res_seq}")
    print(f"  CA: {ca_pos}")
    print(f"  CB: {cb_pos}")
    
    # Place GLU side chain
    cg_pos, cd_pos, oe1_pos, oe2_pos = place_glu_sidechain(ca_pos, cb_pos)
    print(f"  CG: {cg_pos}")
    print(f"  CD: {cd_pos}")
    print(f"  OE1: {oe1_pos}")
    print(f"  OE2: {oe2_pos}")
    
    # Build new residue lines
    new_lines = []
    serial = 1
    
    # First, write all non-target lines
    for line in other_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Update serial number
            new_line = line[:6] + f"{serial:5d}" + line[11:]
            new_lines.append(new_line)
            serial += 1
        else:
            new_lines.append(line)
    
    # Now we need to insert the mutated residue in the correct position
    # Find where to insert (before which line)
    # Actually, simpler approach: rebuild the entire PDB
    
    # Collect all atoms with their positions
    all_atoms = []
    
    for line in lines:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21]
        res_num = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        
        if chain == chain_id and res_num == res_seq and res_name == old_resname:
            # Skip old side chain atoms (OG, HG), keep backbone and CB
            if atom_name in ['OG', 'HG']:
                continue
            # Rename residue
            all_atoms.append((atom_name, new_resname, chain, res_num, np.array([x, y, z])))
        else:
            all_atoms.append((atom_name, res_name, chain, res_num, np.array([x, y, z])))
    
    # Add new GLU side chain atoms
    all_atoms.append(('CG', new_resname, chain_id, res_seq, cg_pos))
    all_atoms.append(('CD', new_resname, chain_id, res_seq, cd_pos))
    all_atoms.append(('OE1', new_resname, chain_id, res_seq, oe1_pos))
    all_atoms.append(('OE2', new_resname, chain_id, res_seq, oe2_pos))
    
    # Sort atoms by chain, res_num, and a reasonable atom order
    # We need to maintain the original order for non-target atoms
    # and insert new atoms after CB of target residue
    
    # Simpler: just write in the order we collected, with new atoms inserted after CB
    ordered_atoms = []
    for atom in all_atoms:
        if atom[0] == 'CB' and atom[2] == chain_id and atom[3] == res_seq:
            ordered_atoms.append(atom)
            ordered_atoms.append(('CG', new_resname, chain_id, res_seq, cg_pos))
            ordered_atoms.append(('CD', new_resname, chain_id, res_seq, cd_pos))
            ordered_atoms.append(('OE1', new_resname, chain_id, res_seq, oe1_pos))
            ordered_atoms.append(('OE2', new_resname, chain_id, res_seq, oe2_pos))
        elif atom[0] not in ['CG', 'CD', 'OE1', 'OE2'] or atom[3] != res_seq:
            ordered_atoms.append(atom)
    
    # Write output PDB (strict PDB format)
    with open(output_pdb, 'w') as f:
        serial = 1
        for atom_name, res_name, chain, res_num, pos in ordered_atoms:
            # Standard PDB ATOM format
            # Col  1- 6: RECORD
            # Col  7-11: serial (5d, right-justified)
            # Col 12   : space
            # Col 13-16: atom name (4 chars, left-justified for 1-char element)
            # Col 17   : altLoc
            # Col 18-20: resName (3 chars, right-justified)
            # Col 21   : space
            # Col 22   : chainID
            # Col 23-26: resSeq (4d, right-justified)
            # Col 27   : iCode
            # Col 28-30: 3 spaces
            # Col 31-38: x (8.3f)
            # Col 39-46: y (8.3f)
            # Col 47-54: z (8.3f)
            # Col 55-60: occupancy (6.2f)
            # Col 61-66: tempFactor (6.2f)
            # Col 67-76: 10 spaces
            # Col 77-78: element (2 chars, right-justified)
            elem = get_element(atom_name)
            # atom name alignment: if 1-char element, left-align in 4 cols; if 2-char, right-align
            if len(elem) == 1 and len(atom_name) <= 3:
                atom_field = f"{atom_name:<4s}"
            else:
                atom_field = f"{atom_name:>4s}"
            
            line = (f"ATOM  {serial:5d} {atom_field} {res_name:>3s} {chain}{res_num:>4d}    "
                    f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                    f"{1.0:6.2f}{0.0:6.2f}          {elem:>2s}\n")
            f.write(line)
            serial += 1
        f.write("END\n")
    
    print(f"Wrote mutated PDB: {output_pdb}")
    print(f"  Total atoms: {serial - 1}")


def get_element(atom_name):
    """Get element from atom name."""
    if atom_name.startswith('C'):
        return 'C'
    elif atom_name.startswith('N'):
        return 'N'
    elif atom_name.startswith('O'):
        return 'O'
    elif atom_name.startswith('S'):
        return 'S'
    elif atom_name.startswith('H'):
        return 'H'
    elif atom_name.startswith('P'):
        return 'P'
    else:
        return atom_name[0]


def main():
    parser = argparse.ArgumentParser(description='Mutate SER to GLU in PDB')
    parser.add_argument('--input', required=True, help='Input PDB')
    parser.add_argument('--output', required=True, help='Output PDB')
    parser.add_argument('--chain', default='B', help='Chain ID')
    parser.add_argument('--resseq', type=int, default=324, help='Residue sequence number')
    parser.add_argument('--old', default='SER', help='Original residue name')
    parser.add_argument('--new', default='GLU', help='New residue name')
    args = parser.parse_args()
    
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    mutate_pdb(args.input, args.output, args.chain, args.resseq, args.old, args.new)


if __name__ == '__main__':
    main()
