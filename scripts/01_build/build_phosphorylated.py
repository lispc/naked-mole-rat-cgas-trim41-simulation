#!/usr/bin/env python3
"""
Build phosphorylated cGAS-TRIM41 system for MD simulation.

Steps:
  1. Read OpenMM NPT-equilibrated PDB
  2. Add PO3 group to target Ser OG (S305 -> prmtop resid 106)
  3. Rename residue SER -> SEP
  4. Save modified PDB
  5. Run pdb4amber --reduce
  6. Run tleap with leaprc.phosaa19SB

Usage:
    python build_phosphorylated_system.py \
        --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
        --dcd data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_npt.dcd \
        --target-resid 106 \
        --name Hsap_WT_S305phos \
        --outdir data/md_runs/Hsap_WT_S305phos
"""
import argparse
import numpy as np
from pathlib import Path
import MDAnalysis as mda


# Phosphorylated SER (SEP) atom parameters
# P-O(OG) bond length ~1.61 A, O-P-O angle ~109.5 deg (tetrahedral)
P_OG_DIST = 1.61  # Angstrom
O_P_O_ANGLE = 109.5  # degrees


def add_phosphate_group(u, target_resid, chain=None):
    """
    Add PO3 group to a Serine residue.
    
    Returns modified universe with new atoms.
    """
    # Select target residue
    selection = f"resid {target_resid}"
    if chain:
        selection += f" and segid {chain}"
    
    target_res = u.select_atoms(selection)
    if len(target_res) == 0:
        raise ValueError(f"Target residue {target_resid} not found")
    
    resname = target_res.residues[0].resname
    if resname != 'SER':
        raise ValueError(f"Target residue is {resname}, not SER")
    
    # Get key atoms
    og = target_res.select_atoms("name OG")
    cb = target_res.select_atoms("name CB")
    ca = target_res.select_atoms("name CA")
    
    if len(og) == 0 or len(cb) == 0 or len(ca) == 0:
        raise ValueError("Missing OG, CB, or CA atom in target Ser")
    
    og_pos = og.positions[0].copy()
    cb_pos = cb.positions[0].copy()
    ca_pos = ca.positions[0].copy()
    
    print(f"Target: resid {target_resid} {resname}")
    print(f"  OG position: {og_pos}")
    print(f"  CB position: {cb_pos}")
    
    # Calculate P position: extend OG along OG-CB direction
    # P is placed at OG + normalize(OG - CB) * P_OG_DIST
    vec_og_cb = og_pos - cb_pos
    vec_og_cb /= np.linalg.norm(vec_og_cb)
    p_pos = og_pos + vec_og_cb * P_OG_DIST
    
    print(f"  P position: {p_pos}")
    
    # Calculate three O positions around P (tetrahedral geometry)
    # One O is bonded to OG (already there), three new O atoms around P
    # The three O atoms should form a trigonal pyramid with P at apex
    # and OG at the base
    
    # Direction from P to OG (this will be one of the tetrahedral directions)
    dir_p_og = og_pos - p_pos
    dir_p_og /= np.linalg.norm(dir_p_og)
    
    # Find two perpendicular directions for the other two O atoms
    # Use arbitrary vector not parallel to dir_p_og
    arbitrary = np.array([1.0, 0.0, 0.0])
    if np.abs(np.dot(dir_p_og, arbitrary)) > 0.9:
        arbitrary = np.array([0.0, 1.0, 0.0])
    
    perp1 = np.cross(dir_p_og, arbitrary)
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(dir_p_og, perp1)
    perp2 /= np.linalg.norm(perp2)
    
    # For tetrahedral: the three O atoms (O1P, O2P, O3P) are at 109.5° from OG direction
    # Place them symmetrically around the OG direction
    angle_rad = np.radians(180 - O_P_O_ANGLE)  # angle from OG direction
    
    # Three O positions at 120° intervals in the plane perpendicular to OG
    o_positions = []
    for i in range(3):
        theta = 2 * np.pi * i / 3
        # Component along OG direction (towards outside)
        comp_og = np.cos(angle_rad) * dir_p_og
        # Component in perpendicular plane
        comp_perp = np.sin(angle_rad) * (np.cos(theta) * perp1 + np.sin(theta) * perp2)
        o_pos = p_pos + P_OG_DIST * (comp_og + comp_perp)
        o_positions.append(o_pos)
    
    for i, o_pos in enumerate(o_positions):
        print(f"  O{i+1}P position: {o_pos}")
    
    return {
        'p_pos': p_pos,
        'o_positions': o_positions,
        'og_pos': og_pos,
        'resid': target_resid,
        'resname': 'SEP',
    }


def write_phosphorylated_pdb(u, phospho_data, out_pdb):
    """Write PDB with phosphorylated residue."""
    target_resid = phospho_data['resid']
    
    # Determine chain ID: cGAS is chain B (resid 1-218 in prmtop, but in MDAnalysis
    # the segid might be 'SYSTEM'). Use proper 1-char chain IDs.
    # In our prmtop, cGAS is the first 218 residues, TRIM41 is the rest.
    # But MDAnalysis may have lost chain info. Let's infer from residue index.
    
    lines = []
    atom_num = 1
    
    def format_atom_line(num, name, resname, chain, resid, x, y, z, element):
        """Proper PDB ATOM format."""
        return (
            f"ATOM  {num:5d} {name:4s} {resname:3s} {chain:1s}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
    
    # Only write protein atoms (water/ions will be re-added by tleap)
    # This avoids resid > 9999 which breaks PDB format
    prot_atoms = u.select_atoms("protein")
    
    for atom in prot_atoms:
        # Determine chain: first 218 residues = cGAS (chain B), rest = TRIM41 (chain A)
        chain = 'B' if atom.resid <= 218 else 'A'
        
        # Skip HG on target residue (OG bonds to P, not H, in phosphorylated SER)
        if atom.resid == target_resid and atom.name == 'HG':
            continue
        
        # Write original atom
        resname = atom.resname
        if atom.resid == target_resid:
            resname = 'SEP'
        
        lines.append(format_atom_line(
            atom_num, atom.name, resname, chain, atom.resid,
            atom.position[0], atom.position[1], atom.position[2],
            atom.element
        ))
        atom_num += 1
        
        # After OG of target residue, insert P and O atoms
        if atom.resid == target_resid and atom.name == 'OG':
            # P atom
            p_pos = phospho_data['p_pos']
            lines.append(format_atom_line(
                atom_num, 'P', 'SEP', chain, atom.resid,
                p_pos[0], p_pos[1], p_pos[2], 'P'
            ))
            atom_num += 1
            
            # O1P, O2P, O3P
            for i, o_pos in enumerate(phospho_data['o_positions']):
                lines.append(format_atom_line(
                    atom_num, f'O{i+1}P', 'SEP', chain, atom.resid,
                    o_pos[0], o_pos[1], o_pos[2], 'O'
                ))
                atom_num += 1
    
    # Write TER and END
    lines.append("TER\n")
    lines.append("END\n")
    
    with open(out_pdb, 'w') as f:
        f.writelines(lines)
    
    print(f"\nSaved phosphorylated PDB: {out_pdb}")
    print(f"Total atoms: {atom_num - 1}")


def main():
    parser = argparse.ArgumentParser(description='Build phosphorylated cGAS system')
    parser.add_argument('--prmtop', required=True, help='AMBER prmtop file')
    parser.add_argument('--dcd', required=True, help='NPT trajectory (last frame used)')
    parser.add_argument('--target-resid', type=int, default=106, help='Target residue ID in prmtop (106 = S305)')
    parser.add_argument('--name', default='phospho', help='System name')
    parser.add_argument('--outdir', default='.', help='Output directory')
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Load NPT final frame
    print(f"Loading {args.dcd}...")
    u = mda.Universe(args.prmtop, args.dcd)
    u.trajectory[-1]
    
    # Add phosphate group
    print(f"\nAdding phosphate to resid {args.target_resid}...")
    phospho_data = add_phosphate_group(u, args.target_resid)
    
    # Write modified PDB
    raw_pdb = outdir / f"{args.name}_raw.pdb"
    write_phosphorylated_pdb(u, phospho_data, raw_pdb)
    
    print(f"\nNext steps:")
    print(f"  1. pdb4amber -i {raw_pdb} -o {outdir}/{args.name}_amber.pdb --reduce")
    print(f"  2. tleap -f tleap_{args.name}.in")
    print(f"\nTleap input template:")
    print(f"""
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.phosaa19SB

complex = loadPdb {args.name}_amber.pdb
check complex
solvateOct complex OPC 12.0
addIonsRand complex Na+ 0
addIonsRand complex Cl- 0
saveAmberParm complex {args.name}.prmtop {args.name}.rst7
savePdb complex {args.name}_solvated.pdb
charge complex
quit
""")


if __name__ == '__main__':
    main()
