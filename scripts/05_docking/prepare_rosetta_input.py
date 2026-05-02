#!/usr/bin/env python3
"""Prepare Rosetta docking input by merging receptor and ligand PDBs."""

import sys
import argparse

def read_pdb_atoms(filepath):
    """Read PDB ATOM/HETATM records."""
    atoms = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(line)
    return atoms

def change_chain(lines, new_chain):
    """Change chain ID in PDB lines."""
    new_lines = []
    for line in lines:
        # PDB format: chain ID at position 22 (0-indexed: 21)
        new_line = line[:21] + new_chain + line[22:]
        new_lines.append(new_line)
    return new_lines

def merge_pdbs(rec_file, lig_file, out_file, lig_chain="B", translate=None):
    """Merge receptor and ligand into single PDB for Rosetta docking."""
    rec_atoms = read_pdb_atoms(rec_file)
    lig_atoms = read_pdb_atoms(lig_file)
    
    # Change ligand chain
    lig_atoms = change_chain(lig_atoms, lig_chain)
    
    with open(out_file, 'w') as f:
        # Write receptor
        f.write("REMARK   1 RECEPTOR\n")
        for atom in rec_atoms:
            f.write(atom)
        f.write("TER\n")
        
        # Write ligand
        f.write("REMARK   1 LIGAND\n")
        for atom in lig_atoms:
            if translate:
                # Parse coordinates and translate
                x = float(atom[30:38]) + translate[0]
                y = float(atom[38:46]) + translate[1]
                z = float(atom[46:54]) + translate[2]
                new_atom = atom[:30] + f"{x:8.3f}" + f"{y:8.3f}" + f"{z:8.3f}" + atom[54:]
                f.write(new_atom)
            else:
                f.write(atom)
        f.write("TER\n")
        f.write("END\n")
    
    print(f"Merged {len(rec_atoms)} receptor atoms + {len(lig_atoms)} ligand atoms -> {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Prepare Rosetta docking input")
    parser.add_argument("--receptor", required=True, help="Receptor PDB")
    parser.add_argument("--ligand", required=True, help="Ligand PDB")
    parser.add_argument("--output", required=True, help="Output PDB")
    parser.add_argument("--ligand-chain", default="B", help="Ligand chain ID")
    parser.add_argument("--translate", nargs=3, type=float, default=None, 
                        help="Translate ligand by (x, y, z)")
    args = parser.parse_args()
    
    merge_pdbs(args.receptor, args.ligand, args.output, 
               args.ligand_chain, args.translate)

if __name__ == "__main__":
    main()
