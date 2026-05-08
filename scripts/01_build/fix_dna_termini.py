#!/usr/bin/env python3
"""Fix DNA terminal residue names for Amber tleap (DA→DA5/DA3 etc)."""
import sys
from pathlib import Path

pdb_in = sys.argv[1] if len(sys.argv) > 1 else None
pdb_out = sys.argv[2] if len(sys.argv) > 2 else None

if not pdb_in or not pdb_out:
    print("Usage: python fix_dna_termini.py input.pdb output.pdb")
    sys.exit(1)

with open(pdb_in) as f:
    lines = f.readlines()

# Parse DNA chains: chain -> {resnum: resname}
dna_chains = {}
for line in lines:
    if line.startswith("ATOM") or line.startswith("HETATM"):
        chain = line[21]
        resname = line[17:20].strip()
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if resname in ("DA", "DT", "DG", "DC"):
            dna_chains.setdefault(chain, {})[resnum] = resname

# Determine terminal residues
fixes = {}
for chain, residues in dna_chains.items():
    sorted_nums = sorted(residues.keys())
    if len(sorted_nums) >= 1:
        first = sorted_nums[0]
        fixes[(chain, first)] = residues[first] + "5"
        print(f"Chain {chain}: 5' {residues[first]}{first} → {fixes[(chain, first)]}")
    if len(sorted_nums) >= 2:
        last = sorted_nums[-1]
        fixes[(chain, last)] = residues[last] + "3"
        print(f"Chain {chain}: 3' {residues[last]}{last} → {fixes[(chain, last)]}")

# Apply fixes
out_lines = []
for line in lines:
    if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
        chain = line[21]
        try:
            resnum = int(line[22:26])
        except ValueError:
            out_lines.append(line)
            continue
        if (chain, resnum) in fixes:
            new_name = fixes[(chain, resnum)]
            line = line[:17] + f"{new_name:<3s}" + line[20:]
    out_lines.append(line)

with open(pdb_out, 'w') as f:
    f.writelines(out_lines)

print(f"Written: {pdb_out} ({len(out_lines)} lines)")
