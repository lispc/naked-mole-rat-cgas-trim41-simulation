#!/usr/bin/env python3
"""Compare cGAS N-terminal conformation: apo (AF3/docked) vs DNA-bound (4LEZ)."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent.parent

u_dna = mda.Universe(str(BASE / "structures/cgas_dna/4LEZ.pdb"))
u_dock = mda.Universe(str(BASE / "structures/docking/rosetta/output_global/input_0001.pdb"))

mobile_ca = u_dock.select_atoms("segid B and name CA")  # cGAS in docked complex
ref_ca = u_dna.select_atoms("segid A and name CA")      # cGAS in DNA-bound

common = sorted(set(mobile_ca.resids) & set(ref_ca.resids))
print(f"Common CA residues: {len(common)}, range {common[0]}-{common[-1]}")

sel = f"name CA and resid {' '.join(map(str, common))}"
rmsd_pre = rms.rmsd(
    u_dock.select_atoms(f"segid B and {sel}").positions,
    u_dna.select_atoms(f"segid A and {sel}").positions)
print(f"RMSD before alignment: {rmsd_pre:.2f} Å")

result = align.alignto(mobile_ca, ref_ca, select=sel)
print(f"Alignment result type: {type(result)}, len: {len(result)}, values: {result}")

for region, (lo, hi) in [("N-term 200-299", (200, 299)),
                           ("Core 300-399", (300, 399)),
                           ("C-term 400-507", (400, 507))]:
    rc = [r for r in range(lo, hi + 1) if r in common]
    if rc:
        s = f"name CA and resid {' '.join(map(str, rc))}"
        rv = rms.rmsd(
            u_dock.select_atoms(f"segid B and {s}").positions,
            u_dna.select_atoms(f"segid A and {s}").positions)
        print(f"  {region}: {rv:.2f} Å ({len(rc)} residues)")

# Key interface residues (paper Table 2: 211-219)
print("\nKey N-term interface residues (211-219):")
for res in range(211, 220):
    if res in common:
        dp = u_dock.select_atoms(f"segid B and resid {res} and name CA").positions
        np_ = u_dna.select_atoms(f"segid A and resid {res} and name CA").positions
        print(f"  Res {res}: {np.linalg.norm(dp[0] - np_[0]):.2f} Å")

# Check: does the N-term in 4LEZ form a different structure due to DNA binding?
# 4LEZ is a dimer - check distance between N-term and the other cGAS chain
print("\n4LEZ dimer N-term across chains:")
nterm_a = u_dna.select_atoms("segid A and resid 200-299 and name CA")
nterm_c = u_dna.select_atoms("segid C and resid 200-299 and name CA")
if len(nterm_c) > 0:
    com_a = nterm_a.center_of_mass()
    com_c = nterm_c.center_of_mass()
    print(f"  Chain A N-term COM -> Chain C N-term COM: {np.linalg.norm(com_a - com_c):.1f} Å")
