#!/usr/bin/env python3
"""
Build E2~Ub-TRIM41-cGAS quaternary complex MVP.

Strategy:
1. Extract RING dimer + E2~Ub from PDB 5FER (TRIM25 RING + UBE2D1-Ub)
2. Use TRIM25 RING as proxy for TRIM41 RING (highly conserved fold)
3. Load cGAS-TRIM41(SPRY) from our best Rosetta docking pose
4. Translate cGAS-SPRY so K315 is ~15 Å from Ub-G76 catalytic center
5. Output PDB for MD preparation

The isopeptide mimic in 5FER uses UBE2D1 K85 linked to Ub C-terminus.
"""

import sys
import copy
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBIO, PDBParser

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/structures/quaternary_mvp"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# STEP 1: Load template structures
# ============================================================================
print("=" * 60)
print("STEP 1: Loading template structures")
print("=" * 60)

parser = PDB.PDBParser(QUIET=True)

# PDB 5FER: TRIM25 RING dimer + 2x(UBE2D1~Ub)
fer = parser.get_structure('5fer', BASE / 'data/structures/quaternary_mvp/5FER.pdb')

# Extract one RING dimer + one E2~Ub
# 5FER chains: A=RING1, D=RING2, B=E2, C=Ub
ring1 = fer[0]['A'].copy()
ring2 = fer[0]['D'].copy()
e2 = fer[0]['B'].copy()
ub = fer[0]['C'].copy()

print(f"RING1 (Chain A): {len([r for r in ring1 if r.id[0]==' '])} aa residues")
print(f"RING2 (Chain D): {len([r for r in ring2 if r.id[0]==' '])} aa residues")
print(f"E2 UBE2D1 (Chain B): {len([r for r in e2 if r.id[0]==' '])} aa residues")
print(f"Ub (Chain C): {len([r for r in ub if r.id[0]==' '])} aa residues")

# Verify key residues
# Ub G76
ub_g76 = None
for res in ub:
    if res.resname == 'GLY' and res.id[1] == 76:
        ub_g76 = res
        break
assert ub_g76 is not None, "Ub G76 not found"
ub_g76_c = ub_g76['C'].coord
print(f"Ub G76 C: {ub_g76_c}")

# E2 K85 (isopeptide mimic linkage point)
e2_k85 = None
for res in e2:
    if res.resname == 'LYS' and res.id[1] == 85:
        e2_k85 = res
        break
if e2_k85 is not None:
    print(f"E2 K85 NZ: {e2_k85['NZ'].coord}")
    catalytic_center = e2_k85['NZ'].coord  # isopeptide N is the transfer point
else:
    print("WARNING: E2 K85 not found, using Ub G76 C as catalytic center")
    catalytic_center = ub_g76_c

# ============================================================================
# STEP 2: Load cGAS-TRIM41(SPRY) docking pose
# ============================================================================
print("\n" + "=" * 60)
print("STEP 2: Loading cGAS-SPRY docking pose")
print("=" * 60)

pose_path = BASE / 'structures/docking/rosetta/output_hsap_WT_global/hsap_WT_input_0081.pdb'
pose = parser.get_structure('pose', pose_path)

spry = pose[0]['A'].copy()
cgas = pose[0]['B'].copy()

print(f"SPRY: {len([r for r in spry if r.id[0]==' '])} aa residues")
print(f"cGAS: {len([r for r in cgas if r.id[0]==' '])} aa residues")

# Verify K315 exists
k315 = None
for res in cgas:
    if res.resname == 'LYS' and res.id[1] == 315:
        k315 = res
        break
assert k315 is not None, "K315 not found in cGAS"

k315_nz = k315['NZ'].coord
k315_ca = k315['CA'].coord
print(f"K315 CA: {k315_ca}")
print(f"K315 NZ: {k315_nz}")

# ============================================================================
# STEP 3: Rename chains
# ============================================================================
print("\n" + "=" * 60)
print("STEP 3: Renaming chains")
print("=" * 60)

for res in ring1.get_residues(): res.parent.id = 'R'
for res in ring2.get_residues(): res.parent.id = 'S'
for res in e2.get_residues(): res.parent.id = 'E'
for res in ub.get_residues(): res.parent.id = 'U'
for res in spry.get_residues(): res.parent.id = 'P'
for res in cgas.get_residues(): res.parent.id = 'C'

print("Chains: R=RING1, S=RING2, E=E2, U=Ub, P=SPRY, C=cGAS")

# ============================================================================
# STEP 4: Position cGAS-SPRY relative to catalytic center
# ============================================================================
print("\n" + "=" * 60)
print("STEP 4: Positioning cGAS-SPRY")
print("=" * 60)

# Current distance
init_dist = np.linalg.norm(k315_nz - catalytic_center)
print(f"Initial K315-NZ -> catalytic center distance: {init_dist:.1f} Å")

# Target distance: ~15 Å (gives room for side chain motion and avoids clashes)
# Literature: catalytic distance for nucleophilic attack is ~3-5 Å, but
# initial modeling target should be larger to allow MD relaxation
TARGET_DIST = 15.0

# Compute translation vector
vec = catalytic_center - k315_nz
current_dist = np.linalg.norm(vec)
translation = vec * (1 - TARGET_DIST / current_dist)

print(f"Translation vector: {translation}")

# Apply translation to cGAS and SPRY
for atom in cgas.get_atoms():
    atom.coord += translation
for atom in spry.get_atoms():
    atom.coord += translation

# Verify new distance
k315_nz_new = k315['NZ'].coord
new_dist = np.linalg.norm(k315_nz_new - catalytic_center)
print(f"After translation: K315-NZ -> catalytic center distance: {new_dist:.1f} Å")

# ============================================================================
# STEP 5: Simple orientation adjustment (optional, gentle rotation)
# ============================================================================
print("\n" + "=" * 60)
print("STEP 5: Optional orientation")
print("=" * 60)

# Compute vector from K315 CA to NZ (side chain direction)
k315_ca_new = k315['CA'].coord
k315_nz_new = k315['NZ'].coord
sc_vec = k315_nz_new - k315_ca_new
sc_vec /= np.linalg.norm(sc_vec)

# Desired direction: toward catalytic center
target_vec = catalytic_center - k315_ca_new
target_vec /= np.linalg.norm(target_vec)

# Angle between current and desired
cos_angle = np.dot(sc_vec, target_vec)
angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
print(f"Angle between K315 side chain and catalytic center: {angle:.1f}°")

# Only rotate if angle > 45° (gentle threshold)
if angle > 45.0:
    print("Applying gentle rotation...")
    # Rotation axis
    rot_axis = np.cross(sc_vec, target_vec)
    rot_norm = np.linalg.norm(rot_axis)
    if rot_norm > 0.001:
        rot_axis /= rot_norm
        
        # Build rotation matrix (Rodrigues)
        # Use scipy if available, otherwise manual
        try:
            from scipy.spatial.transform import Rotation as R
            rotation = R.from_rotvec(rot_axis * np.radians(min(angle, 60)))  # cap at 60°
            
            # Collect all cGAS + SPRY atom coords
            all_atoms = list(cgas.get_atoms()) + list(spry.get_atoms())
            coords = np.array([a.coord for a in all_atoms])
            
            # Rotate around K315 CA
            centered = coords - k315_ca_new
            rotated = rotation.apply(centered)
            new_coords = rotated + k315_ca_new
            
            for i, atom in enumerate(all_atoms):
                atom.coord = new_coords[i]
            
            # Verify
            sc_vec_new = k315['NZ'].coord - k315['CA'].coord
            sc_vec_new /= np.linalg.norm(sc_vec_new)
            new_angle = np.degrees(np.arccos(np.clip(np.dot(sc_vec_new, target_vec), -1.0, 1.0)))
            print(f"After rotation: angle = {new_angle:.1f}°")
            
        except ImportError:
            print("scipy not available, skipping rotation")
else:
    print("Orientation acceptable, no rotation needed")

# ============================================================================
# STEP 6: Merge into single structure
# ============================================================================
print("\n" + "=" * 60)
print("STEP 6: Merging structure")
print("=" * 60)

s = PDB.Structure.Structure('quaternary')
m = PDB.Model.Model(0)
s.add(m)

for chain_obj in [ring1, ring2, e2, ub, spry, cgas]:
    m.add(chain_obj)

# Final geometry check
k315_nz_final = k315['NZ'].coord
final_dist = np.linalg.norm(k315_nz_final - catalytic_center)
print(f"\nFinal geometry:")
print(f"  K315-NZ -> catalytic center: {final_dist:.1f} Å")
print(f"  K315-NZ -> Ub-G76-C: {np.linalg.norm(k315_nz_final - ub_g76_c):.1f} Å")

# Clash check (CA-CA distance < 3.5 Å between different chains)
print(f"\nClash check (inter-chain CA-CA < 3.5 Å):")
chain_residues = {}
for chain_obj in [ring1, ring2, e2, ub, spry, cgas]:
    ca_list = []
    for res in chain_obj:
        if res.id[0] == ' ' and 'CA' in res:
            ca_list.append((res['CA'].coord, res.resname + str(res.id[1])))
    chain_residues[chain_obj.id] = ca_list

n_clashes = 0
clash_pairs = []
chain_ids = list(chain_residues.keys())
for i, c1 in enumerate(chain_ids):
    for c2 in chain_ids[i+1:]:
        for coord1, name1 in chain_residues[c1]:
            for coord2, name2 in chain_residues[c2]:
                dist = np.linalg.norm(coord1 - coord2)
                if dist < 3.5:
                    n_clashes += 1
                    if len(clash_pairs) < 10:
                        clash_pairs.append(f"    {c1}:{name1} - {c2}:{name2}: {dist:.2f} Å")

print(f"  Total inter-chain CA clashes: {n_clashes}")
for cp in clash_pairs:
    print(cp)

# Save raw PDB
io = PDBIO()
io.set_structure(s)
raw_pdb = OUTDIR / 'quaternary_mvp_raw.pdb'
io.save(str(raw_pdb))
print(f"\nSaved raw PDB: {raw_pdb}")

# ============================================================================
# STEP 7: Process with pdb4amber
# ============================================================================
print("\n" + "=" * 60)
print("STEP 7: Processing with pdb4amber")
print("=" * 60)

import subprocess

clean_pdb = OUTDIR / 'quaternary_mvp_clean.pdb'
result = subprocess.run(
    ['pdb4amber', '-i', str(raw_pdb), '-o', str(clean_pdb), '-y'],
    capture_output=True, text=True
)
print(result.stdout)
if result.returncode != 0:
    print(f"pdb4amber stderr: {result.stderr}")
else:
    print(f"Saved clean PDB: {clean_pdb}")

print("\n" + "=" * 60)
print("MVP structure preparation complete!")
print("=" * 60)
print(f"\nNext steps:")
print(f"  1. Run tleap to add hydrogens and generate topology/params")
print(f"  2. Run energy minimization")
print(f"  3. Start 50ns test MD on GPU 2 or 3")
