#!/usr/bin/env python3
"""
Build FULL quaternary complex: RING-E2~Ub + SPRY-cGAS.

Strategy (v3):
1. RING1, RING2, E2, Ub from 5FER (validated in v2)
2. cGAS from AF3, positioned with K315 @ 15 Å from Ub-G76 (validated in v2)
3. SPRY from Rosetta docking best pose, aligned to cGAS via Kabsch on cGAS CA atoms
   (preserves the cGAS-SPRY interface from Rosetta)
4. CC linker restraint: RING C-term CA ↔ SPRY N-term CA (flat-bottom 80-120 Å)

Chains in output: R(RING1), S(RING2), E(E2), U(Ub), P(SPRY), C(cGAS)
"""
import sys
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Structure, Model

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/structures/quaternary_full"
OUTDIR.mkdir(parents=True, exist_ok=True)

parser = PDBParser(QUIET=True)
TARGET_DIST = 15.0  # K315 NZ → Ub G76 C


def norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else v


def kabsch(A, B):
    """Optimal rotation matrix to align A onto B (both N×3)."""
    A_c = A - A.mean(axis=0)
    B_c = B - B.mean(axis=0)
    H = A_c.T @ B_c
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    return R


# ============================================================================
# STEP 1: Load RING-E2~Ub from 5FER
# ============================================================================
print("STEP 1: Loading RING-E2~Ub from 5FER")
fer = parser.get_structure('5fer', str(BASE / 'data/structures/quaternary_mvp/5FER.pdb'))

ring1 = fer[0]['A']
ring2 = fer[0]['D']
e2 = fer[0]['B']
ub = fer[0]['C']

for res in ring1:
    res.parent.id = 'R'
for res in ring2:
    res.parent.id = 'S'
for res in e2:
    res.parent.id = 'E'
for res in ub:
    res.parent.id = 'U'

# Find Ub G76 C
ub_g76_c = None
for res in ub:
    if res.resname == 'GLY' and res.id[1] == 76 and 'C' in res:
        ub_g76_c = res['C'].coord.copy()
        break
assert ub_g76_c is not None
print(f"  Ub G76 C: {ub_g76_c}")

# Find RING1 C-term CA (last residue before chain A ends)
ring1_cas = [(r.id[1], r['CA'].coord.copy()) for r in ring1 if r.id[0] == ' ' and 'CA' in r]
ring1_cterm_ca = ring1_cas[-1][1]
print(f"  RING1 C-term: res {ring1_cas[-1][0]}, CA {ring1_cterm_ca}")

# Collect RING+E2+Ub CAs for clash check
ring_e2_ub_cas = []
for chain_obj in [ring1, ring2, e2, ub]:
    for res in chain_obj:
        if res.id[0] == ' ' and 'CA' in res:
            ring_e2_ub_cas.append(res['CA'].coord.copy())


def count_clashes(cgas_chain, spry_chain=None):
    """Count CA-CA clashes between cGAS/SPRY and RING-E2-Ub."""
    n = 0
    targets = list(ring_e2_ub_cas)
    for chain_obj in [cgas_chain]:
        if chain_obj is None:
            continue
        for res in chain_obj:
            if res.id[0] != ' ' or 'CA' not in res:
                continue
            for ca in targets:
                if np.linalg.norm(res['CA'].coord - ca) < 3.5:
                    n += 1
    # Also check SPRY vs others
    if spry_chain is not None:
        for res in spry_chain:
            if res.id[0] != ' ' or 'CA' not in res:
                continue
            for ca in targets:
                if np.linalg.norm(res['CA'].coord - ca) < 3.5:
                    n += 1
    return n


# ============================================================================
# STEP 2: Load Rosetta SPRY-cGAS pose (for SPRY-cGAS interface)
# ============================================================================
print("\nSTEP 2: Loading Rosetta SPRY-cGAS pose")
pose = parser.get_structure('pose',
                             str(BASE / 'structures/docking/rosetta/output_hsap_WT_global/hsap_WT_input_0081.pdb'))
spry_pose = pose[0]['A']
cgas_pose = pose[0]['B']

# Extract cGAS CA coords from Rosetta pose
cgas_pose_cas = []
for res in cgas_pose:
    if res.id[0] == ' ' and 'CA' in res:
        cgas_pose_cas.append(res['CA'].coord.copy())
cgas_pose_cas = np.array(cgas_pose_cas)
print(f"  Rosetta cGAS: {len(cgas_pose_cas)} CA atoms")

# Extract all SPRY atoms for transformation
spry_atoms = list(spry_pose.get_atoms())
spry_coords = np.array([a.coord.copy() for a in spry_atoms])
print(f"  Rosetta SPRY: {len(spry_atoms)} atoms")

# SPRY N-term CA
spry_nterm_ca = None
for res in spry_pose:
    if res.id[0] == ' ' and 'CA' in res:
        spry_nterm_ca = res['CA'].coord.copy()
        print(f"  SPRY N-term: res {res.id[1]}, CA {spry_nterm_ca}")
        break

# ============================================================================
# STEP 3: Build WT and 4mut models
# ============================================================================
print("\nSTEP 3: Building full quaternary models")

cfgs = {
    "WT": BASE / 'structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb',
    "4mut": BASE / 'data/structures/quaternary_minimal/cgas_4mut.pdb',
}

for label, cgas_path in cfgs.items():
    print(f"\n{'='*60}")
    print(f"  FULL QUATERNARY: {label}")
    print(f"{'='*60}")

    # Load cGAS from AF3
    cgas_struct = parser.get_structure(f'cgas_{label}', str(cgas_path))
    model = cgas_struct[0]
    cgas = model[list(model.child_dict.keys())[0]]
    for res in cgas:
        res.parent.id = 'C'

    # Find K315
    k315_nz = None
    for res in cgas:
        if res.resname == 'LYS' and res.id[1] == 315 and 'NZ' in res:
            k315_nz = res['NZ'].coord.copy()
            break
    assert k315_nz is not None, f"K315 NZ not found in {label}"

    # --- Position cGAS at target distance from Ub-G76 ---
    vec = ub_g76_c - k315_nz
    current_dist = np.linalg.norm(vec)
    print(f"  Initial K315-NZ → Ub-G76-C: {current_dist:.1f} Å")

    vec_unit = vec / current_dist
    translation = vec_unit * (current_dist - TARGET_DIST)
    for atom in cgas.get_atoms():
        atom.coord += translation

    new_dist = np.linalg.norm(
        [r['NZ'].coord for r in cgas if r.resname == 'LYS' and r.id[1] == 315 and 'NZ' in r][0]
        - ub_g76_c)
    print(f"  After positioning: {new_dist:.1f} Å")

    # --- Extract positioned cGAS CA coords ---
    cgas_new_cas = []
    for res in cgas:
        if res.id[0] == ' ' and 'CA' in res:
            cgas_new_cas.append(res['CA'].coord.copy())
    cgas_new_cas = np.array(cgas_new_cas)

    # --- Align SPRY to cGAS ---
    # Compute Kabsch: align Rosetta cGAS CA → positioned cGAS CA
    assert len(cgas_pose_cas) == len(cgas_new_cas), \
        f"CA count mismatch: {len(cgas_pose_cas)} vs {len(cgas_new_cas)}"

    R = kabsch(cgas_pose_cas, cgas_new_cas)
    t_pose = cgas_pose_cas.mean(axis=0)
    t_new = cgas_new_cas.mean(axis=0)

    # Apply to SPRY: center, rotate, translate
    spry_new_coords = (spry_coords - t_pose) @ R.T + t_new
    for i, atom in enumerate(spry_atoms):
        atom.coord = spry_new_coords[i]

    # Rename SPRY chain
    for res in spry_pose:
        res.parent.id = 'P'

    # Verify SPRY-cGAS interface preserved
    spry_new_nterm = spry_new_coords[0]  # first SPRY atom (N-term N)
    cc_dist = np.linalg.norm(ring1_cterm_ca - spry_new_nterm)
    print(f"  RING1 C-term CA → SPRY N-term N: {cc_dist:.0f} Å")

    # Clash check
    n_cl = count_clashes(cgas, spry_pose)
    print(f"  Total CA clashes with RING+E2+Ub: {n_cl}")

    # --- Merge all chains ---
    s = Structure.Structure(f'full_quaternary_{label}')
    m = Model.Model(0)
    s.add(m)
    for chain_obj in [ring1, ring2, e2, ub, spry_pose, cgas]:
        m.add(chain_obj.copy())

    io = PDBIO()
    io.set_structure(s)
    out_pdb = OUTDIR / f'full_quaternary_{label}_raw.pdb'
    io.save(str(out_pdb))

    # Final check
    k315_final = None
    for res in cgas:
        if res.resname == 'LYS' and res.id[1] == 315 and 'NZ' in res:
            k315_final = res['NZ'].coord
            break
    final_k315_ub = np.linalg.norm(k315_final - ub_g76_c) if k315_final is not None else -1
    print(f"  >>> {label} FINAL: K315→Ub={final_k315_ub:.1f}Å, CC={cc_dist:.0f}Å, clashes={n_cl}")
    print(f"  Saved: {out_pdb}")

print("\n" + "=" * 60)
print("BUILD COMPLETE")
print(f"\nChains: R(RING1), S(RING2), E(E2), U(Ub), P(SPRY), C(cGAS)")
print(f"CC linker distance: ~{cc_dist:.0f} Å (target 80-120 Å for restraint)")
print(f"\nNext: pdb4amber → tleap → MD with:")
print(f"  1. Isopeptide restraint (K85 NZ ↔ G76 C, k=5000, r0=1.35 Å)")
print(f"  2. CC linker restraint (RING1 C-term CA ↔ SPRY N-term CA, flat-bottom 80-120 Å)")
print(f"  3. COM flat-bottom (RING+E2+Ub+SPRY ↔ cGAS)")
print("=" * 60)
