#!/usr/bin/env python3
"""
Build E2~Ub-TRIM25_RING + cGAS quaternary model v2.

Key improvements over v1:
1. Include TRIM25 RING dimer (5FER chains A,D) to stabilize E2~Ub closed conformation
2. Isopeptide bond (E2 K85 NZ ↔ Ub G76 C) handled by harmonic restraint in MD
3. Same cGAS backbone for WT and 4mut (4mut = D431S/K479E/L495Y/K498T)

Component ordering in output PDB:
  RING1, RING2, E2, Ub, cGAS
"""
import sys
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Structure, Model

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/structures/quaternary_minimal"
OUTDIR.mkdir(parents=True, exist_ok=True)

parser = PDBParser(QUIET=True)
TARGET_DIST = 15.0


def norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else v


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

n_ring1 = len([r for r in ring1 if r.id[0] == ' '])
n_ring2 = len([r for r in ring2 if r.id[0] == ' '])
n_e2 = len([r for r in e2 if r.id[0] == ' '])
n_ub = len([r for r in ub if r.id[0] == ' '])
print(f"  RING1: {n_ring1}, RING2: {n_ring2}, E2: {n_e2}, Ub: {n_ub}")

# Find Ub G76 C and E2 K85 NZ
ub_g76_c = None
e2_k85_nz = None
for res in ub:
    if res.resname == 'GLY' and res.id[1] == 76 and 'C' in res:
        ub_g76_c = res['C'].coord.copy()
for res in e2:
    if res.resname == 'LYS' and res.id[1] == 85 and 'NZ' in res:
        e2_k85_nz = res['NZ'].coord.copy()
assert ub_g76_c is not None and e2_k85_nz is not None
print(f"  Ub G76 C: {ub_g76_c}")
print(f"  E2 K85 NZ: {e2_k85_nz}")
print(f"  Isopeptide distance: {np.linalg.norm(ub_g76_c - e2_k85_nz):.2f} Å")

catalytic_center = ub_g76_c.copy()

# Collect CA atoms from RING+E2+Ub for clash check
ring_e2_ub_cas = []
for chain_obj in [ring1, ring2, e2, ub]:
    for res in chain_obj:
        if res.id[0] == ' ' and 'CA' in res:
            ring_e2_ub_cas.append(res['CA'].coord.copy())


def count_clashes(cgas_chain):
    n = 0
    for res in cgas_chain:
        if res.id[0] != ' ' or 'CA' not in res:
            continue
        for ca in ring_e2_ub_cas:
            if np.linalg.norm(res['CA'].coord - ca) < 3.5:
                n += 1
    return n


def resolve_clashes(cgas_chain, max_iter=100):
    for it in range(max_iter):
        n = count_clashes(cgas_chain)
        if n == 0:
            return True, it
        vec = np.zeros(3)
        for res in cgas_chain:
            if res.id[0] != ' ' or 'CA' not in res:
                continue
            for ca in ring_e2_ub_cas:
                d = np.linalg.norm(res['CA'].coord - ca)
                if d < 3.5 and d > 0.01:
                    vec += (res['CA'].coord - ca) / d * (3.5 - d)
        nrm = np.linalg.norm(vec)
        if nrm < 0.001:
            return True, it
        vec = vec / nrm * min(nrm, 2.0)
        for atom in cgas_chain.get_atoms():
            atom.coord += vec
    return False, max_iter


# ============================================================================
# STEP 2: Build WT and 4mut
# ============================================================================
print("\nSTEP 2: Building models")

cfgs = {
    "WT": BASE / 'structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb',
    "4mut": OUTDIR / 'cgas_4mut.pdb',
}

for label, cgas_path in cfgs.items():
    print(f"\n{'='*50}")
    print(f"  {label}")
    print(f"{'='*50}")

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
    assert k315_nz is not None

    # Position
    vec = catalytic_center - k315_nz
    current_dist = np.linalg.norm(vec)
    print(f"  Initial K315-NZ → Ub-G76-C: {current_dist:.1f} Å")

    vec_unit = vec / current_dist
    for atom in cgas.get_atoms():
        atom.coord += vec_unit * (current_dist - TARGET_DIST)

    new_dist = np.linalg.norm(
        [r['NZ'].coord for r in cgas if r.resname == 'LYS' and r.id[1] == 315 and 'NZ' in r][0]
        - catalytic_center)
    print(f"  After translation: {new_dist:.1f} Å")

    n = count_clashes(cgas)
    print(f"  CA clashes: {n}")
    if n > 0:
        ok, it = resolve_clashes(cgas)
        print(f"  After resolution ({it} iter): {count_clashes(cgas)} clashes")
        dist = np.linalg.norm(
            [r['NZ'].coord for r in cgas if r.resname == 'LYS' and r.id[1] == 315 and 'NZ' in r][0]
            - catalytic_center)
        print(f"  K315-NZ → Ub-G76-C after resolve: {dist:.1f} Å")

    # Merge
    s = Structure.Structure(f'quaternary_{label}')
    m = Model.Model(0)
    s.add(m)
    for chain_obj in [ring1, ring2, e2, ub, cgas]:
        m.add(chain_obj.copy())

    io = PDBIO()
    io.set_structure(s)
    out_pdb = OUTDIR / f'quaternary_ring_{label}_raw.pdb'
    io.save(str(out_pdb))

    final_dist = np.linalg.norm(
        [r['NZ'].coord for r in cgas if r.resname == 'LYS' and r.id[1] == 315 and 'NZ' in r][0]
        - catalytic_center)
    print(f"  >>> {label}: K315→Ub={final_dist:.1f} Å, clashes={count_clashes(cgas)}")
    print(f"  Saved: {out_pdb}")

print("\n" + "=" * 60)
print("BUILD COMPLETE.")
print("Next: pdb4amber → tleap → MD with isopeptide bond restraint")
print("=" * 60)
