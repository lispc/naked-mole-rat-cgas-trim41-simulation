#!/usr/bin/env python3
"""
Build minimal E2~Ub + cGAS quaternary models for MD.

WT: E2~Ub (5FER) + WT cGAS AF3 (K315 @ 15Å from Ub-G76)
4mut: E2~Ub (5FER) + mutated cGAS (D431S/K479E/L495Y/K498T, K315 @ 15Å)
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


# ---- Load E2~Ub from 5FER ----
print("Loading E2~Ub from 5FER")
fer = parser.get_structure('5fer', str(BASE / 'data/structures/quaternary_mvp/5FER.pdb'))
e2 = fer[0]['B']
ub = fer[0]['C']
for res in e2:
    res.parent.id = 'E'
for res in ub:
    res.parent.id = 'U'

ub_g76_c = None
for res in ub:
    if res.resname == 'GLY' and res.id[1] == 76 and 'C' in res:
        ub_g76_c = res['C'].coord.copy()
        break
assert ub_g76_c is not None
print(f"Ub G76 C: {ub_g76_c}")

# Collect E2/Ub CA for clash check
e2_ub_cas = []
for chain_obj in [e2, ub]:
    for res in chain_obj:
        if res.id[0] == ' ' and 'CA' in res:
            e2_ub_cas.append(res['CA'].coord.copy())


def count_clashes(cgas_chain):
    n = 0
    for res in cgas_chain:
        if res.id[0] != ' ' or 'CA' not in res:
            continue
        for ca in e2_ub_cas:
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
            for ca in e2_ub_cas:
                d = np.linalg.norm(res['CA'].coord - ca)
                if d < 3.5 and d > 0.01:
                    vec += (res['CA'].coord - ca) / d * (3.5 - d)
        nrm = np.linalg.norm(vec)
        if nrm < 0.001:
            return True, it
        step = min(nrm, 2.0)
        vec = vec / nrm * step
        for atom in cgas_chain.get_atoms():
            atom.coord += vec
    return False, max_iter


# ---- Build for each variant ----
cfgs = {
    "WT": BASE / 'structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb',
    "4mut": OUTDIR / 'cgas_4mut.pdb',  # raw mutated PDB (no pdb4amber yet)
}

for label, cgas_path in cfgs.items():
    print(f"\n{'='*50}")
    print(f"Building {label}")
    print(f"{'='*50}")

    cgas_struct = parser.get_structure(f'cgas_{label}', str(cgas_path))
    # Handle both chain A and blank/space chain IDs
    model = cgas_struct[0]
    chain_ids = list(model.child_dict.keys())
    cgas = model[chain_ids[0]]
    for res in cgas:
        res.parent.id = 'C'

    # Find K315
    k315_nz = None
    for res in cgas:
        if res.resname == 'LYS' and res.id[1] == 315 and 'NZ' in res:
            k315_nz = res['NZ'].coord.copy()
            break
    assert k315_nz is not None, f"K315 NZ not found in {label}"

    # Position
    vec = ub_g76_c - k315_nz
    current_dist = np.linalg.norm(vec)
    print(f"Initial K315-NZ -> Ub-G76-C: {current_dist:.1f} Å")

    vec_unit = vec / current_dist
    for atom in cgas.get_atoms():
        atom.coord += vec_unit * (current_dist - TARGET_DIST)

    # Re-read K315 position
    for res in cgas:
        if res.resname == 'LYS' and res.id[1] == 315 and 'NZ' in res:
            new_dist = np.linalg.norm(res['NZ'].coord - ub_g76_c)
            break
    print(f"After translation: {new_dist:.1f} Å")

    # Clash check + resolve
    n = count_clashes(cgas)
    print(f"CA clashes: {n}")
    if n > 0:
        ok, it = resolve_clashes(cgas)
        print(f"After resolution ({it} iter): {count_clashes(cgas)} clashes")
        for res in cgas:
            if res.resname == 'LYS' and res.id[1] == 315 and 'NZ' in res:
                final_dist = np.linalg.norm(res['NZ'].coord - ub_g76_c)
                break
        print(f"K315-NZ -> Ub-G76-C after resolve: {final_dist:.1f} Å")

    # Merge
    e2_c = e2.copy()
    ub_c = ub.copy()
    for res in e2_c:
        res.parent.id = 'E'
    for res in ub_c:
        res.parent.id = 'U'

    s = Structure.Structure(f'minimal_quaternary_{label}')
    m = Model.Model(0)
    s.add(m)
    for chain_obj in [e2_c, ub_c, cgas]:
        m.add(chain_obj)

    io = PDBIO()
    io.set_structure(s)
    out_pdb = OUTDIR / f'minimal_quaternary_{label}_raw.pdb'
    io.save(str(out_pdb))
    print(f"Saved raw: {out_pdb}")

print("\nDONE. Next: pdb4amber on raw PDBs, then tleap, min, MD")
