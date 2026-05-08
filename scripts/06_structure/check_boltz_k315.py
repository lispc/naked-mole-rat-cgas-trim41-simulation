#!/usr/bin/env python3
"""Check K315 position in Boltz-2 cGAS+DNA+TRIM41 predictions (models 0-4)."""
import gemmi
import numpy as np
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent.parent
CIF_DIR = BASE / "data/boltz_cgas_dna_trim41/boltz_results_boltz_cgas_dna_trim41/predictions/boltz_cgas_dna_trim41"

for model_idx in range(5):
    cif_path = CIF_DIR / f"boltz_cgas_dna_trim41_model_{model_idx}.cif"
    st = gemmi.read_structure(str(cif_path))
    ch_a = st[0]["A"]   # cGAS
    ch_b = st[0]["B"]   # TRIM41 SPRY

    # 1. SPRY heavy atoms for distance calculations
    spry_heavy = []
    spry_ca = []
    for res in ch_b:
        for atom in res:
            if atom.name == "CA":
                spry_ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
            if atom.name not in ("H", "HA"):
                spry_heavy.append([atom.pos.x, atom.pos.y, atom.pos.z])
    spry_com = np.array(spry_ca).mean(axis=0)
    spry_heavy = np.array(spry_heavy)

    # 2. All LYS NZ in cGAS
    lys_info = []
    for i, res in enumerate(ch_a):
        if res.name == "LYS":
            for atom in res:
                if atom.name == "NZ":
                    pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                    dist_com = np.linalg.norm(pos - spry_com)
                    dist_min = np.min(np.linalg.norm(spry_heavy - pos, axis=1))
                    lys_info.append((i + 1, dist_com, dist_min, pos))

    # 3. SPRY-cGAS interface (<5A)
    close_cgas = []
    for j, res in enumerate(ch_a):
        min_d = 999.0
        for atom in res:
            if atom.name not in ("H", "HA", "CA"):
                pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                min_d = min(min_d, np.min(np.linalg.norm(spry_heavy - pos, axis=1)))
        if min_d < 5.0:
            close_cgas.append((j + 1, res.name, min_d))

    # 4. K315 = CIF pos 225 (SASKMLSKFRK sequence context)
    k315_dist_com = None
    k315_dist_min = None
    for cpos, dcom, dmin, _ in lys_info:
        if cpos == 225:
            k315_dist_com = dcom
            k315_dist_min = dmin

    lys_info.sort(key=lambda x: x[1])

    print(f"--- Model {model_idx} ---")
    print(f"  SPRY COM: ({spry_com[0]:.1f}, {spry_com[1]:.1f}, {spry_com[2]:.1f})")
    print(f"  K315 (CIF pos 225): COM dist={k315_dist_com:.1f}A, min SPRY dist={k315_dist_min:.1f}A")
    print(f"  Closest LYS to SPRY COM: pos {lys_info[0][0]}, {lys_info[0][1]:.1f}A COM / {lys_info[0][2]:.1f}A min")
    print(f"  2nd closest: pos {lys_info[1][0]}, {lys_info[1][1]:.1f}A COM / {lys_info[1][2]:.1f}A min")

    # Interface residues near K315
    k315_nearby = [(p, n, d) for p, n, d in close_cgas if abs(p - 225) < 25]
    if k315_nearby:
        print(f"  SPRY-contacting cGAS res near K315 (±25 pos): {[(p, n) for p, n, _ in k315_nearby[:5]]}...")
    else:
        print(f"  NO SPRY-contacting cGAS residues within 25 positions of K315")

    # LYS <30A from SPRY COM
    near_lys = [(p, dcom, dmin) for p, dcom, dmin, _ in lys_info if dcom < 30]
    print(f"  LYS <30A from SPRY COM: {[(p, f'{dmin:.1f}') for p, _, dmin in near_lys]}")

    print()
