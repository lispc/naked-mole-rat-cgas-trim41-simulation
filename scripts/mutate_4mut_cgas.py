#!/usr/bin/env python3
"""Mutate cGAS WT PDB to 4mut (D431S, K479E, L495Y, K498T). Places new side chain atoms with idealized geometry."""

import numpy as np
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")


def norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else v


def perp(v):
    arb = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(v, arb)) > 0.99:
        arb = np.array([0.0, 1.0, 0.0])
    return norm(np.cross(v, arb))


def parse_pdb(path):
    """Parse PDB into list of dicts. Returns atoms list."""
    atoms = []
    with open(path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atoms.append({
                'serial': int(line[6:11]),
                'name': line[12:16].strip(),
                'resName': line[17:20].strip(),
                'chainID': line[21].strip(),
                'resSeq': int(line[22:26]),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
                'elem': line[76:78].strip() if len(line) > 76 else '',
            })
    return atoms


def write_pdb(atoms, path):
    """Write atoms list to PDB file with strict PDB column format."""
    with open(path, 'w') as f:
        for i, a in enumerate(atoms):
            serial = i + 1
            elem = a.get('elem', a['name'][0]) if a.get('elem') else a['name'][0]
            # Atom name: cols 13-16 (4 chars). For 1-char element: " EEE" left-justified
            if len(a['name']) <= 3:
                name_field = f" {a['name']:<3s}"  # " N  ", " CA "
            else:
                name_field = f"{a['name']:>4s}"
            # PDB columns: 1-6 record, 7-11 serial, 12 blank, 13-16 atomName,
            # 17 altLoc, 18-20 resName, 21 blank, 22 chainID, 23-26 resSeq
            f.write(f"ATOM  {serial:5d} {name_field} {a['resName']:>3s} {a['chainID']}{a['resSeq']:>4d}    "
                    f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{1.0:6.2f}{0.0:6.2f}          {elem:>2s}\n")
        f.write("END\n")


def get_pos(atoms, chain, res, name):
    for a in atoms:
        if a['chainID'] == chain and a['resSeq'] == res and a['name'] == name:
            return np.array([a['x'], a['y'], a['z']])
    return None


def mutate(input_pdb, output_pdb):
    atoms = parse_pdb(input_pdb)

    # Track removal indices (process in reverse to maintain order)
    to_remove = set()
    to_add = []

    # ---- D431S: ASP → SER (remove CG, OD1, OD2; add OG) ----
    print("D431S: ASP → SER")
    ca = get_pos(atoms, 'A', 431, 'CA')
    cb = get_pos(atoms, 'A', 431, 'CB')
    for i, a in enumerate(atoms):
        if a['chainID'] == 'A' and a['resSeq'] == 431:
            a['resName'] = 'SER'
            if a['name'] in ('CG', 'OD1', 'OD2'):
                to_remove.add(i)
    # Add OG at tetrahedral position from CB
    dir_ca_cb = norm(cb - ca)
    p = perp(dir_ca_cb)
    rad = np.deg2rad(180 - 109.5)
    og = cb + 1.43 * (dir_ca_cb * np.cos(rad) + p * np.sin(rad))
    to_add.append({'name': 'OG', 'resName': 'SER', 'chainID': 'A', 'resSeq': 431,
                   'x': og[0], 'y': og[1], 'z': og[2], 'elem': 'O'})

    # ---- K479E: LYS → GLU (remove CE,NZ; keep CG, rename CD, add OE1,OE2) ----
    print("K479E: LYS → GLU")
    cg = get_pos(atoms, 'A', 479, 'CG')
    cd_lys = get_pos(atoms, 'A', 479, 'CD')
    for i, a in enumerate(atoms):
        if a['chainID'] == 'A' and a['resSeq'] == 479:
            a['resName'] = 'GLU'
            if a['name'] in ('CE', 'NZ'):
                to_remove.add(i)
    # Place CD at standard geometry from CG (staggered relative to LYS CD)
    dir_cg_cd = norm(cd_lys - cg)
    p2 = perp(dir_cg_cd)
    rad2 = np.deg2rad(180 - 109.5)
    cd_glu = cg + 1.52 * (dir_cg_cd * np.cos(rad2) + p2 * np.sin(rad2))
    # OE1, OE2 at 120 deg from CG-CD, in plane perpendicular to CG-CD
    dir_cg_cd2 = norm(cd_glu - cg)
    p3 = perp(dir_cg_cd2)
    # Planar carboxylate
    oe1 = cd_glu + 1.25 * (dir_cg_cd2 * np.cos(np.deg2rad(60)) + p3 * np.sin(np.deg2rad(60)))
    oe2 = cd_glu + 1.25 * (dir_cg_cd2 * np.cos(np.deg2rad(60)) - p3 * np.sin(np.deg2rad(60)))
    # Update CD position
    for i, a in enumerate(atoms):
        if a['chainID'] == 'A' and a['resSeq'] == 479 and a['name'] == 'CD':
            a['x'], a['y'], a['z'] = cd_glu[0], cd_glu[1], cd_glu[2]
    to_add.append({'name': 'OE1', 'resName': 'GLU', 'chainID': 'A', 'resSeq': 479,
                   'x': oe1[0], 'y': oe1[1], 'z': oe1[2], 'elem': 'O'})
    to_add.append({'name': 'OE2', 'resName': 'GLU', 'chainID': 'A', 'resSeq': 479,
                   'x': oe2[0], 'y': oe2[1], 'z': oe2[2], 'elem': 'O'})

    # ---- L495Y: LEU → TYR (keep CG, CD1, CD2; add CE1, CE2, CZ, OH) ----
    print("L495Y: LEU → TYR")
    cd1 = get_pos(atoms, 'A', 495, 'CD1')
    cd2 = get_pos(atoms, 'A', 495, 'CD2')
    for i, a in enumerate(atoms):
        if a['chainID'] == 'A' and a['resSeq'] == 495:
            a['resName'] = 'TYR'
    # LEU: CB-CG-CD1/CD2. TYR: CB-CG-CD1-CD2-CE1-CE2-CZ-OH
    # Build phenyl ring: CE1 beyond CD1, CE2 beyond CD2, CZ on opposite side
    # CD1 and CD2 are the ring's first two carbons
    ring_vec = norm(cd2 - cd1)  # along CD1-CD2 edge
    ring_normal = norm(np.cross(cd1 - get_pos(atoms, 'A', 495, 'CG'), ring_vec))
    # CE1: extend from CD1 opposite to CD2 direction (in ring plane)
    ce1 = cd1 + 1.39 * norm(cd1 - cd2 + ring_normal * 0.5)
    ce2 = cd2 + 1.39 * norm(cd2 - cd1 + ring_normal * 0.5)
    # CZ: opposite side of ring from CD1-CD2
    cz = (ce1 + ce2) / 2 + ring_normal * 1.39
    oh = cz + norm(ring_normal) * 1.36

    to_add.append({'name': 'CE1', 'resName': 'TYR', 'chainID': 'A', 'resSeq': 495,
                   'x': ce1[0], 'y': ce1[1], 'z': ce1[2], 'elem': 'C'})
    to_add.append({'name': 'CE2', 'resName': 'TYR', 'chainID': 'A', 'resSeq': 495,
                   'x': ce2[0], 'y': ce2[1], 'z': ce2[2], 'elem': 'C'})
    to_add.append({'name': 'CZ', 'resName': 'TYR', 'chainID': 'A', 'resSeq': 495,
                   'x': cz[0], 'y': cz[1], 'z': cz[2], 'elem': 'C'})
    to_add.append({'name': 'OH', 'resName': 'TYR', 'chainID': 'A', 'resSeq': 495,
                   'x': oh[0], 'y': oh[1], 'z': oh[2], 'elem': 'O'})

    # ---- K498T: LYS → THR (remove CD,CE,NZ; add OG1,CG2) ----
    print("K498T: LYS → THR")
    ca498 = get_pos(atoms, 'A', 498, 'CA')
    cb498 = get_pos(atoms, 'A', 498, 'CB')
    for i, a in enumerate(atoms):
        if a['chainID'] == 'A' and a['resSeq'] == 498:
            a['resName'] = 'THR'
            if a['name'] in ('CG', 'CD', 'CE', 'NZ'):
                to_remove.add(i)
    # Place OG1 and CG2 from CB
    dir_ca_cb498 = norm(cb498 - ca498)
    p4 = perp(dir_ca_cb498)
    rad_109 = np.deg2rad(180 - 109.5)
    og1 = cb498 + 1.43 * (dir_ca_cb498 * np.cos(rad_109) + p4 * np.sin(rad_109))
    cg2 = cb498 + 1.53 * (dir_ca_cb498 * np.cos(rad_109) - p4 * np.sin(rad_109))
    to_add.append({'name': 'OG1', 'resName': 'THR', 'chainID': 'A', 'resSeq': 498,
                   'x': og1[0], 'y': og1[1], 'z': og1[2], 'elem': 'O'})
    to_add.append({'name': 'CG2', 'resName': 'THR', 'chainID': 'A', 'resSeq': 498,
                   'x': cg2[0], 'y': cg2[1], 'z': cg2[2], 'elem': 'C'})

    # ---- Rebuild atom list (remove deleted, insert new) ----
    new_atoms = [a for i, a in enumerate(atoms) if i not in to_remove]

    # Insert new atoms after the CB of their respective residues
    insert_points = {}
    for idx, a in enumerate(new_atoms):
        if a['chainID'] == 'A' and a['name'] == 'CB' and a['resSeq'] in (431, 479, 495, 498):
            insert_points.setdefault(a['resSeq'], []).append(idx + 1)

    # Build final list with insertions (process in reverse to maintain indices)
    all_new = list(to_add)
    # Sort additions: group by residue, order within each residue
    res_order = {431: ['OG'], 479: ['OE1', 'OE2'], 495: ['CE1', 'CE2', 'CZ', 'OH'], 498: ['OG1', 'CG2']}
    additions_by_res = {r: [] for r in res_order}
    for add in all_new:
        additions_by_res.setdefault(add['resSeq'], []).append(add)

    for res_seq, add_list in sorted(additions_by_res.items(), reverse=True):
        if res_seq not in insert_points:
            continue
        insert_idx = max(insert_points[res_seq])  # insert after last CB
        for add in reversed(sorted(add_list, key=lambda x: res_order.get(res_seq, []).index(x['name']) if x['name'] in res_order.get(res_seq, []) else 99)):
            new_atoms.insert(insert_idx, add)

    write_pdb(new_atoms, output_pdb)
    print(f"\nSaved mutated cGAS to: {output_pdb}")
    print(f"  Total atoms: {len(new_atoms)}")
    print(f"  Removed: {len(to_remove)}, Added: {len(to_add)}")


if __name__ == '__main__':
    import sys
    infile = sys.argv[1] if len(sys.argv) > 1 else str(
        BASE / 'structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb')
    outfile = sys.argv[2] if len(sys.argv) > 2 else str(
        BASE / 'data/structures/quaternary_minimal/cgas_4mut.pdb')
    mutate(infile, outfile)
