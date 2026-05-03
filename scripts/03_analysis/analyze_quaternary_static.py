#!/usr/bin/env python3
"""
Static structure analysis for quaternary MVP complexes (WT and 4mut).

Uses chain-ID-bearing PDB files:
  - WT:  quaternary_mvp_raw.pdb
  - 4mut: quaternary_mvp_4mut_shifted.pdb

Metrics:
  1. Catalytic geometry
  2. Interface contacts
  3. Overall assembly quality
  4. Comparison WT vs 4mut
"""

import numpy as np
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis import distances

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/structures/quaternary_mvp"

def analyze_structure(pdb_path, label):
    """Analyze a single quaternary structure."""
    print(f"\n{'='*60}")
    print(f"Analyzing {label}: {pdb_path.name}")
    print(f"{'='*60}")

    u = mda.Universe(str(pdb_path))

    # =========================================================================
    # 1. CATALYTIC GEOMETRY
    # =========================================================================
    print("\n--- 1. Catalytic Geometry ---")

    # Key atoms
    k315_nz = u.select_atoms("chainID C and resname LYS and resid 315 and name NZ")
    ub_g76_c = u.select_atoms("chainID U and resname GLY and resid 76 and name C")
    e2_k85_nz = u.select_atoms("chainID E and resname LYS and resid 85 and name NZ")
    e2_k85_ca = u.select_atoms("chainID E and resname LYS and resid 85 and name CA")

    if len(k315_nz) == 0:
        print("WARNING: K315 NZ not found!")
        k315_pos = None
    else:
        k315_pos = k315_nz.positions[0]
        print(f"K315 NZ: {k315_pos}")

    if len(ub_g76_c) == 0:
        print("WARNING: Ub G76 C not found!")
        ub_pos = None
    else:
        ub_pos = ub_g76_c.positions[0]
        print(f"Ub G76 C: {ub_pos}")

    if len(e2_k85_nz) == 0:
        print("WARNING: E2 K85 NZ not found!")
        e2_k85_pos = None
    else:
        e2_k85_pos = e2_k85_nz.positions[0]
        print(f"E2 K85 NZ: {e2_k85_pos}")

    # Distances
    if k315_pos is not None and ub_pos is not None:
        d_k315_ub = np.linalg.norm(k315_pos - ub_pos)
        print(f"K315 NZ -> Ub G76 C: {d_k315_ub:.2f} Å")
    else:
        d_k315_ub = None

    if k315_pos is not None and e2_k85_pos is not None:
        d_k315_e2 = np.linalg.norm(k315_pos - e2_k85_pos)
        print(f"K315 NZ -> E2 K85 NZ: {d_k315_e2:.2f} Å")
    else:
        d_k315_e2 = None

    if e2_k85_pos is not None and ub_pos is not None:
        d_e2_ub = np.linalg.norm(e2_k85_pos - ub_pos)
        print(f"E2 K85 NZ -> Ub G76 C: {d_e2_ub:.2f} Å")
    else:
        d_e2_ub = None

    # K315 side chain orientation
    k315_ca = u.select_atoms("chainID C and resname LYS and resid 315 and name CA")
    if len(k315_ca) > 0 and k315_pos is not None and ub_pos is not None:
        sc_vec = k315_pos - k315_ca.positions[0]
        sc_vec /= np.linalg.norm(sc_vec)
        target_vec = ub_pos - k315_ca.positions[0]
        target_vec /= np.linalg.norm(target_vec)
        angle = np.degrees(np.arccos(np.clip(np.dot(sc_vec, target_vec), -1.0, 1.0)))
        print(f"K315 side chain angle toward Ub G76: {angle:.1f}°")

    # =========================================================================
    # 2. INTERFACE CONTACTS
    # =========================================================================
    print("\n--- 2. Interface Contacts ---")

    chains = {
        'RING1': u.select_atoms("chainID R"),
        'RING2': u.select_atoms("chainID S"),
        'E2': u.select_atoms("chainID E"),
        'Ub': u.select_atoms("chainID U"),
        'SPRY': u.select_atoms("chainID P"),
        'cGAS': u.select_atoms("chainID C"),
    }

    # Report chain sizes
    for name, sel in chains.items():
        n_res = len(sel.residues)
        n_ca = len(sel.select_atoms("name CA"))
        print(f"{name}: {n_res} residues, {n_ca} CA atoms")

    # Inter-chain contacts (CA-CA < 8 Å)
    interface_pairs = [
        ('cGAS', 'SPRY'),
        ('cGAS', 'RING1'),
        ('cGAS', 'RING2'),
        ('cGAS', 'E2'),
        ('cGAS', 'Ub'),
        ('SPRY', 'RING1'),
        ('SPRY', 'RING2'),
        ('RING1', 'E2'),
        ('RING2', 'E2'),
        ('E2', 'Ub'),
    ]

    contact_cutoff = 8.0
    interface_contacts = {}
    for name1, name2 in interface_pairs:
        sel1 = chains[name1].select_atoms("name CA")
        sel2 = chains[name2].select_atoms("name CA")
        if len(sel1) == 0 or len(sel2) == 0:
            continue
        dist_matrix = distances.distance_array(sel1.positions, sel2.positions)
        n_contacts = np.sum(dist_matrix < contact_cutoff)
        interface_contacts[f"{name1}-{name2}"] = n_contacts
        print(f"{name1}-{name2} CA contacts (< {contact_cutoff} Å): {n_contacts}")

    # =========================================================================
    # 3. OVERALL ASSEMBLY QUALITY
    # =========================================================================
    print("\n--- 3. Overall Assembly Quality ---")

    # Radius of gyration for each chain and whole complex
    for name, sel in chains.items():
        if len(sel) > 0:
            rg = sel.radius_of_gyration()
            print(f"{name} Rg: {rg:.2f} Å")

    whole = u.select_atoms("protein")
    if len(whole) > 0:
        rg_total = whole.radius_of_gyration()
        print(f"Total complex Rg: {rg_total:.2f} Å")

    # Center-of-mass distances between chains
    print("\nInter-chain COM distances (Å):")
    coms = {name: sel.center_of_mass() for name, sel in chains.items() if len(sel) > 0}
    for i, (n1, c1) in enumerate(coms.items()):
        for n2, c2 in list(coms.items())[i+1:]:
            d = np.linalg.norm(c1 - c2)
            print(f"  {n1}-{n2}: {d:.1f} Å")

    # =========================================================================
    # 4. E2~Ub CONFORMATION (vs 5FER closed)
    # =========================================================================
    print("\n--- 4. E2~Ub Conformation ---")
    fer_path = BASE / "data/structures/quaternary_mvp/5FER.pdb"
    if fer_path.exists():
        u_fer = mda.Universe(str(fer_path))
        # E2 CA alignment
        e2_sel = u.select_atoms("chainID E and name CA")
        fer_e2 = u_fer.select_atoms("chainID B and name CA")
        if len(e2_sel) > 0 and len(fer_e2) > 0:
            # Simple subset alignment (first min(N) CAs)
            n = min(len(e2_sel), len(fer_e2))
            rmsd_e2 = np.sqrt(np.mean(np.sum((e2_sel.positions[:n] - fer_e2.positions[:n])**2, axis=1)))
            print(f"E2 CA RMSD vs 5FER (first {n} CAs): {rmsd_e2:.2f} Å")

        # Ub CA alignment
        ub_sel = u.select_atoms("chainID U and name CA")
        fer_ub = u_fer.select_atoms("chainID C and name CA")
        if len(ub_sel) > 0 and len(fer_ub) > 0:
            n = min(len(ub_sel), len(fer_ub))
            rmsd_ub = np.sqrt(np.mean(np.sum((ub_sel.positions[:n] - fer_ub.positions[:n])**2, axis=1)))
            print(f"Ub CA RMSD vs 5FER (first {n} CAs): {rmsd_ub:.2f} Å")

    return {
        'label': label,
        'd_k315_ub': d_k315_ub,
        'd_k315_e2': d_k315_e2,
        'd_e2_ub': d_e2_ub,
        'interface_contacts': interface_contacts,
        'rg_total': rg_total if len(whole) > 0 else None,
    }


def main():
    results = []

    # WT
    wt_pdb = OUTDIR / "quaternary_mvp_raw.pdb"
    if wt_pdb.exists():
        results.append(analyze_structure(wt_pdb, "WT"))
    else:
        print(f"WARNING: {wt_pdb} not found")

    # 4mut
    mut_pdb = OUTDIR / "quaternary_mvp_4mut_shifted.pdb"
    if mut_pdb.exists():
        results.append(analyze_structure(mut_pdb, "4mut"))
    else:
        print(f"WARNING: {mut_pdb} not found")

    # Summary comparison
    print(f"\n{'='*60}")
    print("SUMMARY: WT vs 4mut")
    print(f"{'='*60}")

    if len(results) == 2:
        wt = results[0]
        mut = results[1]
        print(f"\nK315 -> Ub G76:")
        print(f"  WT:   {wt['d_k315_ub']:.2f} Å" if wt['d_k315_ub'] else "  WT: N/A")
        print(f"  4mut: {mut['d_k315_ub']:.2f} Å" if mut['d_k315_ub'] else "  4mut: N/A")

        print(f"\nK315 -> E2 K85:")
        print(f"  WT:   {wt['d_k315_e2']:.2f} Å" if wt['d_k315_e2'] else "  WT: N/A")
        print(f"  4mut: {mut['d_k315_e2']:.2f} Å" if mut['d_k315_e2'] else "  4mut: N/A")

        print(f"\nE2 K85 -> Ub G76:")
        print(f"  WT:   {wt['d_e2_ub']:.2f} Å" if wt['d_e2_ub'] else "  WT: N/A")
        print(f"  4mut: {mut['d_e2_ub']:.2f} Å" if mut['d_e2_ub'] else "  4mut: N/A")

        print(f"\nInterface contacts (WT -> 4mut delta):")
        for key in wt['interface_contacts']:
            wt_val = wt['interface_contacts'][key]
            mut_val = mut['interface_contacts'].get(key, 0)
            delta = mut_val - wt_val
            print(f"  {key}: {wt_val} -> {mut_val} (Δ={delta:+d})")

        print(f"\nTotal Rg:")
        print(f"  WT:   {wt['rg_total']:.2f} Å" if wt['rg_total'] else "  WT: N/A")
        print(f"  4mut: {mut['rg_total']:.2f} Å" if mut['rg_total'] else "  4mut: N/A")


if __name__ == "__main__":
    main()
