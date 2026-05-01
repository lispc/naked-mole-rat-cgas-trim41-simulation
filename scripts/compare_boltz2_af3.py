#!/usr/bin/env python3
"""
Compare Boltz-2 truncated prediction with AF3 reference structure.

Analysis pipeline:
1. Load both structures (Boltz-2 CIF, AF3 PDB)
2. Align cGAS CT (chain A), measure TRIM41 relative position
3. Per-residue CA-RMSD after alignment
4. Interface contact map comparison
5. Active site distances (4 mutation sites)
6. pLDDT profile comparison
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms
except ImportError:
    print("MDAnalysis not available. Install with: pip install MDAnalysis")
    sys.exit(1)


def kabsch_align(mobile_coords, ref_coords):
    """Kabsch alignment, returns R and translation."""
    # Center
    mobile_center = mobile_coords.mean(axis=0)
    ref_center = ref_coords.mean(axis=0)
    mobile_c = mobile_coords - mobile_center
    ref_c = ref_coords - ref_center
    # Covariance
    H = mobile_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = ref_center - R @ mobile_center
    return R, t


def align_and_rmsd(mobile, ref, selection="name CA"):
    """Align mobile to ref using selection, return RMSD."""
    m_sel = mobile.select_atoms(selection)
    r_sel = ref.select_atoms(selection)
    if len(m_sel) != len(r_sel):
        return None, None
    R, t = kabsch_align(m_sel.positions, r_sel.positions)
    aligned = (mobile.atoms.positions @ R.T) + t
    rmsd = np.sqrt(np.mean((aligned[m_sel.indices] - r_sel.positions)**2))
    return rmsd, aligned


def compute_per_residue_rmsd(mobile_aligned, ref, chain_sel, cgas_resoffset=0):
    """Compute per-residue CA distance after global alignment."""
    m_ca = mobile_aligned.select_atoms(f"{chain_sel} and name CA")
    r_ca = ref.select_atoms(f"{chain_sel} and name CA")
    if len(m_ca) != len(r_ca):
        return None, None
    distances = np.linalg.norm(m_ca.positions - r_ca.positions, axis=1)
    resids = r_ca.resids + cgas_resoffset
    return resids, distances


def compute_interface_contacts(u, chainA_sel, chainB_sel, cutoff=5.0):
    """Return list of (residA, residB) contact pairs."""
    A = u.select_atoms(chainA_sel)
    B = u.select_atoms(chainB_sel)
    contacts = []
    for a in A.residues:
        for b in B.residues:
            dist = np.min(np.linalg.norm(a.atoms.positions[:, None] - b.atoms.positions[None, :], axis=2))
            if dist < cutoff:
                contacts.append((a.resid, b.resid))
    return set(contacts)


def compute_active_site_distances(u, cgas_sel, trim41_sel, sites):
    """
    sites: dict of {name: resid_in_full_numbering}
    Returns dict of {name: min_distance_to_trim41}
    """
    cgas = u.select_atoms(cgas_sel)
    trim41 = u.select_atoms(trim41_sel)
    results = {}
    for name, resid in sites.items():
        site_atoms = cgas.select_atoms(f"resid {resid}")
        if len(site_atoms) == 0:
            results[name] = None
            continue
        dists = np.linalg.norm(site_atoms.positions[:, None] - trim41.positions[None, :], axis=2)
        results[name] = np.min(dists)
    return results


def parse_plddt_from_cif(cif_path):
    """Extract per-atom pLDDT from Boltz-2 CIF."""
    u = mda.Universe(cif_path)
    ca = u.select_atoms("name CA")
    # pLDDT is stored in B-factor column in mmCIF for Boltz
    plddts = ca.tempfactors if hasattr(ca, 'tempfactors') else np.zeros(len(ca))
    return ca.resids, plddts


def main():
    # Paths
    boltz_cif = Path("data/boltz_test_truncated/boltz2_model.pdb")
    af3_pdb_a = Path("structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb")
    af3_pdb_b = Path("structures/af3_raw/job1_Hsap_WT/trim41_SPRY_413-630.pdb")
    outdir = Path("data/boltz_test_truncated/comparison")
    outdir.mkdir(exist_ok=True)

    print("Loading structures...")
    u_boltz = mda.Universe(str(boltz_cif))
    # AF3: need to load both chains and merge
    u_af3_a = mda.Universe(str(af3_pdb_a))
    u_af3_b = mda.Universe(str(af3_pdb_b))

    # Merge AF3 chains into one universe
    # Need to renumber chains to A/B
    u_af3_a.atoms.segids = 'A'
    u_af3_b.atoms.segids = 'B'
    from MDAnalysis.core.universe import Merge
    u_af3 = Merge(u_af3_a.atoms, u_af3_b.atoms)

    # Rename chains for clarity
    # Boltz-2 CIF: chain A = cGAS, chain B = TRIM41
    # AF3 PDB: chain A = cGAS, chain B = TRIM41 (from ranked_0_chain_A/B)

    print(f"Boltz-2: {len(u_boltz.atoms)} atoms")
    print(f"AF3: {len(u_af3.atoms)} atoms")

    # --- 1. Global alignment on cGAS CT, measure TRIM41 RMSD ---
    print("\n=== 1. Global alignment (cGAS CA) ===")
    boltz_cgas_ca = u_boltz.select_atoms("segid A and name CA")
    af3_cgas_ca = u_af3.select_atoms("segid A and name CA")
    print(f"cGAS CA atoms: Boltz={len(boltz_cgas_ca)}, AF3={len(af3_cgas_ca)}")

    boltz_trim41_ca = u_boltz.select_atoms("segid B and name CA")
    af3_trim41_ca = u_af3.select_atoms("segid B and name CA")
    print(f"TRIM41 CA atoms: Boltz={len(boltz_trim41_ca)}, AF3={len(af3_trim41_ca)}")

    # Align Boltz to AF3 using cGAS CA
    R, t = kabsch_align(boltz_cgas_ca.positions, af3_cgas_ca.positions)
    boltz_aligned_pos = (u_boltz.atoms.positions @ R.T) + t
    u_boltz.atoms.positions = boltz_aligned_pos

    # RMSD after alignment
    cgas_rmsd = np.sqrt(np.mean((boltz_cgas_ca.positions - af3_cgas_ca.positions)**2))
    trim41_rmsd = np.sqrt(np.mean((boltz_trim41_ca.positions - af3_trim41_ca.positions)**2))
    print(f"cGAS CA RMSD (after align): {cgas_rmsd:.2f} Å")
    print(f"TRIM41 CA RMSD (after align): {trim41_rmsd:.2f} Å")

    # --- 2. Per-residue CA distance after alignment ---
    print("\n=== 2. Per-residue CA distance ===")
    cgas_resids, cgas_dists = compute_per_residue_rmsd(
        u_boltz, u_af3, "segid A and name CA", cgas_resoffset=0
    )
    trim41_resids, trim41_dists = compute_per_residue_rmsd(
        u_boltz, u_af3, "segid B and name CA", cgas_resoffset=0
    )

    # --- 3. Interface contact comparison ---
    print("\n=== 3. Interface contacts (< 5 Å) ===")
    boltz_contacts = compute_interface_contacts(u_boltz, "segid A", "segid B", cutoff=5.0)
    af3_contacts = compute_interface_contacts(u_af3, "segid A", "segid B", cutoff=5.0)
    shared = boltz_contacts & af3_contacts
    boltz_only = boltz_contacts - af3_contacts
    af3_only = af3_contacts - boltz_contacts
    print(f"Boltz-2 contacts: {len(boltz_contacts)}")
    print(f"AF3 contacts: {len(af3_contacts)}")
    print(f"Shared contacts: {len(shared)}")
    print(f"Jaccard similarity: {len(shared)/len(boltz_contacts | af3_contacts) if (boltz_contacts | af3_contacts) else 0:.3f}")

    # --- 4. Active site distances ---
    print("\n=== 4. Active site distances (Hsap numbering) ===")
    # Hsap sites: D431, K479, L495, K498
    # AF3 truncated file uses original numbering (200-522)
    sites = {"D431": 431, "K479": 479, "L495": 495, "K498": 498}
    boltz_dists = compute_active_site_distances(u_boltz, "segid A", "segid B", sites)
    af3_dists = compute_active_site_distances(u_af3, "segid A", "segid B", sites)
    for name in sites:
        print(f"  {name}: Boltz={boltz_dists[name]:.1f}Å, AF3={af3_dists[name]:.1f}Å, Δ={boltz_dists[name]-af3_dists[name]:+.1f}Å")

    # --- 5. pLDDT profiles ---
    print("\n=== 5. pLDDT profiles ===")
    try:
        b_resids, b_plddt = parse_plddt_from_cif(str(boltz_cif))
        print(f"Boltz-2 pLDDT: mean={b_plddt.mean():.3f}, min={b_plddt.min():.3f}, max={b_plddt.max():.3f}")
    except Exception as e:
        print(f"Could not parse Boltz-2 pLDDT: {e}")
        b_resids, b_plddt = None, None

    # --- Plotting ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Per-residue distance (cGAS)
    if cgas_resids is not None:
        ax = axes[0, 0]
        ax.plot(cgas_resids, cgas_dists, 'b-', linewidth=0.8, label='cGAS CT')
        ax.axhline(2.0, color='g', linestyle='--', alpha=0.5, label='2 Å threshold')
        ax.axhline(5.0, color='r', linestyle='--', alpha=0.5, label='5 Å threshold')
        ax.set_xlabel('Residue (Hsap numbering)')
        ax.set_ylabel('CA distance (Å)')
        ax.set_title('cGAS CT: Boltz-2 vs AF3 per-residue CA distance')
        ax.legend()

    # Plot 2: Per-residue distance (TRIM41)
    if trim41_resids is not None:
        ax = axes[0, 1]
        ax.plot(trim41_resids, trim41_dists, 'r-', linewidth=0.8, label='TRIM41 SPRY')
        ax.axhline(2.0, color='g', linestyle='--', alpha=0.5)
        ax.axhline(5.0, color='r', linestyle='--', alpha=0.5)
        ax.set_xlabel('Residue (Hsap numbering)')
        ax.set_ylabel('CA distance (Å)')
        ax.set_title('TRIM41 SPRY: Boltz-2 vs AF3 per-residue CA distance')
        ax.legend()

    # Plot 3: Contact overlap
    ax = axes[1, 0]
    categories = ['Shared', 'Boltz-only', 'AF3-only']
    values = [len(shared), len(boltz_only), len(af3_only)]
    colors = ['green', 'blue', 'orange']
    ax.bar(categories, values, color=colors, alpha=0.7)
    ax.set_ylabel('Contact pairs')
    ax.set_title(f'Interface contacts (Jaccard={len(shared)/len(boltz_contacts | af3_contacts):.2f})')

    # Plot 4: Active site distance comparison
    ax = axes[1, 1]
    names = list(sites.keys())
    boltz_vals = [boltz_dists[n] for n in names]
    af3_vals = [af3_dists[n] for n in names]
    x = np.arange(len(names))
    width = 0.35
    ax.bar(x - width/2, boltz_vals, width, label='Boltz-2', color='blue', alpha=0.7)
    ax.bar(x + width/2, af3_vals, width, label='AF3', color='orange', alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('Distance to TRIM41 (Å)')
    ax.set_title('Active site distances')
    ax.legend()

    plt.tight_layout()
    plot_path = outdir / "boltz2_vs_af3_comparison.png"
    plt.savefig(plot_path, dpi=300)
    print(f"\nPlot saved: {plot_path}")

    # Save numerical results
    summary = {
        "cgas_ca_rmsd": float(cgas_rmsd),
        "trim41_ca_rmsd": float(trim41_rmsd),
        "boltz_contacts": len(boltz_contacts),
        "af3_contacts": len(af3_contacts),
        "shared_contacts": len(shared),
        "jaccard_similarity": float(len(shared)/len(boltz_contacts | af3_contacts)) if (boltz_contacts | af3_contacts) else 0.0,
        "active_sites": {name: {"boltz": float(boltz_dists[name]), "af3": float(af3_dists[name])} for name in sites},
    }
    import json
    with open(outdir / "comparison_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Summary saved: {outdir / 'comparison_summary.json'}")

    print("\n=== Overall Assessment ===")
    print(f"TRIM41 relative RMSD: {trim41_rmsd:.2f} Å")
    if trim41_rmsd < 5.0:
        print("→ TRIM41 placement is VERY similar between Boltz-2 and AF3")
    elif trim41_rmsd < 10.0:
        print("→ TRIM41 placement is MODERATELY similar")
    else:
        print("→ TRIM41 placement is DIVERGENT")

    print(f"Interface contact Jaccard: {summary['jaccard_similarity']:.2f}")
    if summary['jaccard_similarity'] > 0.7:
        print("→ Interface contact network is HIGHLY conserved")
    elif summary['jaccard_similarity'] > 0.4:
        print("→ Interface contact network is PARTIALLY conserved")
    else:
        print("→ Interface contact network is DIFFERENT")


if __name__ == "__main__":
    main()
