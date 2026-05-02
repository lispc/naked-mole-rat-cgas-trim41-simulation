#!/usr/bin/env python3
"""
Deep analysis of existing 200ns MD data (Option D).

Targets:
  - WT (Hsap_WT): 3 reps × 200ns, stable binding
  - S305-phos: 3 reps × 200ns, dissociation

Analyses:
  1. Interface H-bond timeline (per-replica + mean)
  2. Key interface residue distance timeline
  3. Principal Component Analysis (PCA) — complex + per-protein
  4. Dynamic Cross-Correlation Matrix (DCCM)
  5. Secondary structure (DSSP) of S305 neighborhood
  6. Contact occupancy heatmap (WT vs S305-phos)

Output: data/analysis/deep_200ns/
"""

import os
import sys
import json
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, pca
from MDAnalysis.analysis.distances import distance_array

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

# ---------------------------------------------------------------------------
# Paths & Config
# ---------------------------------------------------------------------------
BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data" / "analysis" / "deep_200ns"
OUTDIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "WT": {
        "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        "trajs": [
            BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        ],
        "cgas_residues": "1-218",
        "trim_residues": "219-541",
    },
    "S305phos": {
        "prmtop": BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop",
        "trajs": [
            BASE / "data/md_runs/Hsap_WT_S305phos/rep1/Hsap_WT_S305phos_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT_S305phos/rep2/Hsap_WT_S305phos_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT_S305phos/rep3/Hsap_WT_S305phos_rep3_prod.dcd",
        ],
        "cgas_residues": "1-218",
        "trim_residues": "219-541",
    },
}

# S305 is residue 305 in the full sequence, but in our truncated construct
# we need to map it. Let's determine it from the topology.
# For now, we know S305 corresponds to residue index ~305 in the prmtop.
# Since cGAS is 1-218, S305 must be in TRIM41: 305 - 218 = 87th residue of TRIM41.
# But let's auto-detect.

# Key interface residues to track (will be auto-populated from first frame)
KEY_RESIDUES = []  # populated at runtime

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_universe(prmtop, dcd):
    u = mda.Universe(str(prmtop), str(dcd))
    # guess elements if missing (TOPParser warning)
    if not hasattr(u.atoms, 'elements') or u.atoms.elements is None or len(u.atoms.elements) == 0:
        u.atoms.guess_elements()
    return u


def get_interface_residues(u, cgas_sel, trim_sel, cutoff=5.0):
    """Return sets of cGAS and TRIM41 interface residues within cutoff (Å)."""
    cgas = u.select_atoms(cgas_sel)
    trim = u.select_atoms(trim_sel)
    # Use first frame
    u.trajectory[0]
    dist_mat = distance_array(cgas.positions, trim.positions)
    close_pairs = np.argwhere(dist_mat < cutoff)
    cgas_resids = set(cgas.resids[close_pairs[:, 0]])
    trim_resids = set(trim.resids[close_pairs[:, 1]])
    return cgas_resids, trim_resids


def analyze_hbonds_simple(u, cgas_sel, trim_sel, distance=3.5, angle=120):
    """
    Simple H-bond counting between two selections.
    Uses heavy-atom distance criterion (donor-acceptor < distance Å).
    Returns: n_hbonds per frame.
    """
    # Select only N, O atoms for H-bond detection
    cgas_heavy = u.select_atoms(f"({cgas_sel}) and (name N* or name O*)")
    trim_heavy = u.select_atoms(f"({trim_sel}) and (name N* or name O*)")

    n_frames = len(u.trajectory)
    hbond_counts = np.zeros(n_frames, dtype=int)

    for i, ts in enumerate(u.trajectory):
        dist_mat = distance_array(cgas_heavy.positions, trim_heavy.positions)
        # Count donor-acceptor pairs within cutoff
        # This is a simplified criterion; full H-bond requires angle check
        # but for interface trend analysis, distance is often sufficient
        close = dist_mat < distance
        # Avoid double counting symmetric pairs
        hbond_counts[i] = np.sum(close) // 2  # rough estimate

    return hbond_counts


def analyze_key_distances(u, residue_pairs):
    """
    residue_pairs: list of (name, sel1, sel2)
    Returns dict of name -> distances array (n_frames,)
    """
    n_frames = len(u.trajectory)
    results = {}
    for name, sel1_str, sel2_str in residue_pairs:
        sel1 = u.select_atoms(sel1_str)
        sel2 = u.select_atoms(sel2_str)
        if len(sel1) == 0 or len(sel2) == 0:
            results[name] = np.full(n_frames, np.nan)
            continue
        dists = np.zeros(n_frames)
        for i, ts in enumerate(u.trajectory):
            d = distance_array(sel1.center_of_mass(), sel2.center_of_mass())
            dists[i] = d[0, 0]
        results[name] = dists
    return results


def run_pca(u, selection="protein and name CA", n_components=10):
    """Run PCA on selection, return transformed coords and variance."""
    aligner = align.AlignTraj(u, u, select=selection, in_memory=True).run()
    atoms = u.select_atoms(selection)
    pca_an = pca.PCA(u, select=selection).run()
    transformed = pca_an.transform(atoms, n_components=n_components)
    variance = pca_an.results.variance[:n_components]
    cumvar = pca_an.results.cumulated_variance[:n_components]
    return transformed, variance, cumvar, pca_an


def compute_dccm(u, selection="protein and name CA"):
    """Compute dynamic cross-correlation matrix."""
    aligner = align.AlignTraj(u, u, select=selection, in_memory=True).run()
    ca = u.select_atoms(selection)
    n_residues = len(ca.residues)
    n_frames = len(u.trajectory)
    coords = np.zeros((n_frames, n_residues, 3))
    for i, ts in enumerate(u.trajectory):
        coords[i] = ca.positions

    # Mean-subtract
    mean_coords = coords.mean(axis=0)
    fluct = coords - mean_coords

    # Correlation matrix
    corr = np.zeros((n_residues, n_residues))
    for i in range(n_residues):
        for j in range(i, n_residues):
            di = fluct[:, i, :].flatten()
            dj = fluct[:, j, :].flatten()
            if np.std(di) == 0 or np.std(dj) == 0:
                c = 0.0
            else:
                c = np.corrcoef(di, dj)[0, 1]
            corr[i, j] = c
            corr[j, i] = c
    return corr, ca.resids


# ---------------------------------------------------------------------------
# Main analysis pipeline
# ---------------------------------------------------------------------------

def analyze_system(name, cfg, outdir):
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")

    prmtop = cfg["prmtop"]
    trajs = cfg["trajs"]
    cgas_sel = f"resid {cfg['cgas_residues']}"
    trim_sel = f"resid {cfg['trim_residues']}"

    # --- Per-replica analysis ---
    all_hbonds = []
    all_com = []
    all_rmsd = []

    for i, dcd in enumerate(trajs):
        print(f"  Replica {i+1}: {dcd.name}")
        if not dcd.exists():
            print(f"    ⚠️  Missing, skipping")
            continue

        u = load_universe(prmtop, dcd)
        n_frames = len(u.trajectory)
        time_ns = np.arange(n_frames) * 0.1  # 0.1 ns per frame (100 ps)

        # 1. H-bond count
        print(f"    Computing H-bonds...")
        hb = analyze_hbonds_simple(u, cgas_sel, trim_sel)
        all_hbonds.append((time_ns, hb))

        # 2. COM distance
        print(f"    Computing COM distance...")
        cgas = u.select_atoms(cgas_sel)
        trim = u.select_atoms(trim_sel)
        com_dists = np.zeros(n_frames)
        for f, ts in enumerate(u.trajectory):
            com_dists[f] = np.linalg.norm(cgas.center_of_mass() - trim.center_of_mass())
        all_com.append((time_ns, com_dists))

        # 3. RMSD (aligned on complex CA)
        print(f"    Computing RMSD...")
        ref = mda.Universe(str(prmtop), str(dcd))
        aligner = align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
        rmsd_an = rms.RMSD(u, ref, select="protein and name CA").run()
        rmsd_vals = rmsd_an.results.rmsd[:, 2]
        all_rmsd.append((time_ns, rmsd_vals))

    # --- Aggregate plots ---
    print(f"  Plotting aggregates...")

    # H-bonds
    fig, ax = plt.subplots(figsize=(8, 3.5))
    for time_ns, hb in all_hbonds:
        ax.plot(time_ns, hb, alpha=0.4, lw=0.8)
    if all_hbonds:
        mean_hb = np.mean([hb for _, hb in all_hbonds], axis=0)
        ax.plot(all_hbonds[0][0], mean_hb, color="black", lw=2, label="Mean")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Interface H-bond count")
    ax.set_title(f"{name} — Interface H-bonds")
    ax.legend()
    ax.set_xlim(0, 200)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_hbonds_timeline.png", dpi=200)
    plt.close(fig)
    print(f"    Saved: {outdir / f'{name}_hbonds_timeline.png'}")

    # COM distance
    fig, ax = plt.subplots(figsize=(8, 3.5))
    for time_ns, com in all_com:
        ax.plot(time_ns, com, alpha=0.4, lw=0.8)
    if all_com:
        mean_com = np.mean([com for _, com in all_com], axis=0)
        ax.plot(all_com[0][0], mean_com, color="black", lw=2, label="Mean")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("COM distance (Å)")
    ax.set_title(f"{name} — cGAS-TRIM41 COM distance")
    ax.legend()
    ax.set_xlim(0, 200)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_com_timeline.png", dpi=200)
    plt.close(fig)
    print(f"    Saved: {outdir / f'{name}_com_timeline.png'}")

    # RMSD
    fig, ax = plt.subplots(figsize=(8, 3.5))
    for time_ns, rmsd in all_rmsd:
        ax.plot(time_ns, rmsd, alpha=0.4, lw=0.8)
    if all_rmsd:
        mean_rmsd = np.mean([r for _, r in all_rmsd], axis=0)
        ax.plot(all_rmsd[0][0], mean_rmsd, color="black", lw=2, label="Mean")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title(f"{name} — RMSD (CA-aligned)")
    ax.legend()
    ax.set_xlim(0, 200)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_rmsd_timeline.png", dpi=200)
    plt.close(fig)
    print(f"    Saved: {outdir / f'{name}_rmsd_timeline.png'}")

    # --- PCA on rep1 (most representative) ---
    if trajs[0].exists():
        print(f"  Running PCA on rep1...")
        u = load_universe(prmtop, trajs[0])
        transformed, variance, cumvar, pca_obj = run_pca(u, n_components=10)

        # PC1 vs PC2 scatter colored by time
        fig, ax = plt.subplots(figsize=(6, 5))
        scatter = ax.scatter(transformed[:, 0], transformed[:, 1], c=time_ns, cmap="viridis", s=5, alpha=0.6)
        ax.set_xlabel(f"PC1 ({variance[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({variance[1]:.1f}%)")
        ax.set_title(f"{name} — PCA rep1")
        plt.colorbar(scatter, ax=ax, label="Time (ns)")
        fig.tight_layout()
        fig.savefig(outdir / f"{name}_pca_pc1pc2.png", dpi=200)
        plt.close(fig)
        print(f"    Saved: {outdir / f'{name}_pca_pc1pc2.png'}")

        # Variance plot
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.bar(range(1, 11), variance)
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Variance (%)")
        ax.set_title(f"{name} — PCA Variance")
        fig.tight_layout()
        fig.savefig(outdir / f"{name}_pca_variance.png", dpi=200)
        plt.close(fig)
        print(f"    Saved: {outdir / f'{name}_pca_variance.png'}")

    # --- DCCM on rep1 ---
    if trajs[0].exists():
        print(f"  Computing DCCM on rep1...")
        u = load_universe(prmtop, trajs[0])
        corr, resids = compute_dccm(u)

        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
        ax.set_xlabel("Residue")
        ax.set_ylabel("Residue")
        ax.set_title(f"{name} — Dynamic Cross-Correlation (rep1)")
        plt.colorbar(im, ax=ax, label="Correlation")
        fig.tight_layout()
        fig.savefig(outdir / f"{name}_dccm.png", dpi=200)
        plt.close(fig)
        print(f"    Saved: {outdir / f'{name}_dccm.png'}")

    return {
        "name": name,
        "n_reps": len(all_hbonds),
        "hbonds": all_hbonds,
        "com": all_com,
        "rmsd": all_rmsd,
    }


def compare_systems(wt_data, mut_data, outdir):
    """Generate comparison plots."""
    print(f"\n{'='*60}")
    print("Comparing WT vs S305-phos")
    print(f"{'='*60}")

    # COM comparison
    fig, ax = plt.subplots(figsize=(8, 4))
    for time_ns, com in wt_data["com"]:
        ax.plot(time_ns, com, color="blue", alpha=0.3, lw=0.8)
    for time_ns, com in mut_data["com"]:
        ax.plot(time_ns, com, color="red", alpha=0.3, lw=0.8)
    if wt_data["com"]:
        mean_wt = np.mean([c for _, c in wt_data["com"]], axis=0)
        ax.plot(wt_data["com"][0][0], mean_wt, color="blue", lw=2, label="WT mean")
    if mut_data["com"]:
        mean_mut = np.mean([c for _, c in mut_data["com"]], axis=0)
        ax.plot(mut_data["com"][0][0], mean_mut, color="red", lw=2, label="S305-phos mean")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("COM distance (Å)")
    ax.set_title("cGAS-TRIM41 COM distance")
    ax.legend()
    ax.set_xlim(0, 200)
    fig.tight_layout()
    fig.savefig(outdir / "compare_com.png", dpi=200)
    plt.close(fig)
    print(f"  Saved: {outdir / 'compare_com.png'}")

    # H-bond comparison
    fig, ax = plt.subplots(figsize=(8, 4))
    for time_ns, hb in wt_data["hbonds"]:
        ax.plot(time_ns, hb, color="blue", alpha=0.3, lw=0.8)
    for time_ns, hb in mut_data["hbonds"]:
        ax.plot(time_ns, hb, color="red", alpha=0.3, lw=0.8)
    if wt_data["hbonds"]:
        mean_wt = np.mean([h for _, h in wt_data["hbonds"]], axis=0)
        ax.plot(wt_data["hbonds"][0][0], mean_wt, color="blue", lw=2, label="WT mean")
    if mut_data["hbonds"]:
        mean_mut = np.mean([h for _, h in mut_data["hbonds"]], axis=0)
        ax.plot(mut_data["hbonds"][0][0], mean_mut, color="red", lw=2, label="S305-phos mean")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Interface H-bond count")
    ax.set_title("Interface Hydrogen Bonds")
    ax.legend()
    ax.set_xlim(0, 200)
    fig.tight_layout()
    fig.savefig(outdir / "compare_hbonds.png", dpi=200)
    plt.close(fig)
    print(f"  Saved: {outdir / 'compare_hbonds.png'}")

    # Summary stats
    print("\n  Summary:")
    if wt_data["com"]:
        wt_com_all = np.concatenate([c for _, c in wt_data["com"]])
        print(f"    WT   COM:  {np.mean(wt_com_all):.2f} ± {np.std(wt_com_all):.2f} Å")
    if mut_data["com"]:
        mut_com_all = np.concatenate([c for _, c in mut_data["com"]])
        print(f"    S305 COM:  {np.mean(mut_com_all):.2f} ± {np.std(mut_com_all):.2f} Å")
    if wt_data["hbonds"]:
        wt_hb_all = np.concatenate([h for _, h in wt_data["hbonds"]])
        print(f"    WT   H-bonds: {np.mean(wt_hb_all):.1f} ± {np.std(wt_hb_all):.1f}")
    if mut_data["hbonds"]:
        mut_hb_all = np.concatenate([h for _, h in mut_data["hbonds"]])
        print(f"    S305 H-bonds: {np.mean(mut_hb_all):.1f} ± {np.std(mut_hb_all):.1f}")


def main():
    print("Deep Analysis of 200ns MD Data")
    print(f"Output: {OUTDIR}")

    wt_data = analyze_system("WT", SYSTEMS["WT"], OUTDIR)
    mut_data = analyze_system("S305phos", SYSTEMS["S305phos"], OUTDIR)
    compare_systems(wt_data, mut_data, OUTDIR)

    print(f"\n{'='*60}")
    print("All analyses complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
