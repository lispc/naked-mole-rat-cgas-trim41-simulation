#!/usr/bin/env python3
"""
Joint PCA comparison of Hsap_WT vs Hsap_4mut.

Computes PCA on merged trajectories, then projects each system separately
onto the joint PC space. Generates comparison plots.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/pca"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_and_align(prmtop, dcds, name):
    """Load trajectories, align to average, return CA coords."""
    print(f"[{name}] Loading {len(dcds)} trajectories...")
    u = mda.Universe(str(prmtop), [str(d) for d in dcds])
    ca = u.select_atoms("protein and name CA")
    print(f"  {ca.n_atoms} CA atoms, {len(u.trajectory)} total frames")
    
    # Align to first frame
    ref = u.copy()
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    
    # Collect coordinates
    coords = np.array([ca.positions.copy() for _ in u.trajectory])
    resids = ca.resids
    return coords, resids


def compute_pca(coords, n_components=10):
    """PCA on flattened coordinate matrix."""
    n_frames, n_atoms, n_dims = coords.shape
    X = coords.reshape(n_frames, -1)  # (frames, atoms*3)
    
    # Mean-center
    mean = X.mean(axis=0)
    Xc = X - mean
    
    # Covariance and eigendecomposition
    print(f"[PCA] Computing covariance ({X.shape[1]}D)...")
    cov = np.cov(Xc, rowvar=False)
    
    print("[PCA] Eigendecomposition...")
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # Projections
    projections = Xc @ eigenvectors
    
    return eigenvalues, eigenvectors, projections, mean


def plot_variance(eigenvalues, name, outdir):
    """Plot eigenvalue spectrum and cumulative variance."""
    n = len(eigenvalues)
    total = eigenvalues.sum()
    explained = eigenvalues / total * 100
    cumsum = np.cumsum(explained)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    x = np.arange(1, min(n + 1, 21))
    axes[0].bar(x, eigenvalues[:len(x)], color="steelblue", alpha=0.7)
    axes[0].set_xlabel("Principal Component")
    axes[0].set_ylabel("Eigenvalue (Å²)")
    axes[0].set_title(f"PCA Eigenvalue Spectrum\n{name}")
    
    axes[1].plot(np.arange(1, n + 1), cumsum, "o-", color="darkgreen", markersize=3)
    axes[1].axhline(80, color="red", ls="--", alpha=0.5, label="80%")
    axes[1].axhline(90, color="orange", ls="--", alpha=0.5, label="90%")
    axes[1].set_xlabel("Principal Component")
    axes[1].set_ylabel("Cumulative Variance (%)")
    axes[1].set_title(f"Cumulative Variance\n{name}")
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(outdir / f"{name}_pca_variance.png", dpi=300)
    print(f"  Saved: {outdir / f'{name}_pca_variance.png'}")


def plot_comparison(projections_wt, projections_mut, name, outdir):
    """Plot PC1/PC2 projection comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # PC1 vs PC2
    ax = axes[0]
    ax.scatter(projections_wt[:, 0], projections_wt[:, 1], 
               alpha=0.3, s=5, color="tab:blue", label="WT")
    ax.scatter(projections_mut[:, 0], projections_mut[:, 1], 
               alpha=0.3, s=5, color="tab:red", label="4mut")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PC1 vs PC2")
    ax.legend()
    
    # PC1 vs PC3
    ax = axes[1]
    ax.scatter(projections_wt[:, 0], projections_wt[:, 2], 
               alpha=0.3, s=5, color="tab:blue", label="WT")
    ax.scatter(projections_mut[:, 0], projections_mut[:, 2], 
               alpha=0.3, s=5, color="tab:red", label="4mut")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC3")
    ax.set_title("PC1 vs PC3")
    ax.legend()
    
    # PC2 vs PC3
    ax = axes[2]
    ax.scatter(projections_wt[:, 1], projections_wt[:, 2], 
               alpha=0.3, s=5, color="tab:blue", label="WT")
    ax.scatter(projections_mut[:, 1], projections_mut[:, 2], 
               alpha=0.3, s=5, color="tab:red", label="4mut")
    ax.set_xlabel("PC2")
    ax.set_ylabel("PC3")
    ax.set_title("PC2 vs PC3")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(outdir / f"{name}_pca_comparison.png", dpi=300)
    print(f"  Saved: {outdir / f'{name}_pca_comparison.png'}")


def plot_pc1_loadings(eigenvectors, resids, name, outdir):
    """Plot per-residue loading on PC1."""
    # PC1 loadings: magnitude of eigenvector components per residue
    n_atoms = len(resids)
    pc1_load = np.zeros(n_atoms)
    for i in range(n_atoms):
        pc1_load[i] = np.linalg.norm(eigenvectors[i*3:(i+1)*3, 0])
    
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(resids, pc1_load, color="steelblue", alpha=0.7, width=1)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("PC1 Loading (Å)")
    ax.set_title(f"Per-Residue Contribution to PC1\n{name}")
    
    plt.tight_layout()
    plt.savefig(outdir / f"{name}_pc1_loadings.png", dpi=300)
    print(f"  Saved: {outdir / f'{name}_pc1_loadings.png'}")


def main():
    # Load both systems
    wt_coords, resids = load_and_align(
        BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        [
            BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        ],
        "Hsap_WT"
    )
    
    mut_coords, _ = load_and_align(
        BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
        [
            BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
        ],
        "Hsap_4mut"
    )
    
    # Merge for joint PCA
    print(f"\n[Joint PCA] Merging {wt_coords.shape[0]} + {mut_coords.shape[0]} = {wt_coords.shape[0] + mut_coords.shape[0]} frames")
    merged_coords = np.vstack([wt_coords, mut_coords])
    
    # Compute joint PCA
    eigenvalues, eigenvectors, projections, mean = compute_pca(merged_coords, n_components=10)
    
    # Split projections back
    n_wt = wt_coords.shape[0]
    proj_wt = projections[:n_wt]
    proj_mut = projections[n_wt:]
    
    # Statistics
    explained = eigenvalues / eigenvalues.sum() * 100
    print(f"\n{'='*50}")
    print("Joint PCA Summary — Hsap_WT + Hsap_4mut")
    print(f"{'='*50}")
    for i in range(min(5, len(eigenvalues))):
        print(f"  PC{i+1}: λ={eigenvalues[i]:.2f} ({explained[i]:.1f}% variance)")
    print(f"{'='*50}")
    
    # Centroid separation
    wt_centroid = proj_wt[:, :2].mean(axis=0)
    mut_centroid = proj_mut[:, :2].mean(axis=0)
    separation = np.linalg.norm(wt_centroid - mut_centroid)
    print(f"\n[Comparison] PC1/PC2 centroid separation: {separation:.2f} Å")
    
    # Overlap (fraction of WT points within 1 std of mut centroid)
    mut_std = proj_mut[:, :2].std(axis=0)
    distances = np.linalg.norm((proj_wt[:, :2] - mut_centroid) / mut_std, axis=1)
    overlap = np.mean(distances < 2.0) * 100
    print(f"[Comparison] WT frames within 2σ of 4mut centroid: {overlap:.1f}%")
    
    # Plots
    print("\n[Plots] Generating...")
    plot_variance(eigenvalues, "Hsap_joint", OUTDIR)
    plot_comparison(proj_wt, proj_mut, "Hsap_joint", OUTDIR)
    plot_pc1_loadings(eigenvectors, resids, "Hsap_joint", OUTDIR)
    
    # Save data
    np.savez(OUTDIR / "pca_joint_data.npz",
             eigenvalues=eigenvalues,
             eigenvectors=eigenvectors,
             mean=mean,
             proj_wt=proj_wt,
             proj_mut=proj_mut,
             resids=resids)
    print(f"\nSaved: {OUTDIR / 'pca_joint_data.npz'}")


if __name__ == "__main__":
    main()
