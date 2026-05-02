#!/usr/bin/env python3
"""
Principal Component Analysis (PCA) of MD trajectory.

Usage:
  python scripts/run_pca.py \
      --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
      --trajectories data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd \
      --selection "protein and backbone" \
      --n-components 10 \
      --name Hsap_WT \
      --outdir data/analysis/pca

Multiple replicas:
  python scripts/run_pca.py \
      --trajectories rep1.dcd rep2.dcd rep3.dcd \
      --name Hsap_WT_merged
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_trajectory(prmtop, trajectories, selection):
    """Load and align trajectory using MDAnalysis."""
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    
    print(f"[PCA] Loading topology: {prmtop}")
    print(f"[PCA] Loading {len(trajectories)} trajectory file(s)")
    
    u = mda.Universe(prmtop, trajectories)
    atoms = u.select_atoms(selection)
    print(f"[PCA] Selected {atoms.n_atoms} atoms ({selection})")
    
    # Align to first frame
    ref = u.copy()
    align.AlignTraj(u, ref, select=selection, in_memory=True).run()
    
    # Collect coordinates
    coords = []
    for ts in u.trajectory:
        coords.append(atoms.positions.copy().flatten())
    
    coords = np.array(coords, dtype=np.float32)
    print(f"[PCA] Collected {len(coords)} frames × {atoms.n_atoms} atoms")
    return coords, atoms


def run_pca(coords, n_components=10):
    """Perform PCA on coordinate matrix."""
    # Mean-center
    mean = coords.mean(axis=0)
    X = coords - mean
    
    # Covariance matrix
    print(f"[PCA] Computing covariance matrix ({coords.shape[1]}D)...")
    cov = np.cov(X, rowvar=False)
    
    # Eigendecomposition
    print(f"[PCA] Eigendecomposition...")
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # Projections
    projections = X @ eigenvectors
    
    return eigenvalues, eigenvectors, projections, mean


def plot_variance(eigenvalues, name, outdir):
    """Plot eigenvalue spectrum and cumulative variance."""
    n = len(eigenvalues)
    total = eigenvalues.sum()
    explained = eigenvalues / total * 100
    cumsum = np.cumsum(explained)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # Eigenvalue bar plot
    ax = axes[0]
    x = np.arange(1, min(n + 1, 21))
    ax.bar(x, eigenvalues[:len(x)], color="steelblue", alpha=0.7)
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Eigenvalue (Å²)")
    ax.set_title(f"PCA Eigenvalue Spectrum\n{name}")
    
    # Cumulative variance
    ax = axes[1]
    ax.plot(np.arange(1, n + 1), cumsum, "o-", color="darkgreen", markersize=3)
    ax.axhline(80, color="red", ls="--", alpha=0.5, label="80%")
    ax.axhline(90, color="orange", ls="--", alpha=0.5, label="90%")
    ax.set_xlabel("Number of PCs")
    ax.set_ylabel("Cumulative Variance (%)")
    ax.set_title(f"Explained Variance\n{name}")
    ax.legend()
    ax.set_xlim(1, min(n, 50))
    
    fig.tight_layout()
    path = outdir / f"{name}_pca_variance.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path.name}")


def plot_pc_projections(projections, name, outdir, max_pc=4):
    """Plot 2D projections of PC pairs."""
    n_plots = min(max_pc, projections.shape[1])
    n_pairs = n_plots * (n_plots - 1) // 2
    if n_pairs == 0:
        return
    
    ncols = min(n_pairs, 3)
    nrows = (n_pairs + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    if n_pairs == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    pair_idx = 0
    for i in range(n_plots):
        for j in range(i + 1, n_plots):
            ax = axes[pair_idx]
            ax.scatter(projections[:, i], projections[:, j], s=1, alpha=0.3, c=np.arange(len(projections)), cmap="viridis")
            ax.set_xlabel(f"PC{i+1}")
            ax.set_ylabel(f"PC{j+1}")
            ax.set_title(f"PC{i+1} vs PC{j+1}")
            pair_idx += 1
    
    # Hide unused subplots
    for k in range(pair_idx, len(axes)):
        axes[k].axis("off")
    
    fig.suptitle(f"PCA Projections — {name}", fontsize=12)
    fig.tight_layout()
    path = outdir / f"{name}_pca_projections.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path.name}")


def write_pc_modes(eigenvectors, mean, atoms, name, outdir, n_modes=3):
    """Write extreme conformations along top PCs as PDB."""
    import MDAnalysis as mda
    
    n_atoms = atoms.n_atoms
    mean_coords = mean.reshape(n_atoms, 3)
    
    # Build residue-level attributes
    resids = []
    resnames = []
    segids = []
    for res in atoms.residues:
        for _ in res.atoms:
            resids.append(res.resid)
            resnames.append(res.resname)
            segids.append(res.segid if hasattr(res, 'segid') else '')
    
    for pc in range(min(n_modes, eigenvectors.shape[1])):
        vec = eigenvectors[:, pc].reshape(n_atoms, 3)
        scale = 2.0  # Angstroms
        
        # Write +/- extremes
        for sign, suffix in [(+1, "plus"), (-1, "minus")]:
            coords = mean_coords + sign * scale * vec
            u_tmp = mda.Universe.empty(n_atoms, n_residues=atoms.n_residues, 
                                       atom_resindex=[a.residue.ix for a in atoms],
                                       trajectory=True)
            u_tmp.add_TopologyAttr("names", atoms.names)
            u_tmp.add_TopologyAttr("resnames", resnames[:atoms.n_residues])
            u_tmp.add_TopologyAttr("resids", resids[:atoms.n_residues])
            u_tmp.atoms.positions = coords
            
            pdb_path = outdir / f"{name}_PC{pc+1}_{suffix}.pdb"
            u_tmp.atoms.write(str(pdb_path))
    
    print(f"  Saved: {n_modes} PC mode PDBs (±2Å displacement)")


def main():
    parser = argparse.ArgumentParser(description="PCA of MD trajectory")
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--trajectories", nargs="+", required=True)
    parser.add_argument("--selection", default="protein and backbone",
                        help="Atom selection for PCA (default: protein backbone)")
    parser.add_argument("--n-components", type=int, default=10)
    parser.add_argument("--name", required=True)
    parser.add_argument("--outdir", default="data/analysis/pca")
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print(f"[PCA] Starting analysis: {args.name}")
    
    # Load
    coords, atoms = load_trajectory(args.prmtop, args.trajectories, args.selection)
    
    # PCA
    eigenvalues, eigenvectors, projections, mean = run_pca(coords, args.n_components)
    
    # Plots
    print(f"\n[PCA] Generating plots...")
    plot_variance(eigenvalues, args.name, outdir)
    plot_pc_projections(projections, args.name, outdir)
    
    # Mode PDBs
    write_pc_modes(eigenvectors, mean, atoms, args.name, outdir)
    
    # Save data
    np.savez(
        outdir / f"{args.name}_pca.npz",
        eigenvalues=eigenvalues,
        projections=projections,
        mean=mean,
    )
    print(f"\n[PCA] Saved: {args.name}_pca.npz")
    
    # Summary
    total_var = eigenvalues.sum()
    print(f"\n{'='*50}")
    print(f"PCA Summary — {args.name}")
    print(f"{'='*50}")
    for i in range(min(args.n_components, 5)):
        pct = eigenvalues[i] / total_var * 100
        print(f"  PC{i+1}: λ={eigenvalues[i]:.2f} ({pct:.1f}% variance)")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
