#!/usr/bin/env python3
"""
Dynamical Cross-Correlation Map (DCCM) from MD trajectory.

Usage:
  python scripts/run_dccm.py \
      --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
      --trajectories data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd \
      --selection "protein and name CA" \
      --name Hsap_WT \
      --outdir data/analysis/dccm

Multiple replicas:
  python scripts/run_dccm.py \
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


def load_ca_trajectory(prmtop, trajectories, selection="protein and name CA"):
    """Load trajectory and extract CA coordinates per residue."""
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    
    print(f"[DCCM] Loading: {prmtop}")
    print(f"[DCCM] Trajectories: {len(trajectories)}")
    
    u = mda.Universe(prmtop, trajectories)
    ca_atoms = u.select_atoms(selection)
    print(f"[DCCM] Selected {ca_atoms.n_atoms} CA atoms from {len(set(ca_atoms.resids))} residues")
    
    # Align to average structure
    ref = u.copy()
    align.AlignTraj(u, ref, select=selection, in_memory=True).run()
    
    # Collect CA coordinates: (n_frames, n_residues, 3)
    coords = []
    resids = ca_atoms.resids
    resnames = ca_atoms.resnames
    for ts in u.trajectory:
        coords.append(ca_atoms.positions.copy())
    
    coords = np.array(coords, dtype=np.float32)
    print(f"[DCCM] Collected {len(coords)} frames")
    return coords, resids, resnames


def compute_dccm(coords):
    """Compute dynamical cross-correlation map."""
    n_frames, n_res, _ = coords.shape
    
    # Mean positions
    mean_pos = coords.mean(axis=0)  # (n_res, 3)
    
    # Displacements from mean
    disp = coords - mean_pos  # (n_frames, n_res, 3)
    
    print(f"[DCCM] Computing correlation matrix ({n_res}×{n_res})...")
    
    # Compute DCCM
    dccm = np.zeros((n_res, n_res), dtype=np.float32)
    
    # Vectorized computation
    for i in range(n_res):
        di = disp[:, i, :]  # (n_frames, 3)
        di_norm = np.sqrt(np.mean((di ** 2).sum(axis=1)))
        if di_norm == 0:
            continue
        
        for j in range(i, n_res):
            dj = disp[:, j, :]  # (n_frames, 3)
            dj_norm = np.sqrt(np.mean((dj ** 2).sum(axis=1)))
            if dj_norm == 0:
                continue
            
            # <Δr_i · Δr_j>
            corr = np.mean((di * dj).sum(axis=1))
            dccm[i, j] = corr / (di_norm * dj_norm)
            dccm[j, i] = dccm[i, j]
    
    return dccm


def plot_dccm(dccm, resids, resnames, name, outdir, threshold=0.5):
    """Plot DCCM heatmap."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(dccm, cmap="RdBu_r", vmin=-1, vmax=1, 
                   aspect="auto", origin="lower")
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Cross-correlation", fontsize=11)
    
    # Labels
    n = len(resids)
    tick_step = max(1, n // 10)
    ticks = np.arange(0, n, tick_step)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([f"{resids[i]}\n{resnames[i]}" for i in ticks], 
                       rotation=45, ha="right", fontsize=6)
    ax.set_yticklabels([f"{resids[i]}\n{resnames[i]}" for i in ticks], 
                       fontsize=6)
    
    ax.set_xlabel("Residue", fontsize=11)
    ax.set_ylabel("Residue", fontsize=11)
    ax.set_title(f"Dynamical Cross-Correlation Map\n{name} ({n} residues)", fontsize=12)
    
    fig.tight_layout()
    path = outdir / f"{name}_dccm.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path.name}")


def extract_high_correlations(dccm, resids, resnames, threshold=0.7, min_sep=5):
    """Extract residue pairs with |correlation| > threshold."""
    n = len(resids)
    pairs = []
    
    for i in range(n):
        for j in range(i + min_sep, n):
            if abs(dccm[i, j]) >= threshold:
                pairs.append({
                    "res_i": int(resids[i]),
                    "name_i": resnames[i],
                    "res_j": int(resids[j]),
                    "name_j": resnames[j],
                    "corr": float(dccm[i, j]),
                    "sep": int(resids[j] - resids[i]),
                })
    
    pairs.sort(key=lambda x: abs(x["corr"]), reverse=True)
    return pairs


def plot_correlation_distribution(dccm, name, outdir):
    """Plot histogram of correlation values."""
    # Extract upper triangle (i < j)
    n = dccm.shape[0]
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(dccm[i, j])
    vals = np.array(vals)
    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(vals, bins=100, range=(-1, 1), color="steelblue", alpha=0.7, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", lw=0.5)
    ax.axvline(vals.mean(), color="red", ls="--", label=f"mean={vals.mean():.3f}")
    ax.set_xlabel("Cross-correlation")
    ax.set_ylabel("Count")
    ax.set_title(f"DCCM Distribution — {name}\n{n*(n-1)//2} residue pairs")
    ax.legend()
    ax.set_xlim(-1, 1)
    fig.tight_layout()
    path = outdir / f"{name}_dccm_distribution.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path.name}")


def main():
    parser = argparse.ArgumentParser(description="DCCM from MD trajectory")
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--trajectories", nargs="+", required=True)
    parser.add_argument("--selection", default="protein and name CA",
                        help="Atom selection for DCCM (default: CA)")
    parser.add_argument("--threshold", type=float, default=0.7,
                        help="Correlation threshold for pair extraction")
    parser.add_argument("--min-separation", type=int, default=5,
                        help="Minimum residue separation for pair extraction")
    parser.add_argument("--name", required=True)
    parser.add_argument("--outdir", default="data/analysis/dccm")
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print(f"[DCCM] Starting: {args.name}")
    
    # Load
    coords, resids, resnames = load_ca_trajectory(
        args.prmtop, args.trajectories, args.selection
    )
    
    # Compute
    dccm = compute_dccm(coords)
    
    # Plots
    print(f"\n[DCCM] Generating plots...")
    plot_dccm(dccm, resids, resnames, args.name, outdir)
    plot_correlation_distribution(dccm, args.name, outdir)
    
    # Extract high correlations
    pairs = extract_high_correlations(
        dccm, resids, resnames, 
        threshold=args.threshold, min_sep=args.min_separation
    )
    
    # Save
    np.savez(
        outdir / f"{args.name}_dccm.npz",
        dccm=dccm,
        resids=resids,
        resnames=resnames,
    )
    print(f"\n  Saved: {args.name}_dccm.npz")
    
    import json
    with open(outdir / f"{args.name}_highcorr.json", "w") as f:
        json.dump(pairs[:50], f, indent=2)
    print(f"  Saved: {args.name}_highcorr.json (top {min(50, len(pairs))} pairs)")
    
    # Summary
    print(f"\n{'='*50}")
    print(f"DCCM Summary — {args.name}")
    print(f"{'='*50}")
    print(f"  Residues: {len(resids)}")
    print(f"  Frames:   {len(coords)}")
    print(f"  Mean |C|: {np.abs(dccm).mean():.3f}")
    print(f"  |C|>0.7:  {sum(1 for p in pairs if abs(p['corr']) >= 0.7)}")
    print(f"  |C|>0.5:  {sum(1 for p in pairs if abs(p['corr']) >= 0.5)}")
    if pairs:
        print(f"\n  Top correlated pairs:")
        for p in pairs[:5]:
            print(f"    {p['name_i']}{p['res_i']} — {p['name_j']}{p['res_j']} "
                  f"(sep={p['sep']}): C={p['corr']:.3f}")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
