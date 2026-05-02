#!/usr/bin/env python3
"""
Compute ΔDCCM (4mut - WT) for Hsap systems.
Merges 3 reps per system, aligns, computes DCCM, then difference map.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/delta_dccm"
OUTDIR.mkdir(parents=True, exist_ok=True)


def compute_dccm(prmtop, dcds, name):
    print(f"[{name}] Loading {len(dcds)} trajectories...")
    u = mda.Universe(str(prmtop), [str(d) for d in dcds])
    ca = u.select_atoms("protein and name CA")
    n_res = ca.n_atoms
    print(f"  {n_res} CA atoms, {len(u.trajectory)} total frames")
    
    # Align to first frame
    ref = u.copy()
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    
    # Collect coordinates: (n_frames, n_residues, 3)
    coords = np.array([ca.positions.copy() for _ in u.trajectory])
    
    # Compute displacement from mean
    mean_pos = coords.mean(axis=0)
    disp = coords - mean_pos  # (frames, residues, 3)
    
    # Compute covariance matrix
    n_frames = coords.shape[0]
    cov = np.zeros((n_res, n_res))
    for i in range(n_res):
        for j in range(i, n_res):
            # dot product of displacement vectors across all frames and xyz
            val = np.sum(disp[:, i, :] * disp[:, j, :]) / n_frames
            cov[i, j] = val
            cov[j, i] = val
    
    # Convert to correlation matrix
    std = np.sqrt(np.diag(cov))
    std_mat = np.outer(std, std)
    std_mat[std_mat == 0] = 1  # avoid division by zero
    corr = cov / std_mat
    
    resids = ca.resids
    print(f"[{name}] DCCM computed. Mean |C| = {np.abs(corr).mean():.3f}")
    return resids, corr


def main():
    systems = {
        "Hsap_WT": {
            "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
            "dcds": [
                BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
                BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
                BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
            ],
        },
        "Hsap_4mut": {
            "prmtop": BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
            "dcds": [
                BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
                BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
                BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
            ],
        },
    }
    
    results = {}
    for sys_name, paths in systems.items():
        resids, corr = compute_dccm(paths["prmtop"], paths["dcds"], sys_name)
        results[sys_name] = {"resids": resids, "corr": corr}
    
    # ΔDCCM = 4mut - WT
    wt_corr = results["Hsap_WT"]["corr"]
    mut_corr = results["Hsap_4mut"]["corr"]
    delta = mut_corr - wt_corr
    resids = results["Hsap_WT"]["resids"]
    
    print(f"\n[ΔDCCM] Max increase: +{delta.max():.3f}")
    print(f"[ΔDCCM] Max decrease: {delta.min():.3f}")
    print(f"[ΔDCCM] |Δ| > 0.3: {np.sum(np.abs(delta) > 0.3)} pairs")
    print(f"[ΔDCCM] |Δ| > 0.5: {np.sum(np.abs(delta) > 0.5)} pairs")
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    vmax = 1.0
    im0 = axes[0].imshow(wt_corr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
    axes[0].set_title("WT DCCM")
    axes[0].set_xlabel("Residue")
    axes[0].set_ylabel("Residue")
    plt.colorbar(im0, ax=axes[0], fraction=0.046)
    
    im1 = axes[1].imshow(mut_corr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
    axes[1].set_title("4mut DCCM")
    axes[1].set_xlabel("Residue")
    axes[1].set_ylabel("Residue")
    plt.colorbar(im1, ax=axes[1], fraction=0.046)
    
    dmax = max(np.abs(delta).max(), 0.5)
    im2 = axes[2].imshow(delta, cmap="RdBu_r", vmin=-dmax, vmax=dmax, origin="lower")
    axes[2].set_title("ΔDCCM (4mut − WT)")
    axes[2].set_xlabel("Residue")
    axes[2].set_ylabel("Residue")
    plt.colorbar(im2, ax=axes[2], fraction=0.046)
    
    plt.tight_layout()
    plt.savefig(OUTDIR / "delta_dccm_Hsap_WT_vs_4mut.png", dpi=300)
    print(f"\nSaved: {OUTDIR / 'delta_dccm_Hsap_WT_vs_4mut.png'}")
    
    # Save data
    np.savez(OUTDIR / "delta_dccm_data.npz",
             resids=resids,
             wt_corr=wt_corr,
             mut_corr=mut_corr,
             delta=delta)
    print(f"Saved: {OUTDIR / 'delta_dccm_data.npz'}")


if __name__ == "__main__":
    main()
