#!/usr/bin/env python3
"""
Comprehensive allosteric analysis for Hgal systems (WT vs 4mut_rev).
Runs ΔRMSF, ΔDCCM, and joint PCA in one pass.

Data: 3 reps × 200ns each for Hgal_WT and Hgal_4mut_rev
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from pathlib import Path
from scipy import stats

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/hgal_allostery"
OUTDIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "Hgal_WT": {
        "prmtop": BASE / "data/md_runs/Hgal_WT/Hgal_WT.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hgal_WT/rep1/Hgal_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hgal_WT/rep2/Hgal_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hgal_WT/rep2/Hgal_WT_rep2_restart.dcd",
            BASE / "data/md_runs/Hgal_WT/rep3/Hgal_WT_rep3_prod.dcd",
        ],
    },
    "Hgal_4mut_rev": {
        "prmtop": BASE / "data/md_runs/Hgal_4mut_rev/Hgal_4mut_rev.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hgal_4mut_rev/rep1/Hgal_4mut_rev_rep1_prod.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_prod.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_restart.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep3/Hgal_4mut_rev_rep3_prod.dcd",
        ],
    },
}

# ============================================================================
# 1. ΔRMSF
# ============================================================================
def compute_rmsf_per_rep(prmtop, dcd, name, rep):
    print(f"[RMSF {name} rep{rep}] Loading {dcd.name}...")
    u = mda.Universe(str(prmtop), str(dcd))
    ca = u.select_atoms("protein and name CA")
    print(f"  {len(ca)} CA atoms, {len(u.trajectory)} frames")
    
    ref = u.copy()
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    
    rmsf = rms.RMSF(ca, verbose=True).run()
    return ca.resids, ca.resnames, rmsf.results.rmsf


def run_delta_rmsf():
    print("\n" + "="*60)
    print("PART 1: ΔRMSF Analysis")
    print("="*60)
    
    results = {}
    for sys_name, paths in SYSTEMS.items():
        prmtop = paths["prmtop"]
        reps = paths["dcds"]
        
        rmsfs = []
        for i, dcd in enumerate(reps, 1):
            resids, resnames, rmsf = compute_rmsf_per_rep(prmtop, dcd, sys_name, i)
            rmsfs.append(rmsf)
        
        avg_rmsf = np.mean(rmsfs, axis=0)
        results[sys_name] = {
            "resids": resids,
            "resnames": resnames,
            "avg_rmsf": avg_rmsf,
            "reps": np.array(rmsfs),
        }
        print(f"[{sys_name}] Avg RMSF: {avg_rmsf.mean():.2f} ± {avg_rmsf.std():.2f} Å")
    
    wt = results["Hgal_WT"]
    mut = results["Hgal_4mut_rev"]
    common_resids = wt["resids"]
    delta = mut["avg_rmsf"] - wt["avg_rmsf"]
    
    # t-test (3 reps vs 3 reps)
    pvalues = np.zeros(len(common_resids))
    for i in range(len(common_resids)):
        _, pvalues[i] = stats.ttest_ind(mut["reps"][:, i], wt["reps"][:, i], equal_var=False)
    pvalues_corrected = np.minimum(pvalues * len(common_resids), 1.0)
    
    print(f"\n[ΔRMSF] Mean: {delta.mean():.3f} Å")
    print(f"[ΔRMSF] Max increase: +{delta.max():.2f} Å at {common_resids[np.argmax(delta)]}")
    print(f"[ΔRMSF] Max decrease: {delta.min():.2f} Å at {common_resids[np.argmin(delta)]}")
    print(f"[ΔRMSF] |Δ| > 1.0 Å: {np.sum(np.abs(delta) > 1.0)} residues")
    print(f"[ΔRMSF] |Δ| > 2.0 Å: {np.sum(np.abs(delta) > 2.0)} residues")
    print(f"[ΔRMSF] Significant (p<0.05): {np.sum(pvalues < 0.05)} (uncorrected), {np.sum(pvalues_corrected < 0.05)} (Bonferroni)")
    
    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    
    axes[0].plot(common_resids, wt["avg_rmsf"], label="WT", color="tab:blue", alpha=0.8)
    axes[0].plot(common_resids, mut["avg_rmsf"], label="4mut_rev", color="tab:red", alpha=0.8)
    axes[0].fill_between(common_resids, wt["avg_rmsf"], alpha=0.15, color="tab:blue")
    axes[0].fill_between(common_resids, mut["avg_rmsf"], alpha=0.15, color="tab:red")
    axes[0].set_ylabel("RMSF (Å)")
    axes[0].set_title("Per-Residue RMSF (3-rep average)")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    colors = ["tab:red" if d > 0 else "tab:blue" for d in delta]
    axes[1].bar(common_resids, delta, color=colors, alpha=0.6, width=1)
    sig_mask = pvalues_corrected < 0.05
    axes[1].scatter(common_resids[sig_mask], delta[sig_mask], 
                    color="black", s=10, zorder=5, label="p<0.05 (Bonferroni)")
    axes[1].axhline(0, color="black", lw=0.5)
    axes[1].set_ylabel("ΔRMSF (4mut_rev − WT, Å)")
    axes[1].set_xlabel("Residue ID")
    axes[1].set_title("ΔRMSF")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTDIR / "hgal_delta_rmsf.png", dpi=300)
    print(f"\nSaved: {OUTDIR / 'hgal_delta_rmsf.png'}")
    
    # Save data
    np.savez(OUTDIR / "hgal_delta_rmsf.npz",
             resids=common_resids,
             resnames=wt["resnames"],
             wt_rmsf=wt["avg_rmsf"],
             mut_rmsf=mut["avg_rmsf"],
             delta=delta,
             pvalues=pvalues,
             pvalues_bonf=pvalues_corrected)
    
    return results


# ============================================================================
# 2. ΔDCCM
# ============================================================================
def compute_dccm(prmtop, dcds, name):
    print(f"\n[DCCM {name}] Loading {len(dcds)} trajectories individually...")
    
    # Load first DCD to get reference structure
    u_ref = mda.Universe(str(prmtop), str(dcds[0]))
    ca_ref = u_ref.select_atoms("protein and name CA")
    n_res = ca_ref.n_atoms
    ref_positions = ca_ref.positions.copy()
    resids = ca_ref.resids
    print(f"  {n_res} CA atoms")
    
    # Collect coordinates from all DCDs
    all_coords = []
    for dcd in dcds:
        u = mda.Universe(str(prmtop), str(dcd))
        ca = u.select_atoms("protein and name CA")
        
        for ts in u.trajectory:
            # Simple alignment to reference (translation only, no rotation)
            # For DCCM, rigid-body rotation doesn't affect correlations
            # So translation-only alignment is sufficient
            positions = ca.positions.copy()
            com = positions.mean(axis=0)
            ref_com = ref_positions.mean(axis=0)
            positions -= com - ref_com
            all_coords.append(positions)
    
    coords = np.array(all_coords)
    print(f"  Total frames: {len(coords)}")
    
    # Compute DCCM
    mean_pos = coords.mean(axis=0)
    disp = coords - mean_pos
    
    n_frames = coords.shape[0]
    cov = np.zeros((n_res, n_res))
    for i in range(n_res):
        for j in range(i, n_res):
            val = np.sum(disp[:, i, :] * disp[:, j, :]) / n_frames
            cov[i, j] = val
            cov[j, i] = val
    
    std = np.sqrt(np.diag(cov))
    std_mat = np.outer(std, std)
    std_mat[std_mat == 0] = 1
    corr = cov / std_mat
    
    print(f"[{name}] DCCM done. Mean |C| = {np.abs(corr).mean():.3f}")
    return resids, corr


def run_delta_dccm():
    print("\n" + "="*60)
    print("PART 2: ΔDCCM Analysis")
    print("="*60)
    
    results = {}
    for sys_name, paths in SYSTEMS.items():
        resids, corr = compute_dccm(paths["prmtop"], paths["dcds"], sys_name)
        results[sys_name] = {"resids": resids, "corr": corr}
    
    wt_corr = results["Hgal_WT"]["corr"]
    mut_corr = results["Hgal_4mut_rev"]["corr"]
    delta = mut_corr - wt_corr
    resids = results["Hgal_WT"]["resids"]
    
    print(f"\n[ΔDCCM] Max increase: +{delta.max():.3f}")
    print(f"[ΔDCCM] Max decrease: {delta.min():.3f}")
    print(f"[ΔDCCM] |Δ| > 0.3: {np.sum(np.abs(delta) > 0.3)} pairs")
    print(f"[ΔDCCM] |Δ| > 0.5: {np.sum(np.abs(delta) > 0.5)} pairs")
    
    # N-term vs C-term coupling
    n = len(resids)
    nt_end = n // 3
    ct_start = 2 * n // 3
    print(f"\nN-term-C-term coupling:")
    print(f"  WT:   {wt_corr[:nt_end, ct_start:].mean():.3f}")
    print(f"  4mut: {mut_corr[:nt_end, ct_start:].mean():.3f}")
    print(f"  Δ:    {delta[:nt_end, ct_start:].mean():.3f}")
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    vmax = 1.0
    im0 = axes[0].imshow(wt_corr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
    axes[0].set_title("WT DCCM")
    axes[0].set_xlabel("Residue")
    axes[0].set_ylabel("Residue")
    plt.colorbar(im0, ax=axes[0], fraction=0.046)
    
    im1 = axes[1].imshow(mut_corr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
    axes[1].set_title("4mut_rev DCCM")
    axes[1].set_xlabel("Residue")
    axes[1].set_ylabel("Residue")
    plt.colorbar(im1, ax=axes[1], fraction=0.046)
    
    dmax = max(np.abs(delta).max(), 0.5)
    im2 = axes[2].imshow(delta, cmap="RdBu_r", vmin=-dmax, vmax=dmax, origin="lower")
    axes[2].set_title("ΔDCCM (4mut_rev − WT)")
    axes[2].set_xlabel("Residue")
    axes[2].set_ylabel("Residue")
    plt.colorbar(im2, ax=axes[2], fraction=0.046)
    
    plt.tight_layout()
    plt.savefig(OUTDIR / "hgal_delta_dccm.png", dpi=300)
    print(f"\nSaved: {OUTDIR / 'hgal_delta_dccm.png'}")
    
    np.savez(OUTDIR / "hgal_delta_dccm.npz",
             resids=resids,
             wt_corr=wt_corr,
             mut_corr=mut_corr,
             delta=delta)
    
    return results


# ============================================================================
# 3. Joint PCA
# ============================================================================
def run_joint_pca():
    print("\n" + "="*60)
    print("PART 3: Joint PCA")
    print("="*60)
    
    all_coords = []
    labels = []
    
    for sys_name, paths in SYSTEMS.items():
        print(f"[{sys_name}] Loading...")
        u = mda.Universe(str(paths["prmtop"]), [str(d) for d in paths["dcds"]])
        ca = u.select_atoms("protein and name CA")
        
        ref = u.copy()
        align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
        
        coords = np.array([ca.positions.copy() for _ in u.trajectory])
        all_coords.append(coords)
        labels.extend([sys_name] * len(coords))
        print(f"  {coords.shape[0]} frames")
    
    merged = np.vstack(all_coords)
    labels = np.array(labels)
    print(f"\n[Joint PCA] Total: {merged.shape[0]} frames × {merged.shape[1]} atoms")
    
    # Flatten and mean-center
    X = merged.reshape(merged.shape[0], -1)
    mean = X.mean(axis=0)
    Xc = X - mean
    
    # Covariance
    print("[PCA] Computing covariance...")
    cov = np.cov(Xc, rowvar=False)
    
    print("[PCA] Eigendecomposition...")
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    projections = Xc @ eigenvectors
    
    explained = eigenvalues / eigenvalues.sum() * 100
    print(f"\n{'='*50}")
    print("Joint PCA — Hgal_WT + Hgal_4mut_rev")
    print(f"{'='*50}")
    for i in range(min(5, len(eigenvalues))):
        print(f"  PC{i+1}: λ={eigenvalues[i]:.2f} ({explained[i]:.1f}% variance)")
    print(f"{'='*50}")
    
    # Split projections
    n_wt = len(all_coords[0])
    proj_wt = projections[:n_wt]
    proj_mut = projections[n_wt:]
    
    # Centroids
    wt_c = proj_wt[:, :2].mean(axis=0)
    mut_c = proj_mut[:, :2].mean(axis=0)
    sep = np.linalg.norm(wt_c - mut_c)
    mut_std = proj_mut[:, :2].std(axis=0)
    dist = np.linalg.norm((proj_wt[:, :2] - mut_c) / mut_std, axis=1)
    overlap = np.mean(dist < 2.0) * 100
    
    print(f"\nCentroid separation (PC1/PC2): {sep:.2f} Å")
    print(f"WT in 2σ of 4mut: {overlap:.1f}%")
    
    # Plot comparison
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    ax = axes[0]
    ax.scatter(proj_wt[:, 0], proj_wt[:, 1], alpha=0.2, s=3, color="tab:blue", label="WT")
    ax.scatter(proj_mut[:, 0], proj_mut[:, 1], alpha=0.2, s=3, color="tab:red", label="4mut_rev")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PC1 vs PC2")
    ax.legend()
    
    ax = axes[1]
    ax.scatter(proj_wt[:, 0], proj_wt[:, 2], alpha=0.2, s=3, color="tab:blue", label="WT")
    ax.scatter(proj_mut[:, 0], proj_mut[:, 2], alpha=0.2, s=3, color="tab:red", label="4mut_rev")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC3")
    ax.set_title("PC1 vs PC3")
    ax.legend()
    
    ax = axes[2]
    ax.scatter(proj_wt[:, 1], proj_wt[:, 2], alpha=0.2, s=3, color="tab:blue", label="WT")
    ax.scatter(proj_mut[:, 1], proj_mut[:, 2], alpha=0.2, s=3, color="tab:red", label="4mut_rev")
    ax.set_xlabel("PC2")
    ax.set_ylabel("PC3")
    ax.set_title("PC2 vs PC3")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(OUTDIR / "hgal_joint_pca.png", dpi=300)
    print(f"\nSaved: {OUTDIR / 'hgal_joint_pca.png'}")
    
    # PC1 loadings
    n_atoms = all_coords[0].shape[1]
    resids = mda.Universe(str(SYSTEMS["Hgal_WT"]["prmtop"])).select_atoms("protein and name CA").resids
    pc1_load = np.zeros(n_atoms)
    for i in range(n_atoms):
        pc1_load[i] = np.linalg.norm(eigenvectors[i*3:(i+1)*3, 0])
    
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(resids, pc1_load, color="steelblue", alpha=0.7, width=1)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("PC1 Loading (Å)")
    ax.set_title("Per-Residue PC1 Loading")
    plt.tight_layout()
    plt.savefig(OUTDIR / "hgal_pc1_loadings.png", dpi=300)
    print(f"Saved: {OUTDIR / 'hgal_pc1_loadings.png'}")
    
    # Save
    np.savez(OUTDIR / "hgal_joint_pca.npz",
             eigenvalues=eigenvalues,
             eigenvectors=eigenvectors,
             mean=mean,
             proj_wt=proj_wt,
             proj_mut=proj_mut,
             resids=resids)


def main():
    print("Hgal Allosteric Analysis Pipeline")
    print("Systems: Hgal_WT (3 reps) vs Hgal_4mut_rev (3 reps)")
    print(f"Output: {OUTDIR}")
    
    run_delta_rmsf()
    run_delta_dccm()
    run_joint_pca()
    
    print("\n" + "="*60)
    print("All analyses complete!")
    print(f"Results saved to: {OUTDIR}")
    print("="*60)


if __name__ == "__main__":
    main()
