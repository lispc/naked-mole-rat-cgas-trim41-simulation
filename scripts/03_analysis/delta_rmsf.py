#!/usr/bin/env python3
"""
Compute ΔRMSF (4mut - WT) for Hsap systems across 3 reps each.
Includes per-residue Welch t-test for statistical significance.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from pathlib import Path
from scipy import stats

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/delta_rmsf"
OUTDIR.mkdir(parents=True, exist_ok=True)

def compute_rmsf_per_rep(prmtop, dcd, name, rep):
    print(f"[{name} rep{rep}] Loading...")
    u = mda.Universe(str(prmtop), str(dcd))
    ca = u.select_atoms("protein and name CA")
    print(f"  {len(ca)} CA atoms, {len(u.trajectory)} frames")
    
    # Align to first frame
    ref = u.copy()
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    
    # RMSF
    rmsf = rms.RMSF(ca, verbose=True).run()
    return ca.resids, ca.resnames, rmsf.results.rmsf

def main():
    systems = {
        "Hsap_WT": {
            "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
            "rep1": BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            "rep2": BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            "rep3": BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        },
        "Hsap_4mut": {
            "prmtop": BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
            "rep1": BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
            "rep2": BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
            "rep3": BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
        },
    }
    
    results = {}
    for sys_name, paths in systems.items():
        prmtop = paths["prmtop"]
        reps = [paths[f"rep{i}"] for i in [1, 2, 3]]
        
        rmsfs = []
        for i, dcd in enumerate(reps, 1):
            resids, resnames, rmsf = compute_rmsf_per_rep(prmtop, dcd, sys_name, i)
            rmsfs.append(rmsf)
        
        avg_rmsf = np.mean(rmsfs, axis=0)
        results[sys_name] = {
            "resids": resids,
            "resnames": resnames,
            "avg_rmsf": avg_rmsf,
            "reps": np.array(rmsfs),  # (3, n_residues)
        }
        print(f"[{sys_name}] Average RMSF: {avg_rmsf.mean():.2f} ± {avg_rmsf.std():.2f} Å")
    
    # ΔRMSF = 4mut - WT
    wt = results["Hsap_WT"]
    mut = results["Hsap_4mut"]
    
    common_resids = wt["resids"]
    delta = mut["avg_rmsf"] - wt["avg_rmsf"]
    
    # Welch t-test per residue (3 reps vs 3 reps)
    wt_reps = wt["reps"]   # (3, n_res)
    mut_reps = mut["reps"] # (3, n_res)
    
    pvalues = np.zeros(len(common_resids))
    for i in range(len(common_resids)):
        _, pvalues[i] = stats.ttest_ind(mut_reps[:, i], wt_reps[:, i], equal_var=False)
    
    # Bonferroni correction
    pvalues_corrected = np.minimum(pvalues * len(common_resids), 1.0)
    
    print(f"\n[ΔRMSF] Mean: {delta.mean():.3f} Å")
    print(f"[ΔRMSF] Max increase: +{delta.max():.2f} Å at resid {common_resids[np.argmax(delta)]}")
    print(f"[ΔRMSF] Max decrease: {delta.min():.2f} Å at resid {common_resids[np.argmin(delta)]}")
    print(f"[ΔRMSF] |Δ| > 1.0 Å: {np.sum(np.abs(delta) > 1.0)} residues")
    print(f"[ΔRMSF] |Δ| > 2.0 Å: {np.sum(np.abs(delta) > 2.0)} residues")
    print(f"[ΔRMSF] Significant (p<0.05, uncorrected): {np.sum(pvalues < 0.05)} residues")
    print(f"[ΔRMSF] Significant (p<0.05, Bonferroni): {np.sum(pvalues_corrected < 0.05)} residues")
    
    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    
    # Top: individual RMSF
    axes[0].plot(common_resids, wt["avg_rmsf"], label="WT", color="tab:blue", alpha=0.8)
    axes[0].plot(common_resids, mut["avg_rmsf"], label="4mut", color="tab:red", alpha=0.8)
    axes[0].fill_between(common_resids, wt["avg_rmsf"], alpha=0.2, color="tab:blue")
    axes[0].fill_between(common_resids, mut["avg_rmsf"], alpha=0.2, color="tab:red")
    axes[0].set_ylabel("RMSF (Å)")
    axes[0].set_title("Hsap WT vs 4mut — CA RMSF (3 reps average)")
    axes[0].legend()
    axes[0].axvline(x=218, color="gray", linestyle="--", alpha=0.5, label="TRIM41|cGAS")
    
    # Bottom: ΔRMSF with significance
    colors = ["tab:red" if d > 0 else "tab:blue" for d in delta]
    axes[1].bar(common_resids, delta, color=colors, width=1.0, alpha=0.7)
    # Mark significant residues
    sig_mask = pvalues_corrected < 0.05
    if np.any(sig_mask):
        axes[1].scatter(common_resids[sig_mask], delta[sig_mask], 
                       color="black", s=10, zorder=5, label="p<0.05 (Bonferroni)")
    axes[1].axhline(y=0, color="black", linewidth=0.5)
    axes[1].axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="+1.0 Å")
    axes[1].axhline(y=-1.0, color="blue", linestyle="--", alpha=0.5, label="-1.0 Å")
    axes[1].set_xlabel("Residue ID")
    axes[1].set_ylabel("ΔRMSF (4mut − WT, Å)")
    axes[1].set_title("ΔRMSF — Positive = 4mut more flexible")
    axes[1].legend(fontsize=8)
    axes[1].axvline(x=218, color="gray", linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(OUTDIR / "delta_rmsf_Hsap_WT_vs_4mut.png", dpi=300)
    print(f"\nSaved: {OUTDIR / 'delta_rmsf_Hsap_WT_vs_4mut.png'}")
    
    # Save data
    np.savez(OUTDIR / "delta_rmsf_data.npz",
             resids=common_resids,
             resnames=np.array(wt["resnames"]),
             wt_rmsf=wt["avg_rmsf"],
             mut_rmsf=mut["avg_rmsf"],
             delta=delta,
             pvalues=pvalues,
             pvalues_corrected=pvalues_corrected,
             wt_reps=wt_reps,
             mut_reps=mut_reps)
    print(f"Saved: {OUTDIR / 'delta_rmsf_data.npz'}")
    
    # Save CSV
    csv_lines = ["resid,resname,wt_rmsf_mean,wt_rmsf_std,mut_rmsf_mean,mut_rmsf_std,delta,pvalue,pvalue_bonferroni,significant"]
    for i in range(len(common_resids)):
        sig = "YES" if pvalues_corrected[i] < 0.05 else "NO"
        csv_lines.append(
            f"{int(common_resids[i])},{wt['resnames'][i]},"
            f"{wt_reps[:,i].mean():.4f},{wt_reps[:,i].std():.4f},"
            f"{mut_reps[:,i].mean():.4f},{mut_reps[:,i].std():.4f},"
            f"{delta[i]:.4f},{pvalues[i]:.6f},{pvalues_corrected[i]:.6f},{sig}"
        )
    
    with open(OUTDIR / "delta_rmsf.csv", "w") as f:
        f.write("\n".join(csv_lines) + "\n")
    print(f"Saved: {OUTDIR / 'delta_rmsf.csv'}")
    
    # Save significant regions
    sig_indices = np.where(pvalues_corrected < 0.05)[0]
    with open(OUTDIR / "significant_regions.txt", "w") as f:
        f.write(f"Significant ΔRMSF regions (Bonferroni-corrected p < 0.05)\n")
        f.write(f"Total significant residues: {len(sig_indices)} / {len(common_resids)}\n")
        f.write("=" * 60 + "\n")
        f.write(f"{'Resid':>6} {'ResName':>8} {'WT(Å)':>8} {'4mut(Å)':>8} {'Δ(Å)':>8} {'pval':>12}\n")
        f.write("-" * 60 + "\n")
        for i in sig_indices:
            f.write(f"{int(common_resids[i]):>6} {wt['resnames'][i]:>8} "
                   f"{wt_reps[:,i].mean():>8.2f} {mut_reps[:,i].mean():>8.2f} "
                   f"{delta[i]:>8.2f} {pvalues_corrected[i]:>12.2e}\n")
    print(f"Saved: {OUTDIR / 'significant_regions.txt'}")

if __name__ == "__main__":
    main()
