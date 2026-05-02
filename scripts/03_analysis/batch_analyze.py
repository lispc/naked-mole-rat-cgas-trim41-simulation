#!/usr/bin/env python3
"""
Fast GROMACS batch analysis using streaming + numpy.
Avoids MDAnalysis in_memory=True bottleneck for XTC format.
"""
import numpy as np
import MDAnalysis as mda
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture

outdir = Path("data/analysis/hsap_batch_gmx")
outdir.mkdir(parents=True, exist_ok=True)

replicas = [
    ("Hsap_WT_rep1", "data/md_runs_gmx/Hsap_WT/rep1/analysis.gro", "data/md_runs_gmx/Hsap_WT/rep1/prod.xtc", "WT"),
    ("Hsap_WT_rep2", "data/md_runs_gmx/Hsap_WT/rep2/analysis.gro", "data/md_runs_gmx/Hsap_WT/rep2/prod.xtc", "WT"),
    ("Hsap_4mut_rep1", "data/md_runs_gmx/Hsap_4mut/rep1/analysis.gro", "data/md_runs_gmx/Hsap_4mut/rep1/prod.xtc", "4mut"),
    ("Hsap_4mut_rep2", "data/md_runs_gmx/Hsap_4mut/rep2/analysis.gro", "data/md_runs_gmx/Hsap_4mut/rep2/prod.xtc", "4mut"),
]

def kabsch_align(mobile, reference):
    """Kabsch alignment. mobile, reference are Nx3 arrays."""
    # Center
    mobile_centered = mobile - mobile.mean(axis=0)
    ref_centered = reference - reference.mean(axis=0)
    # Covariance
    H = mobile_centered.T @ ref_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    aligned = (R @ mobile_centered.T).T + reference.mean(axis=0)
    return aligned

def compute_metrics_streaming(u, cgas, trim41, prot_ca, cgas_ca, trim41_ca):
    n_frames = len(u.trajectory)
    n_ca = len(prot_ca)
    
    # Stream all CA coordinates
    print(f"  Streaming {n_frames} frames...")
    coords = np.zeros((n_frames, n_ca, 3))
    com_dists = np.zeros(n_frames)
    rg = np.zeros(n_frames)
    min_dists = np.zeros(n_frames)
    
    for i, ts in enumerate(u.trajectory):
        if i % 5000 == 0:
            print(f"    {i}/{n_frames}")
        coords[i] = prot_ca.positions.copy()
        com_dists[i] = np.linalg.norm(cgas.center_of_mass() - trim41.center_of_mass())
        # Rg
        com = prot_ca.center_of_mass()
        rg[i] = np.sqrt(np.mean(np.sum((prot_ca.positions - com)**2, axis=1)))
        # min CA-CA distance
        min_dists[i] = np.min(np.linalg.norm(
            cgas_ca.positions[:, None] - trim41_ca.positions[None, :], axis=2
        ))
    
    # Align all frames to frame 0 using Kabsch
    print(f"  Aligning {n_frames} frames...")
    ref = coords[0]
    rmsd = np.zeros(n_frames)
    for i in range(n_frames):
        if i % 5000 == 0:
            print(f"    {i}/{n_frames}")
        aligned = kabsch_align(coords[i], ref)
        rmsd[i] = np.sqrt(np.mean(np.sum((aligned - ref)**2, axis=1)))
        coords[i] = aligned
    
    # RMSF
    rmsf = np.sqrt(np.mean(np.sum((coords - coords.mean(axis=0))**2, axis=2), axis=0))
    cgas_mask = np.isin(prot_ca.indices, cgas_ca.indices)
    trim_mask = np.isin(prot_ca.indices, trim41_ca.indices)
    cgas_rmsf = rmsf[cgas_mask]
    trim_rmsf = rmsf[trim_mask]
    
    return rmsd, rmsf, cgas_rmsf, trim_rmsf, com_dists, rg, min_dists

def classify_states(com_dists, min_dists):
    states = np.empty(len(com_dists), dtype=object)
    bound = (com_dists < 40) & (min_dists < 5)
    trans = ((com_dists >= 40) & (com_dists < 48)) | ((min_dists >= 5) & (min_dists < 8))
    unbound = (com_dists >= 48) | (min_dists >= 8)
    states[bound] = "bound"
    states[trans] = "transition"
    states[unbound] = "unbound"
    return states

def main():
    all_data = {}
    all_features = []
    all_names = []
    
    for name, gro, xtc, label in replicas:
        print(f"\n=== {name} ===")
        if not Path(xtc).exists():
            print(f"  Missing {xtc}, skipping")
            continue
        
        u = mda.Universe(gro, xtc)
        print(f"  Frames: {len(u.trajectory)}")
        
        prot = u.select_atoms("protein")
        cgas = u.select_atoms("resindex 0-217")
        trim41 = u.select_atoms("resindex 218-540")
        prot_ca = prot.select_atoms("name CA")
        cgas_ca = cgas.select_atoms("name CA")
        trim41_ca = trim41.select_atoms("name CA")
        
        rmsd, rmsf, cgas_rmsf, trim_rmsf, com_dists, rg, min_dists = compute_metrics_streaming(
            u, cgas, trim41, prot_ca, cgas_ca, trim41_ca
        )
        
        states = classify_states(com_dists, min_dists)
        
        print(f"  RMSD: {rmsd.mean():.2f} ± {rmsd.std():.2f} Å")
        print(f"  COM: {com_dists.mean():.2f} ± {com_dists.std():.2f} Å")
        print(f"  Rg: {rg.mean():.2f} ± {rg.std():.2f} Å")
        print(f"  Min CA: {min_dists.mean():.2f} ± {min_dists.std():.2f} Å")
        print(f"  Bound: {(states=='bound').sum()} ({(states=='bound').mean()*100:.1f}%)")
        print(f"  Transition: {(states=='transition').sum()} ({(states=='transition').mean()*100:.1f}%)")
        print(f"  Unbound: {(states=='unbound').sum()} ({(states=='unbound').mean()*100:.1f}%)")
        
        np.savez(outdir / f"{name}_data.npz",
                 rmsd=rmsd, rmsf=rmsf, cgas_rmsf=cgas_rmsf, trim_rmsf=trim_rmsf,
                 com_dists=com_dists, rg=rg, min_ca_dist=min_dists,
                 states=states)
        
        feats = np.column_stack([com_dists, min_dists, rg, rmsd])
        valid = ~np.isnan(feats).any(axis=1)
        all_features.append(feats[valid])
        all_names.extend([name] * valid.sum())
        all_data[name] = {"label": label, "com": com_dists, "rmsd": rmsd, "rg": rg, "states": states}
    
    if len(all_features) == 0:
        print("No data loaded!")
        return
    
    X = np.vstack(all_features)
    
    # GMM clustering
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    bic = []
    for k in range(2, 8):
        gmm = GaussianMixture(n_components=k, random_state=42, n_init=10)
        gmm.fit(X_scaled)
        bic.append(gmm.bic(X_scaled))
    optimal_k = range(2, 8)[np.argmin(bic)]
    
    gmm = GaussianMixture(n_components=optimal_k, random_state=42, n_init=10)
    cluster_ids = gmm.fit_predict(X_scaled)
    
    print(f"\n=== GMM Clustering (k={optimal_k}) ===")
    idx = 0
    for name, _, _, label in replicas:
        if name not in all_data:
            continue
        n = len(all_data[name]["com"])
        cids = cluster_ids[idx:idx+n]
        print(f"{name}: {dict(zip(*np.unique(cids, return_counts=True)))}")
        idx += n
    
    # Plot comparison
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    
    wt_names = [n for n, _, _, l in replicas if l == "WT" and n in all_data]
    mut_names = [n for n, _, _, l in replicas if l == "4mut" and n in all_data]
    
    for ax, metric in zip(axes.flat, ["rmsd", "com", "rg"]):
        for name in wt_names:
            ax.plot(all_data[name][metric], alpha=0.7, label=name)
        for name in mut_names:
            ax.plot(all_data[name][metric], alpha=0.7, linestyle='--', label=name)
        ax.set_title(metric.upper())
        ax.legend(fontsize=6)
    
    # State occupancy
    ax = axes[1, 2]
    wt_bound = np.mean([np.mean(all_data[n]["states"] == "bound") for n in wt_names])
    wt_unbound = np.mean([np.mean(all_data[n]["states"] == "unbound") for n in wt_names])
    mut_bound = np.mean([np.mean(all_data[n]["states"] == "bound") for n in mut_names])
    mut_unbound = np.mean([np.mean(all_data[n]["states"] == "unbound") for n in mut_names])
    ax.bar([0, 1], [wt_bound, mut_bound], width=0.3, label="bound")
    ax.bar([0.3, 1.3], [wt_unbound, mut_unbound], width=0.3, label="unbound")
    ax.set_xticks([0.15, 1.15])
    ax.set_xticklabels(["WT", "4mut"])
    ax.set_ylabel("Fraction")
    ax.set_title("State Occupancy")
    ax.legend()
    
    fig.tight_layout()
    fig.savefig(outdir / "gmx_comparison.png", dpi=300)
    plt.close(fig)
    print(f"\nSaved: {outdir / 'gmx_comparison.png'}")

if __name__ == "__main__":
    main()
