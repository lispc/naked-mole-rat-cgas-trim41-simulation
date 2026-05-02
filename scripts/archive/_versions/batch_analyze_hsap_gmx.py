#!/usr/bin/env python3
"""
Batch analysis of Hsap GROMACS replicas (WT + 4mut).
Adapts batch_analyze_hsap.py for GROMACS gro+xtc format.
Uses resindex instead of resid for cGAS/TRIM41 selection due to GROMACS resid reset.
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import distance_array
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture

outdir = Path("data/analysis/hsap_batch_gmx")
outdir.mkdir(parents=True, exist_ok=True)

replicas = [
    ("Hsap_WT_rep1", "data/md_runs_gmx/Hsap_WT/rep1/Hsap_WT_rep1.gro", "data/md_runs_gmx/Hsap_WT/rep1/prod.xtc", "WT"),
    ("Hsap_WT_rep2", "data/md_runs_gmx/Hsap_WT/rep2/Hsap_WT_rep2.gro", "data/md_runs_gmx/Hsap_WT/rep2/prod.xtc", "WT"),
    ("Hsap_4mut_rep1", "data/md_runs_gmx/Hsap_4mut/rep1/Hsap_4mut_rep1.gro", "data/md_runs_gmx/Hsap_4mut/rep1/prod.xtc", "4mut"),
    ("Hsap_4mut_rep2", "data/md_runs_gmx/Hsap_4mut/rep2/Hsap_4mut_rep2.gro", "data/md_runs_gmx/Hsap_4mut/rep2/prod.xtc", "4mut"),
]

def compute_metrics(u, cgas, trim41, prot_ca, cgas_ca, trim41_ca):
    # Align
    aligner = align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    
    # RMSD
    rmsd = rms.RMSD(prot_ca, prot_ca, ref_frame=0).run().results.rmsd[:, 2]
    
    # RMSF
    rmsf = rms.RMSF(prot_ca).run().results.rmsf
    cgas_rmsf = rmsf[np.isin(prot_ca.indices, cgas_ca.indices)]
    trim_rmsf = rmsf[np.isin(prot_ca.indices, trim41_ca.indices)]
    
    # COM distance
    com_dists = np.array([np.linalg.norm(cgas.center_of_mass() - trim41.center_of_mass()) for ts in u.trajectory])
    
    # Rg
    rg = np.array([prot_ca.radius_of_gyration() for ts in u.trajectory])
    
    # Interface min CA-CA distance
    min_dists = np.array([distance_array(cgas_ca.positions, trim41_ca.positions).min() for ts in u.trajectory])
    
    # Mutation sites (resindex-based mapping)
    # In AMBER: D431->res19, K479->res67, L495->res83, K498->res86
    # These are cGAS residue indices within the cGAS chain (0-based)
    mut_cgas_indices = [18, 66, 82, 85]  # 0-based within cGAS
    mut_dists = []
    for ts in u.trajectory:
        frame_dists = []
        for idx in mut_cgas_indices:
            mut_ca = cgas.residues[idx].atoms.select_atoms("name CA")
            if len(mut_ca) == 0:
                frame_dists.append(np.nan)
            else:
                d = distance_array(mut_ca.positions, trim41_ca.positions).min()
                frame_dists.append(d)
        mut_dists.append(frame_dists)
    mut_dists = np.array(mut_dists)
    
    return rmsd, rmsf, cgas_rmsf, trim_rmsf, com_dists, rg, min_dists, mut_dists

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
        
        # Use resindex for selection (GROMACS resets resid per molecule)
        prot = u.select_atoms("protein")
        cgas = u.select_atoms("resindex 0-217")
        trim41 = u.select_atoms("resindex 218-540")
        prot_ca = prot.select_atoms("name CA")
        cgas_ca = cgas.select_atoms("name CA")
        trim41_ca = trim41.select_atoms("name CA")
        
        print(f"  Protein CA: {len(prot_ca)}, cGAS CA: {len(cgas_ca)}, TRIM41 CA: {len(trim41_ca)}")
        
        rmsd, rmsf, cgas_rmsf, trim_rmsf, com_dists, rg, min_dists, mut_dists = compute_metrics(
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
                 com_dists=com_dists, rg=rg, min_ca_dist=min_dists, mut_dists=mut_dists,
                 states=states)
        
        feats = np.column_stack([com_dists, min_dists, rg, rmsd, mut_dists])
        valid = ~np.isnan(feats).any(axis=1)
        all_features.append(feats[valid])
        all_names.extend([name] * valid.sum())
        all_data[name] = {"label": label, "com": com_dists, "rmsd": rmsd, "rg": rg, "states": states}
    
    if len(all_features) == 0:
        print("No data loaded!")
        return
    
    X = np.vstack(all_features)
    labels = np.array(all_names)
    
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
