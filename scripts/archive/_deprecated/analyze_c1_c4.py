#!/usr/bin/env python3
"""
Deep conformational analysis of Cluster C1 vs C4.
C1: most stable bound state (mainly 4mut_rep3)
C4: stable bound state (mainly 4mut_rep2 + WT_rep2)
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

outdir = Path("data/analysis/hsap_batch/c1_c4_analysis")
outdir.mkdir(parents=True, exist_ok=True)

# ---- Load cluster assignments ----
cluster_dir = Path("data/analysis/hsap_batch")

# For each replica, load data + cluster assignment
# We'll reload from the full analysis to get aligned coordinates
# Actually we need to re-read trajectories because we didn't save coordinates

replicas = {
    "Hsap_WT_rep2": ("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd", "WT"),
    "Hsap_4mut_rep2": ("data/md_runs/Hsap_4mut/Hsap_4mut.prmtop", "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd", "4mut"),
    "Hsap_4mut_rep3": ("data/md_runs/Hsap_4mut/Hsap_4mut.prmtop", "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd", "4mut"),
}

# We need cluster assignments. Since we didn't save them per-replica,
# we'll re-run GMM on the saved features just for these 3 replicas.
# Simpler: load npz and recompute cluster ids using saved GMM params.
# Even simpler: the cluster assignments were sequential in the order of replicas.
# Let's reconstruct.

# Load all data again in the same order
all_names = [
    "Hsap_WT_rep1", "Hsap_WT_rep2", "Hsap_WT_rep3",
    "Hsap_4mut_rep1", "Hsap_4mut_rep2", "Hsap_4mut_rep3",
]

feature_names = ["COM", "min_CA_dist", "Rg", "RMSD"] + [f"mut_{m}" for m in ["D431", "K479", "L495", "K498"]]

X_list = []
name_list = []
for name in all_names:
    d = np.load(cluster_dir / f"{name}_data.npz")
    feats = np.column_stack([d["com_dists"], d["min_ca_dist"], d["rg"], d["rmsd"]])
    for mut in ["D431", "K479", "L495", "K498"]:
        dists = np.load(cluster_dir / f"{name}_mut_{mut}_dist.npy")
        feats = np.column_stack([feats, dists])
    valid = ~np.isnan(feats).any(axis=1)
    X_list.append(feats[valid])
    name_list.extend([name] * valid.sum())

X = np.vstack(X_list)
labels = np.array(name_list)

# Re-fit GMM with k=7 to get same cluster ids
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
gmm = GaussianMixture(n_components=7, random_state=42, max_iter=300)
cluster_ids = gmm.fit_predict(X_scaled)

# Map frame indices back to each replica
replica_clusters = {}
idx = 0
for name in all_names:
    n = len(X_list[all_names.index(name)])
    replica_clusters[name] = cluster_ids[idx:idx+n]
    idx += n

print("Cluster distributions for target replicas:")
for name in ["Hsap_WT_rep2", "Hsap_4mut_rep2", "Hsap_4mut_rep3"]:
    cids, counts = np.unique(replica_clusters[name], return_counts=True)
    print(f"  {name}: {dict(zip(cids, counts))}")

# ---- Extract C1 and C4 frames and analyze ----

def analyze_cluster_frames(rep_name, top, traj, system_label, target_cid):
    """Extract frames belonging to target cluster and compute metrics."""
    print(f"\n--- Analyzing {rep_name} C{target_cid} ---")
    
    u = mda.Universe(top, traj)
    ref = mda.Universe(top, traj)
    
    # Align to first frame
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=False).run()
    
    cids = replica_clusters[rep_name]
    mask = cids == target_cid
    frame_indices = np.where(mask)[0]
    n_frames = len(frame_indices)
    
    if n_frames == 0:
        print(f"  No frames in C{target_cid}")
        return None
    
    print(f"  {n_frames} frames")
    
    prot = u.select_atoms("protein")
    cgas = prot.select_atoms("resid 1-218")
    trim = prot.select_atoms("resid 219-541")
    cgas_ca = cgas.select_atoms("name CA")
    trim_ca = trim.select_atoms("name CA")
    
    # Interface residue pairs with min distance < 5A
    interface_pairs = defaultdict(list)
    com_dists = []
    rg_vals = []
    
    # Collect all CA coordinates for average structure
    cgas_ca_coords = []
    trim_ca_coords = []
    
    for fi in frame_indices:
        u.trajectory[fi]
        com_dists.append(np.linalg.norm(cgas.center_of_mass() - trim.center_of_mass()))
        rg_vals.append(prot.radius_of_gyration())
        
        cgas_ca_coords.append(cgas_ca.positions.copy())
        trim_ca_coords.append(trim_ca.positions.copy())
        
        # Compute interface CA-CA distances
        if len(cgas_ca) > 0 and len(trim_ca) > 0:
            dist_mat = np.linalg.norm(cgas_ca.positions[:, None] - trim_ca.positions[None, :], axis=2)
            close_pairs = np.argwhere(dist_mat < 5.0)
            for i, j in close_pairs:
                resid_cgas = cgas_ca.resids[i]
                resid_trim = trim_ca.resids[j]
                pair_key = (int(resid_cgas), int(resid_trim))
                interface_pairs[pair_key].append(dist_mat[i, j])
    
    cgas_avg = np.mean(cgas_ca_coords, axis=0)
    trim_avg = np.mean(trim_ca_coords, axis=0)
    
    # Average interface distance per pair
    avg_interface = {}
    for pair, dists in interface_pairs.items():
        avg_interface[pair] = np.mean(dists)
    
    # Top interface pairs
    top_pairs = sorted(avg_interface.items(), key=lambda x: x[1])[:20]
    
    # Mutation site analysis
    MUT_SITES = {"D431": 19, "K479": 67, "L495": 83, "K498": 86}
    mut_env = {}
    for mut_name, mut_resid in MUT_SITES.items():
        mut_ca = prot.select_atoms(f"resid {mut_resid} and name CA")
        if len(mut_ca) == 0:
            continue
        # Distance to all TRIM41 CA atoms
        dists_to_trim = []
        for fi in frame_indices:
            u.trajectory[fi]
            d = np.linalg.norm(trim_ca.positions - mut_ca.positions[0], axis=1)
            dists_to_trim.append(d)
        dists_to_trim = np.array(dists_to_trim)
        avg_dists = dists_to_trim.mean(axis=0)
        # Nearest 5 TRIM41 residues
        nearest_idx = np.argsort(avg_dists)[:5]
        nearest = [(int(trim_ca.resids[i]), float(avg_dists[i])) for i in nearest_idx]
        mut_env[mut_name] = nearest
    
    return {
        "name": rep_name,
        "system": system_label,
        "cluster": target_cid,
        "n_frames": n_frames,
        "com_mean": float(np.mean(com_dists)),
        "com_std": float(np.std(com_dists)),
        "rg_mean": float(np.mean(rg_vals)),
        "top_interface_pairs": top_pairs,
        "mut_env": mut_env,
        "cgas_avg_ca": cgas_avg,
        "trim_avg_ca": trim_avg,
        "cgas_ca_indices": cgas_ca.indices,
        "trim_ca_indices": trim_ca.indices,
    }

# Analyze C1 (mainly 4mut_rep3) and C4 (mainly 4mut_rep2 + WT_rep2)
c1_data = analyze_cluster_frames("Hsap_4mut_rep3", *replicas["Hsap_4mut_rep3"][:2], replicas["Hsap_4mut_rep3"][2], 1)
c4_4mut = analyze_cluster_frames("Hsap_4mut_rep2", *replicas["Hsap_4mut_rep2"][:2], replicas["Hsap_4mut_rep2"][2], 4)
c4_wt = analyze_cluster_frames("Hsap_WT_rep2", *replicas["Hsap_WT_rep2"][:2], replicas["Hsap_WT_rep2"][2], 4)

# Combine C4 from both systems
c4_combined = {
    "name": "C4_combined",
    "n_frames": c4_4mut["n_frames"] + c4_wt["n_frames"],
    "systems": f"4mut({c4_4mut['n_frames']}) + WT({c4_wt['n_frames']})",
    "com_mean": (c4_4mut["com_mean"] * c4_4mut["n_frames"] + c4_wt["com_mean"] * c4_wt["n_frames"]) / (c4_4mut["n_frames"] + c4_wt["n_frames"]),
    "top_interface_pairs": c4_4mut["top_interface_pairs"],  # Use 4mut as representative
    "mut_env": c4_4mut["mut_env"],  # Use 4mut
}

# ---- Print comparison ----
print("\n" + "="*70)
print("C1 vs C4 COMPARISON")
print("="*70)
print(f"\nC1 ({c1_data['system']} {c1_data['name']}): {c1_data['n_frames']} frames, COM={c1_data['com_mean']:.2f}Å")
print(f"C4 ({c4_combined['systems']}): {c4_combined['n_frames']} frames, COM={c4_combined['com_mean']:.2f}Å")

print(f"\n--- Top Interface Pairs (CA-CA < 5Å) ---")
print(f"{'C1 (4mut_rep3)':<30} | {'C4 (4mut_rep2)':<30}")
print("-"*70)
for i in range(min(15, len(c1_data['top_interface_pairs']), len(c4_4mut['top_interface_pairs']))):
    p1, d1 = c1_data['top_interface_pairs'][i]
    p2, d2 = c4_4mut['top_interface_pairs'][i]
    print(f"cGAS{p1[0]:>3}-TRIM{p1[1]:>3}: {d1:5.2f}Å | cGAS{p2[0]:>3}-TRIM{p2[1]:>3}: {d2:5.2f}Å")

print(f"\n--- Mutation Site Nearest TRIM41 Neighbors ---")
for mut in ["D431", "K479", "L495", "K498"]:
    print(f"\n{mut}:")
    print(f"  C1 nearest: {c1_data['mut_env'].get(mut, [])}")
    print(f"  C4 nearest: {c4_4mut['mut_env'].get(mut, [])}")

# ---- Visualization: Interface pair distance comparison ----
fig, ax = plt.subplots(figsize=(10, 6))

# Get common pairs between C1 and C4
c1_dict = dict(c1_data['top_interface_pairs'])
c4_dict = dict(c4_4mut['top_interface_pairs'])
common_pairs = sorted(set(c1_dict.keys()) & set(c4_dict.keys()), key=lambda p: c1_dict[p])

if len(common_pairs) > 0:
    c1_vals = [c1_dict[p] for p in common_pairs]
    c4_vals = [c4_dict[p] for p in common_pairs]
    labels = [f"c{p[0]}-t{p[1]}" for p in common_pairs]
    
    x = np.arange(len(common_pairs))
    width = 0.35
    ax.bar(x - width/2, c1_vals, width, label="C1 (4mut_rep3)", color="tab:blue")
    ax.bar(x + width/2, c4_vals, width, label="C4 (4mut_rep2)", color="tab:orange")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_ylabel("Avg CA-CA distance (Å)")
    ax.set_title("Interface Residue Pair Distances: C1 vs C4")
    ax.legend()
    ax.axhline(5.0, color="gray", ls="--", alpha=0.5, label="Contact cutoff")
    fig.tight_layout()
    fig.savefig(outdir / "c1_c4_interface_pairs.png", dpi=300)
    plt.close(fig)
    print(f"\n  Saved interface pair comparison ({len(common_pairs)} common pairs)")

# ---- Visualization: Mutation site distance distributions ----
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()

# We need to recompute distance distributions for mutation sites
MUT_SITES = {"D431": 19, "K479": 67, "L495": 83, "K498": 86}

for idx, (mut_name, mut_resid) in enumerate(MUT_SITES.items()):
    ax = axes[idx]
    
    # C1: 4mut_rep3
    u1 = mda.Universe(replicas["Hsap_4mut_rep3"][0], replicas["Hsap_4mut_rep3"][1])
    mut_ca1 = u1.select_atoms(f"resid {mut_resid} and name CA")
    trim_ca1 = u1.select_atoms("resid 219-541 and name CA")
    c1_frames = np.where(replica_clusters["Hsap_4mut_rep3"] == 1)[0]
    c1_dists = []
    for fi in c1_frames:
        u1.trajectory[fi]
        c1_dists.append(np.min(np.linalg.norm(trim_ca1.positions - mut_ca1.positions[0], axis=1)))
    
    # C4: 4mut_rep2
    u2 = mda.Universe(replicas["Hsap_4mut_rep2"][0], replicas["Hsap_4mut_rep2"][1])
    mut_ca2 = u2.select_atoms(f"resid {mut_resid} and name CA")
    trim_ca2 = u2.select_atoms("resid 219-541 and name CA")
    c4_frames = np.where(replica_clusters["Hsap_4mut_rep2"] == 4)[0]
    c4_dists = []
    for fi in c4_frames:
        u2.trajectory[fi]
        c4_dists.append(np.min(np.linalg.norm(trim_ca2.positions - mut_ca2.positions[0], axis=1)))
    
    ax.hist(c1_dists, bins=30, alpha=0.5, label=f"C1 (n={len(c1_dists)})", color="tab:blue", density=True)
    ax.hist(c4_dists, bins=30, alpha=0.5, label=f"C4 (n={len(c4_dists)})", color="tab:orange", density=True)
    ax.axvline(np.mean(c1_dists), color="tab:blue", ls="--")
    ax.axvline(np.mean(c4_dists), color="tab:orange", ls="--")
    ax.set_xlabel("Min distance to TRIM41 CA (Å)")
    ax.set_ylabel("Density")
    ax.set_title(f"{mut_name} (cGAS res {mut_resid})")
    ax.legend()

fig.tight_layout()
fig.savefig(outdir / "c1_c4_mutation_distances.png", dpi=300)
plt.close(fig)

# ---- RMSF within C1 vs within C4 ----
print("\nComputing RMSF within each cluster...")

def compute_cluster_rmsf(rep_name, top, traj, target_cid, system_label):
    u = mda.Universe(top, traj)
    ref = mda.Universe(top, traj)
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=False).run()
    
    ca = u.select_atoms("protein and name CA")
    cids = replica_clusters[rep_name]
    frame_indices = np.where(cids == target_cid)[0]
    
    # Collect coordinates
    coords = []
    for fi in frame_indices:
        u.trajectory[fi]
        coords.append(ca.positions.copy())
    coords = np.array(coords)
    
    mean_pos = coords.mean(axis=0)
    rmsf = np.sqrt(((coords - mean_pos) ** 2).mean(axis=0).sum(axis=1))
    return ca.resids, rmsf

resids_c1, rmsf_c1 = compute_cluster_rmsf("Hsap_4mut_rep3", *replicas["Hsap_4mut_rep3"][:2], 1, "4mut")
resids_c4, rmsf_c4 = compute_cluster_rmsf("Hsap_4mut_rep2", *replicas["Hsap_4mut_rep2"][:2], 4, "4mut")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(resids_c1, rmsf_c1, label="C1 (4mut_rep3)", color="tab:blue", lw=1)
ax.plot(resids_c4, rmsf_c4, label="C4 (4mut_rep2)", color="tab:orange", lw=1)
ax.axvline(218.5, color="gray", ls="--", alpha=0.5)
ax.set_xlabel("Residue ID")
ax.set_ylabel("CA RMSF (Å)")
ax.set_title("RMSF within C1 vs within C4 (aligned to cluster mean)")
ax.legend()
fig.tight_layout()
fig.savefig(outdir / "c1_c4_rmsf.png", dpi=300)
plt.close(fig)

# ---- Summary output ----
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"C1 (4mut_rep3): {c1_data['n_frames']} frames, COM={c1_data['com_mean']:.2f}±{c1_data['com_std']:.2f}Å")
print(f"C4 (4mut_rep2): {c4_4mut['n_frames']} frames, COM={c4_4mut['com_mean']:.2f}±{c4_4mut['com_std']:.2f}Å")
print(f"\nOutputs saved to {outdir}")
print("  c1_c4_interface_pairs.png")
print("  c1_c4_mutation_distances.png")
print("  c1_c4_rmsf.png")
