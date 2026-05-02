#!/usr/bin/env python3
"""
Cluster analysis of all 6 Hsap replicas to identify metastable states.
Uses features extracted from pre-computed npz files.
"""
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import json

outdir = Path("data/analysis/hsap_batch")

replicas = [
    "Hsap_WT_rep1", "Hsap_WT_rep2", "Hsap_WT_rep3",
    "Hsap_4mut_rep1", "Hsap_4mut_rep2", "Hsap_4mut_rep3",
]

# ---- Load pre-computed features ----
print("Loading features from npz files...")
all_data = {}
for name in replicas:
    fpath = outdir / f"{name}_data.npz"
    if not fpath.exists():
        print(f"  WARNING: {fpath} not found, skipping")
        continue
    d = np.load(fpath)
    all_data[name] = {
        "time_ns": d["time_ns"],
        "rmsd": d["rmsd"],
        "com_dists": d["com_dists"],
        "rg": d["rg"],
        "min_ca_dist": d["min_ca_dist"],
    }
    print(f"  {name}: {len(d['time_ns'])} frames")

# ---- Load mutation site distances ----
MUTATION_SITES = ["D431", "K479", "L495", "K498"]
for name in all_data:
    for mut in MUTATION_SITES:
        fpath = outdir / f"{name}_mut_{mut}_dist.npy"
        if fpath.exists():
            all_data[name][f"mut_{mut}"] = np.load(fpath)

# ---- Build feature matrix ----
# For each frame, extract: [COM, min_ca_dist, Rg, RMSD, mut_D431, mut_K479, mut_L495, mut_K498]
feature_names = ["COM", "min_CA_dist", "Rg", "RMSD"] + [f"mut_{m}" for m in MUTATION_SITES]

X_list = []
labels_list = []
rep_labels = []

for name in all_data:
    d = all_data[name]
    n = len(d["time_ns"])
    feats = np.column_stack([
        d["com_dists"],
        d["min_ca_dist"],
        d["rg"],
        d["rmsd"],
    ])
    for mut in MUTATION_SITES:
        key = f"mut_{mut}"
        if key in d:
            feats = np.column_stack([feats, d[key]])
        else:
            feats = np.column_stack([feats, np.full(n, np.nan)])
    
    # Only keep frames where all features are valid
    valid = ~np.isnan(feats).any(axis=1)
    feats = feats[valid]
    
    X_list.append(feats)
    labels_list.extend([name] * len(feats))
    rep_labels.append({"name": name, "n_frames": len(feats), "valid": valid})

X = np.vstack(X_list)
labels = np.array(labels_list)
print(f"\nTotal frames for clustering: {len(X)}")
print(f"Features: {feature_names}")

# ---- Standardize features ----
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ---- Determine optimal number of clusters ----
print("\nEvaluating cluster numbers...")
bics = []
ks = range(2, 8)
for k in ks:
    gmm = GaussianMixture(n_components=k, random_state=42, max_iter=200)
    gmm.fit(X_scaled)
    bics.append(gmm.bic(X_scaled))
    print(f"  k={k}: BIC={gmm.bic(X_scaled):.1f}")

optimal_k = ks[np.argmin(bics)]
print(f"\nOptimal clusters by BIC: {optimal_k}")

# ---- Fit GMM with optimal k ----
print(f"\nFitting GMM (k={optimal_k})...")
gmm = GaussianMixture(n_components=optimal_k, random_state=42, max_iter=300)
cluster_ids = gmm.fit_predict(X_scaled)
probs = gmm.predict_proba(X_scaled)

# ---- Characterize each cluster ----
print("\n" + "="*70)
print("CLUSTER CHARACTERIZATION (raw feature means)")
print("="*70)
cluster_summary = {}
for cid in range(optimal_k):
    mask = cluster_ids == cid
    n = mask.sum()
    print(f"\nCluster {cid}: {n} frames ({n/len(X)*100:.1f}%)")
    means = X[mask].mean(axis=0)
    stds = X[mask].std(axis=0)
    for i, fname in enumerate(feature_names):
        print(f"  {fname:<15}: {means[i]:>8.2f} ± {stds[i]:>6.2f}")
    cluster_summary[cid] = {
        "n_frames": int(n),
        "fraction": float(n/len(X)),
        "means": {fname: float(means[i]) for i, fname in enumerate(feature_names)},
        "stds": {fname: float(stds[i]) for i, fname in enumerate(feature_names)},
    }

# ---- Cluster occupancy per replica ----
print("\n" + "="*70)
print("CLUSTER OCCUPANCY PER REPLICA (%)")
print("="*70)
header = f"{'Replica':<18}" + "".join([f"C{i:>8}" for i in range(optimal_k)])
print(header)
print("-"*70)

occupancy = {}
for name in all_data:
    mask = labels == name
    occ = {}
    for cid in range(optimal_k):
        frac = np.sum((cluster_ids[mask] == cid)) / mask.sum() * 100
        occ[cid] = float(frac)
    occupancy[name] = occ
    row = f"{name:<18}" + "".join([f"{occ[cid]:>7.1f}%" for cid in range(optimal_k)])
    print(row)

# ---- Transition matrix (pooled across replicas) ----
print("\n" + "="*70)
print("TRANSITION COUNTS (pooled)")
print("="*70)
# Count transitions between clusters
from collections import Counter
transitions = Counter()
for i in range(len(cluster_ids) - 1):
    if labels[i] == labels[i+1]:  # only within same replica
        transitions[(cluster_ids[i], cluster_ids[i+1])] += 1

print("From -> To : count")
for (c1, c2), count in sorted(transitions.items()):
    print(f"  {c1} -> {c2} : {count}")

# ---- Visualizations ----
print("\nGenerating plots...")

# 1. BIC plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(ks, bics, "o-")
ax.axvline(optimal_k, color="red", ls="--", alpha=0.5, label=f"Optimal: {optimal_k}")
ax.set_xlabel("Number of clusters")
ax.set_ylabel("BIC")
ax.set_title("GMM Model Selection")
ax.legend()
fig.tight_layout()
fig.savefig(outdir / "cluster_bic.png", dpi=300)
plt.close(fig)

# 2. PCA projection colored by cluster
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# By cluster
ax = axes[0]
colors = plt.cm.tab10(np.linspace(0, 1, optimal_k))
for cid in range(optimal_k):
    mask = cluster_ids == cid
    ax.scatter(X_pca[mask, 0], X_pca[mask, 1], c=[colors[cid]], s=1, alpha=0.5, label=f"C{cid}")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA projection colored by cluster")
ax.legend(markerscale=5)

# By system (WT vs 4mut)
ax = axes[1]
wt_mask = np.array(["WT" in l for l in labels])
mut_mask = np.array(["4mut" in l for l in labels])
ax.scatter(X_pca[wt_mask, 0], X_pca[wt_mask, 1], c="tab:blue", s=1, alpha=0.4, label="WT")
ax.scatter(X_pca[mut_mask, 0], X_pca[mut_mask, 1], c="tab:orange", s=1, alpha=0.4, label="4mut")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA projection colored by system")
ax.legend()

fig.tight_layout()
fig.savefig(outdir / "cluster_pca.png", dpi=300)
plt.close(fig)

# 3. Feature distributions per cluster
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()
for i, fname in enumerate(feature_names[:6]):  # first 6 features
    ax = axes[i]
    for cid in range(optimal_k):
        mask = cluster_ids == cid
        ax.hist(X[mask, i], bins=50, alpha=0.4, label=f"C{cid}", color=colors[cid])
    ax.set_xlabel(fname)
    ax.set_ylabel("Count")
    ax.set_title(f"Distribution: {fname}")
    ax.legend()
fig.tight_layout()
fig.savefig(outdir / "cluster_feature_distributions.png", dpi=300)
plt.close(fig)

# 4. Cluster occupancy bar chart
fig, ax = plt.subplots(figsize=(10, 5))
x = np.arange(len(all_data))
width = 0.8 / optimal_k
for cid in range(optimal_k):
    vals = [occupancy[name][cid] for name in all_data]
    ax.bar(x + cid * width, vals, width, label=f"C{cid}", color=colors[cid])
ax.set_xticks(x + width * (optimal_k - 1) / 2)
ax.set_xticklabels(all_data, rotation=45, ha="right")
ax.set_ylabel("Occupancy (%)")
ax.set_title("Cluster Occupancy per Replica")
ax.legend()
fig.tight_layout()
fig.savefig(outdir / "cluster_occupancy.png", dpi=300)
plt.close(fig)

# 5. Time evolution of cluster assignment per replica
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()
for idx, name in enumerate(all_data):
    ax = axes[idx]
    mask = labels == name
    time = all_data[name]["time_ns"]
    # time may have NaN-filtered frames removed; align by index
    valid = rep_labels[idx]["valid"]
    time_valid = time[valid]
    clusters = cluster_ids[mask]
    
    ax.scatter(time_valid, clusters, c=[colors[c] for c in clusters], s=2, alpha=0.6)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Cluster")
    ax.set_title(name)
    ax.set_yticks(range(optimal_k))
    ax.set_ylim(-0.5, optimal_k - 0.5)
fig.tight_layout()
fig.savefig(outdir / "cluster_time_evolution.png", dpi=300)
plt.close(fig)

# ---- Save results ----
results = {
    "optimal_k": int(optimal_k),
    "feature_names": feature_names,
    "cluster_summary": cluster_summary,
    "occupancy": occupancy,
    "pca_explained_variance": pca.explained_variance_ratio_.tolist(),
}
with open(outdir / "cluster_results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\n✅ Cluster analysis complete. Outputs in {outdir}")
print(f"  cluster_bic.png")
print(f"  cluster_pca.png")
print(f"  cluster_feature_distributions.png")
print(f"  cluster_occupancy.png")
print(f"  cluster_time_evolution.png")
print(f"  cluster_results.json")
