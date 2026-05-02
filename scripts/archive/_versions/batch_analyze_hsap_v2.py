#!/usr/bin/env python3
"""
Lightweight batch analysis of all 6 Hsap replicas.
Optimized: uses CA-only for interface metrics, avoids heavy O(N^2) loops.
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
import time

outdir = Path("data/analysis/hsap_batch")
outdir.mkdir(parents=True, exist_ok=True)

systems = {
    "Hsap_WT_rep1": ("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd"),
    "Hsap_WT_rep2": ("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd"),
    "Hsap_WT_rep3": ("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd"),
    "Hsap_4mut_rep1": ("data/md_runs/Hsap_4mut/Hsap_4mut.prmtop", "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd"),
    "Hsap_4mut_rep2": ("data/md_runs/Hsap_4mut/Hsap_4mut.prmtop", "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd"),
    "Hsap_4mut_rep3": ("data/md_runs/Hsap_4mut/Hsap_4mut.prmtop", "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd"),
}

MUTATION_SITES = {"D431": 19, "K479": 67, "L495": 83, "K498": 86}

def analyze_one(name, top, traj):
    t0 = time.time()
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    u = mda.Universe(top, traj)
    ref = mda.Universe(top, traj)
    
    prot = u.select_atoms("protein")
    cgas = prot.select_atoms("resid 1-218")
    trim = prot.select_atoms("resid 219-541")
    ca_atoms = prot.select_atoms("name CA")
    cgas_ca = prot.select_atoms("name CA and resid 1-218")
    trim_ca = prot.select_atoms("name CA and resid 219-541")
    
    n_frames = len(u.trajectory)
    dt_ns = u.trajectory.dt / 1000.0
    time_ns = np.arange(n_frames) * dt_ns
    
    print(f"  {n_frames} frames, loading+align...")
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=False).run()
    
    print(f"  RMSD...")
    rmsd_vals = rms.RMSD(u, ref, select="protein and name CA").run().results.rmsd[:, 2]
    
    print(f"  RMSF...")
    rmsf_vals = rms.RMSF(ca_atoms).run().results.rmsf
    ca_resids = ca_atoms.resids
    
    print(f"  COM & Rg...")
    com_dists = np.empty(n_frames)
    rg_vals = np.empty(n_frames)
    for i, ts in enumerate(u.trajectory):
        com_dists[i] = np.linalg.norm(cgas.center_of_mass() - trim.center_of_mass())
        rg_vals[i] = prot.radius_of_gyration()
    
    # Fast interface contact metric: min CA-CA distance between cGAS and TRIM41
    print(f"  Interface CA contacts...")
    min_ca_dist = np.empty(n_frames)
    for i, ts in enumerate(u.trajectory):
        dmat = np.linalg.norm(cgas_ca.positions[:, None] - trim_ca.positions[None, :], axis=2)
        min_ca_dist[i] = dmat.min()
    
    # Mutation site distances (CA to nearest TRIM41 CA)
    print(f"  Mutation sites...")
    mut_distances = {}
    for mut_name, mut_resid in MUTATION_SITES.items():
        mut_ca = prot.select_atoms(f"resid {mut_resid} and name CA")
        if len(mut_ca) == 0:
            mut_distances[mut_name] = np.full(n_frames, np.nan)
            continue
        dists = np.empty(n_frames)
        for i, ts in enumerate(u.trajectory):
            dists[i] = np.min(np.linalg.norm(trim_ca.positions - mut_ca.positions[0], axis=1))
        mut_distances[mut_name] = dists
    
    # State classification
    states = np.full(n_frames, "", dtype=object)
    states[com_dists < 40] = "bound"
    states[(com_dists >= 40) & (com_dists < 48)] = "transition"
    states[com_dists >= 48] = "unbound"
    
    state_counts = {
        "bound": int(np.sum(states == "bound")),
        "transition": int(np.sum(states == "transition")),
        "unbound": int(np.sum(states == "unbound")),
    }
    unbound_frames = np.where(states == "unbound")[0]
    time_to_unbind = float(time_ns[unbound_frames[0]]) if len(unbound_frames) > 0 else None
    
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s | RMSD={rmsd_vals.mean():.2f} | COM={com_dists.mean():.2f} | Unbind={time_to_unbind}")
    
    return {
        "name": name,
        "time_ns": time_ns,
        "rmsd": rmsd_vals,
        "rmsf": rmsf_vals,
        "ca_resids": ca_resids,
        "com_dists": com_dists,
        "rg": rg_vals,
        "min_ca_dist": min_ca_dist,
        "mut_distances": mut_distances,
        "states": state_counts,
        "time_to_unbind": time_to_unbind,
        "n_frames": n_frames,
    }

all_results = []
for name, (top, traj) in systems.items():
    try:
        result = analyze_one(name, top, traj)
        all_results.append(result)
    except Exception as e:
        print(f"ERROR: {name}: {e}")
        import traceback; traceback.print_exc()

# Save data
print("\nSaving data...")
for r in all_results:
    np.savez(outdir / f"{r['name']}_data.npz",
        time_ns=r["time_ns"], rmsd=r["rmsd"], rmsf=r["rmsf"],
        ca_resids=r["ca_resids"], com_dists=r["com_dists"],
        rg=r["rg"], min_ca_dist=r["min_ca_dist"])
    for mut_name, dists in r["mut_distances"].items():
        np.save(outdir / f"{r['name']}_mut_{mut_name}_dist.npy", dists)

summary = {}
for r in all_results:
    summary[r["name"]] = {
        "rmsd_mean": float(r["rmsd"].mean()),
        "rmsd_std": float(r["rmsd"].std()),
        "com_mean": float(r["com_dists"].mean()),
        "com_std": float(r["com_dists"].std()),
        "com_min": float(r["com_dists"].min()),
        "com_max": float(r["com_dists"].max()),
        "rg_mean": float(r["rg"].mean()),
        "rg_std": float(r["rg"].std()),
        "min_ca_dist_mean": float(r["min_ca_dist"].mean()),
        "states": r["states"],
        "time_to_unbind_ns": r["time_to_unbind"],
    }

with open(outdir / "summary.json", "w") as f:
    json.dump(summary, f, indent=2)

# Plots
print("Generating plots...")

fig, axes = plt.subplots(3, 2, figsize=(12, 14))
for r in all_results:
    color = "tab:blue" if "WT" in r["name"] else "tab:orange"
    axes[0,0].plot(r["time_ns"], r["rmsd"], color=color, alpha=0.6, lw=0.5)
    axes[0,1].plot(r["time_ns"], r["com_dists"], color=color, alpha=0.6, lw=0.5)
    axes[1,0].plot(r["time_ns"], r["rg"], color=color, alpha=0.6, lw=0.5)
    axes[1,1].plot(r["time_ns"], r["min_ca_dist"], color=color, alpha=0.6, lw=0.5)
    for mut_name, dists in r["mut_distances"].items():
        axes[2,0].plot(r["time_ns"], dists, color=color, alpha=0.3, lw=0.4)

axes[0,0].set_title("CA RMSD")
axes[0,0].set_ylabel("RMSD (Å)")
axes[0,1].set_title("COM Distance")
axes[0,1].set_ylabel("Å")
axes[1,0].set_title("Radius of Gyration")
axes[1,0].set_ylabel("Å")
axes[1,1].set_title("Min CA-CA Interface Distance")
axes[1,1].set_ylabel("Å")
axes[2,0].set_title("Mutation Site Distances")
axes[2,0].set_ylabel("Å")
axes[2,1].axis("off")

for ax in axes.flat:
    if ax.has_data():
        ax.set_xlabel("Time (ns)")

from matplotlib.lines import Line2D
fig.legend([Line2D([0],[0],color="tab:blue",lw=2), Line2D([0],[0],color="tab:orange",lw=2)],
           ["WT (3 reps)", "4mut (3 reps)"], loc="upper center", ncol=2, bbox_to_anchor=(0.5,0.98))
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig(outdir / "all_replicas_comparison.png", dpi=300)
plt.close(fig)

# State occupancy
fig, ax = plt.subplots(figsize=(10,4))
names = [r["name"] for r in all_results]
bound = [r["states"]["bound"]/r["n_frames"]*100 for r in all_results]
trans = [r["states"]["transition"]/r["n_frames"]*100 for r in all_results]
unb = [r["states"]["unbound"]/r["n_frames"]*100 for r in all_results]
x = np.arange(len(names))
ax.bar(x-0.25, bound, 0.25, label="Bound (<40Å)", color="tab:green")
ax.bar(x, trans, 0.25, label="Transition", color="tab:orange")
ax.bar(x+0.25, unb, 0.25, label="Unbound (>48Å)", color="tab:red")
ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_ylabel("Occupancy (%)"); ax.set_title("State Occupancy"); ax.legend()
fig.tight_layout(); fig.savefig(outdir / "state_occupancy.png", dpi=300); plt.close(fig)

# Time to unbind
fig, ax = plt.subplots(figsize=(8,4))
ttb = [r["time_to_unbind"] if r["time_to_unbind"] else 200 for r in all_results]
colors = ["tab:blue" if "WT" in n else "tab:orange" for n in names]
ax.bar(names, ttb, color=colors, alpha=0.7)
ax.set_ylabel("Time to unbind (ns)"); ax.set_title("First Entry into Unbound State")
ax.set_xticklabels(names, rotation=45, ha="right")
fig.tight_layout(); fig.savefig(outdir / "time_to_unbind.png", dpi=300); plt.close(fig)

# Average RMSF
fig, ax = plt.subplots(figsize=(10,4))
wt_rmsf = [r["rmsf"] for r in all_results if "WT" in r["name"]]
mut_rmsf = [r["rmsf"] for r in all_results if "4mut" in r["name"]]
if len(wt_rmsf)==len(mut_rmsf)==3:
    wt_m, wt_s = np.mean(wt_rmsf,axis=0), np.std(wt_rmsf,axis=0)
    mut_m, mut_s = np.mean(mut_rmsf,axis=0), np.std(mut_rmsf,axis=0)
    resids = all_results[0]["ca_resids"]
    ax.plot(resids, wt_m, color="tab:blue", label="WT mean±std")
    ax.fill_between(resids, wt_m-wt_s, wt_m+wt_s, color="tab:blue", alpha=0.2)
    ax.plot(resids, mut_m, color="tab:orange", label="4mut mean±std")
    ax.fill_between(resids, mut_m-mut_s, mut_m+mut_s, color="tab:orange", alpha=0.2)
    ax.axvline(218.5, color="gray", ls="--", alpha=0.5)
    ax.set_xlabel("Residue"); ax.set_ylabel("CA RMSF (Å)")
    ax.set_title("Average RMSF"); ax.legend()
    fig.tight_layout(); fig.savefig(outdir / "rmsf_avg_comparison.png", dpi=300); plt.close(fig)

# Summary table
print("\n" + "="*80)
print("SUMMARY: All 6 Replicas")
print("="*80)
print(f"{'Replica':<18} {'RMSD':>8} {'COM':>8} {'Rg':>7} {'MinCA':>7} {'Unbind':>8}")
print("-"*80)
for r in all_results:
    s = summary[r["name"]]
    unb = f"{s['time_to_unbind_ns']:.1f}ns" if s['time_to_unbind_ns'] else "N/A"
    print(f"{r['name']:<18} {s['rmsd_mean']:>7.2f} {s['com_mean']:>7.2f} {s['rg_mean']:>6.2f} {s['min_ca_dist_mean']:>6.2f} {unb:>8}")

print(f"\n✅ All outputs in {outdir}")
