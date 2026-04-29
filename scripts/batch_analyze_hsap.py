#!/usr/bin/env python3
"""
Batch analysis of all 6 Hsap replicas (WT rep1-3 + 4mut rep1-3).
Outputs: per-replica metrics, cross-replica comparison, dissociation kinetics,
interface breaking analysis, and mutation site distance tracking.
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json

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

# ---- Mutation site mapping (prmtop residue numbering) ----
# Hsap cGAS C-terminal residues in prmtop: 1-218
# Original cGAS numbering: 413-630
# D431 -> prmtop resid = 431 - 412 = 19
# K479 -> prmtop resid = 479 - 412 = 67
# L495 -> prmtop resid = 495 - 412 = 83
# K498 -> prmtop resid = 498 - 412 = 86
MUTATION_SITES = {
    "D431": 19,
    "K479": 67,
    "L495": 83,
    "K498": 86,
}

def analyze_one(name, top, traj):
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")
    
    u = mda.Universe(top, traj)
    ref = mda.Universe(top, traj)
    
    prot = u.select_atoms("protein")
    cgas = prot.select_atoms("resid 1-218")
    trim = prot.select_atoms("resid 219-541")
    ca_atoms = prot.select_atoms("name CA")
    cgas_ca = prot.select_atoms("name CA and resid 1-218")
    trim_ca = prot.select_atoms("name CA and resid 219-541")
    
    n_frames = len(u.trajectory)
    dt_ns = u.trajectory.dt / 1000.0  # ps -> ns
    time_ns = np.arange(n_frames) * dt_ns
    
    print(f"  Frames: {n_frames}, dt: {dt_ns:.3f} ns, total: {time_ns[-1]:.1f} ns")
    print(f"  Protein atoms: {len(prot)}, cGAS res: {len(cgas.residues)}, TRIM41 res: {len(trim.residues)}")
    
    # Align and RMSD
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    rmsd_analysis = rms.RMSD(u, ref, select="protein and name CA").run()
    rmsd_vals = rmsd_analysis.results.rmsd[:, 2]
    
    # RMSF
    rmsf_analysis = rms.RMSF(ca_atoms).run()
    rmsf_vals = rmsf_analysis.results.rmsf
    ca_resids = ca_atoms.resids
    
    # COM distance
    com_dists = []
    for ts in u.trajectory:
        com_cgas = cgas.center_of_mass()
        com_trim = trim.center_of_mass()
        com_dists.append(np.linalg.norm(com_cgas - com_trim))
    com_dists = np.array(com_dists)
    
    # Radius of gyration
    rg_vals = []
    for ts in u.trajectory:
        rg_vals.append(prot.radius_of_gyration())
    rg_vals = np.array(rg_vals)
    
    # Interface contacts (< 5A between any cGAS-TRIM41 heavy atoms)
    # Sample every frame but use a faster method
    contact_counts = []
    for ts in u.trajectory:
        cgas_pos = cgas.positions
        trim_pos = trim.positions
        dist_mat = np.linalg.norm(cgas_pos[:, np.newaxis] - trim_pos[np.newaxis, :], axis=2)
        n_contacts = np.sum(dist_mat < 5.0)
        contact_counts.append(n_contacts)
    contact_counts = np.array(contact_counts)
    
    # Interface H-bonds (N/O atoms within 3.5A)
    hbond_counts = []
    cgas_NO = cgas.select_atoms("name N* or name O*")
    trim_NO = trim.select_atoms("name N* or name O*")
    for ts in u.trajectory:
        if len(cgas_NO) > 0 and len(trim_NO) > 0:
            dist_mat = np.linalg.norm(cgas_NO.positions[:, np.newaxis] - trim_NO.positions[np.newaxis, :], axis=2)
            n_hb = np.sum(dist_mat < 3.5)
            hbond_counts.append(n_hb)
        else:
            hbond_counts.append(0)
    hbond_counts = np.array(hbond_counts)
    
    # Mutation site distances to nearest TRIM41 atom
    mut_distances = {}
    for mut_name, mut_resid in MUTATION_SITES.items():
        mut_atom = prot.select_atoms(f"resid {mut_resid} and name CA")
        if len(mut_atom) == 0:
            mut_distances[mut_name] = np.full(n_frames, np.nan)
            continue
        dists = []
        for ts in u.trajectory:
            d = np.min(np.linalg.norm(trim.positions - mut_atom.positions[0], axis=1))
            dists.append(d)
        mut_distances[mut_name] = np.array(dists)
    
    # State classification based on COM distance
    # Bound: COM < 40A, Transition: 40-48A, Unbound: > 48A
    states = np.full(n_frames, "", dtype=object)
    states[com_dists < 40] = "bound"
    states[(com_dists >= 40) & (com_dists < 48)] = "transition"
    states[com_dists >= 48] = "unbound"
    
    state_counts = {
        "bound": np.sum(states == "bound"),
        "transition": np.sum(states == "transition"),
        "unbound": np.sum(states == "unbound"),
    }
    
    # Time to unbind (first frame entering unbound)
    unbound_frames = np.where(states == "unbound")[0]
    time_to_unbind = time_ns[unbound_frames[0]] if len(unbound_frames) > 0 else None
    
    print(f"  RMSD: {rmsd_vals.mean():.2f} ± {rmsd_vals.std():.2f} Å")
    print(f"  COM: {com_dists.mean():.2f} ± {com_dists.std():.2f} Å, range [{com_dists.min():.2f}, {com_dists.max():.2f}]")
    print(f"  Rg: {rg_vals.mean():.2f} ± {rg_vals.std():.2f} Å")
    print(f"  Contacts: {contact_counts.mean():.0f} ± {contact_counts.std():.0f}")
    print(f"  H-bonds: {hbond_counts.mean():.1f} ± {hbond_counts.std():.1f}")
    print(f"  States: bound={state_counts['bound']}, transition={state_counts['transition']}, unbound={state_counts['unbound']}")
    if time_to_unbind:
        print(f"  Time to unbind: {time_to_unbind:.1f} ns")
    else:
        print(f"  Time to unbind: N/A (never unbound)")
    
    return {
        "name": name,
        "time_ns": time_ns,
        "rmsd": rmsd_vals,
        "rmsf": rmsf_vals,
        "ca_resids": ca_resids,
        "com_dists": com_dists,
        "rg": rg_vals,
        "contacts": contact_counts,
        "hbonds": hbond_counts,
        "mut_distances": mut_distances,
        "states": state_counts,
        "time_to_unbind": time_to_unbind,
        "n_frames": n_frames,
    }

# ---- Run analysis for all 6 replicas ----
all_results = []
for name, (top, traj) in systems.items():
    try:
        result = analyze_one(name, top, traj)
        all_results.append(result)
    except Exception as e:
        print(f"ERROR analyzing {name}: {e}")
        import traceback
        traceback.print_exc()

# ---- Save raw data ----
print("\nSaving raw data...")
for r in all_results:
    np.savez(
        outdir / f"{r['name']}_data.npz",
        time_ns=r["time_ns"],
        rmsd=r["rmsd"],
        rmsf=r["rmsf"],
        ca_resids=r["ca_resids"],
        com_dists=r["com_dists"],
        rg=r["rg"],
        contacts=r["contacts"],
        hbonds=r["hbonds"],
    )
    for mut_name, dists in r["mut_distances"].items():
        np.save(outdir / f"{r['name']}_mut_{mut_name}_dist.npy", dists)

# ---- Summary JSON ----
summary = {}
for r in all_results:
    summary[r["name"]] = {
        "rmsd_mean": float(r["rmsd"].mean()),
        "rmsd_std": float(r["rmsd"].std()),
        "rmsd_max": float(r["rmsd"].max()),
        "com_mean": float(r["com_dists"].mean()),
        "com_std": float(r["com_dists"].std()),
        "com_min": float(r["com_dists"].min()),
        "com_max": float(r["com_dists"].max()),
        "rg_mean": float(r["rg"].mean()),
        "rg_std": float(r["rg"].std()),
        "contacts_mean": float(r["contacts"].mean()),
        "contacts_std": float(r["contacts"].std()),
        "hbonds_mean": float(r["hbonds"].mean()),
        "hbonds_std": float(r["hbonds"].std()),
        "states": r["states"],
        "time_to_unbind_ns": r["time_to_unbind"],
    }

with open(outdir / "summary.json", "w") as f:
    json.dump(summary, f, indent=2)

# ---- Cross-replica comparison plots ----
print("\nGenerating comparison plots...")

fig, axes = plt.subplots(3, 2, figsize=(12, 14))

def plot_metric(ax, metric_key, ylabel, title):
    for r in all_results:
        color = "tab:blue" if "WT" in r["name"] else "tab:orange"
        alpha = 0.7
        lw = 0.6
        ax.plot(r["time_ns"], r[metric_key], color=color, alpha=alpha, lw=lw, label=r["name"])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

plot_metric(axes[0, 0], "rmsd", "CA RMSD (Å)", "Backbone RMSD")
plot_metric(axes[0, 1], "com_dists", "COM distance (Å)", "cGAS-TRIM41 COM Distance")
plot_metric(axes[1, 0], "rg", "Rg (Å)", "Radius of Gyration")
plot_metric(axes[1, 1], "contacts", "Contact count", "Interface Contacts (<5Å)")
plot_metric(axes[2, 0], "hbonds", "H-bond count", "Interface H-bonds (<3.5Å)")

# Mutation site distances
ax = axes[2, 1]
for r in all_results:
    color = "tab:blue" if "WT" in r["name"] else "tab:orange"
    for mut_name, dists in r["mut_distances"].items():
        ax.plot(r["time_ns"], dists, color=color, alpha=0.4, lw=0.5)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Distance to nearest TRIM41 atom (Å)")
ax.set_title("Mutation Site Distances (all 4 sites)")

# Add legend for WT vs 4mut
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color="tab:blue", lw=2, label="Hsap_WT (3 reps)"),
    Line2D([0], [0], color="tab:orange", lw=2, label="Hsap_4mut (3 reps)"),
]
fig.legend(handles=legend_elements, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.98))

fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(outdir / "all_replicas_comparison.png", dpi=300)
plt.close(fig)

# ---- State occupancy bar plot ----
fig, ax = plt.subplots(figsize=(10, 4))
names = [r["name"] for r in all_results]
bound = [r["states"]["bound"] / r["n_frames"] * 100 for r in all_results]
transition = [r["states"]["transition"] / r["n_frames"] * 100 for r in all_results]
unbound = [r["states"]["unbound"] / r["n_frames"] * 100 for r in all_results]

x = np.arange(len(names))
width = 0.25
ax.bar(x - width, bound, width, label="Bound (<40Å)", color="tab:green")
ax.bar(x, transition, width, label="Transition (40-48Å)", color="tab:orange")
ax.bar(x + width, unbound, width, label="Unbound (>48Å)", color="tab:red")
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=45, ha="right")
ax.set_ylabel("Occupancy (%)")
ax.set_title("State Occupancy per Replica")
ax.legend()
fig.tight_layout()
fig.savefig(outdir / "state_occupancy.png", dpi=300)
plt.close(fig)

# ---- Time to unbind comparison ----
fig, ax = plt.subplots(figsize=(8, 4))
ttb = [r["time_to_unbind"] if r["time_to_unbind"] is not None else 200 for r in all_results]
colors = ["tab:blue" if "WT" in n else "tab:orange" for n in names]
ax.bar(names, ttb, color=colors, alpha=0.7)
ax.set_ylabel("Time to unbind (ns)")
ax.set_title("First Entry into Unbound State (COM > 48Å)")
ax.axhline(200, color="gray", ls="--", alpha=0.5, label="Never unbound")
ax.set_xticklabels(names, rotation=45, ha="right")
fig.tight_layout()
fig.savefig(outdir / "time_to_unbind.png", dpi=300)
plt.close(fig)

# ---- RMSF comparison (WT avg vs 4mut avg) ----
fig, ax = plt.subplots(figsize=(10, 4))
wt_rmsf = [r["rmsf"] for r in all_results if "WT" in r["name"]]
mut_rmsf = [r["rmsf"] for r in all_results if "4mut" in r["name"]]
wt_resids = all_results[0]["ca_resids"]

if len(wt_rmsf) == len(mut_rmsf) == 3:
    wt_mean = np.mean(wt_rmsf, axis=0)
    wt_std = np.std(wt_rmsf, axis=0)
    mut_mean = np.mean(mut_rmsf, axis=0)
    mut_std = np.std(mut_rmsf, axis=0)
    
    ax.plot(wt_resids, wt_mean, color="tab:blue", label="WT mean ± std", lw=1)
    ax.fill_between(wt_resids, wt_mean - wt_std, wt_mean + wt_std, color="tab:blue", alpha=0.2)
    ax.plot(wt_resids, mut_mean, color="tab:orange", label="4mut mean ± std", lw=1)
    ax.fill_between(wt_resids, mut_mean - mut_std, mut_mean + mut_std, color="tab:orange", alpha=0.2)
    ax.axvline(218.5, color="gray", ls="--", alpha=0.5)
    ax.text(220, ax.get_ylim()[1]*0.9, "cGAS | TRIM41", fontsize=9)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("CA RMSF (Å)")
    ax.set_title("Average RMSF: WT vs 4mut (3 replicas each)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "rmsf_avg_comparison.png", dpi=300)
    plt.close(fig)

# ---- Print summary table ----
print("\n" + "="*80)
print("SUMMARY TABLE: All 6 Replicas")
print("="*80)
print(f"{'Replica':<18} {'RMSD':>10} {'COM':>10} {'Rg':>8} {'Contacts':>10} {'HBonds':>8} {'Unbind':>8}")
print("-"*80)
for r in all_results:
    s = summary[r["name"]]
    unbind_str = f"{s['time_to_unbind_ns']:.1f}ns" if s['time_to_unbind_ns'] else "N/A"
    print(f"{r['name']:<18} {s['rmsd_mean']:>9.2f}Å {s['com_mean']:>9.2f}Å {s['rg_mean']:>7.2f}Å {s['contacts_mean']:>9.0f} {s['hbonds_mean']:>7.1f} {unbind_str:>8}")

print("\nPlots saved to:")
for f in ["all_replicas_comparison.png", "state_occupancy.png", "time_to_unbind.png", "rmsf_avg_comparison.png"]:
    print(f"  {outdir / f}")
print(f"\nData saved to: {outdir}")
print("✅ Batch analysis complete.")
