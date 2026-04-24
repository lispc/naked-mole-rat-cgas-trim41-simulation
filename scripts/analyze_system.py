#!/usr/bin/env python3
"""
Generic single-system MD analysis pipeline.
Analyzes one or more replicas for a given system, producing:
  - RMSD time series per replica
  - RMSF per residue (cross-replica average)
  - Interface contact occupancy (heatmap + top list)
  - COM distance time series
  - Active site distances (if configured)
Outputs: PNG plots + JSON summary to a specified directory.
"""
import argparse
import json
import sys
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.distances import distance_array


def analyze_replica(prmtop, dcd, trim_range, cgas_range, active_sites=None,
                    dt_ns=0.1, contact_cutoff=5.0):
    """
    Analyze a single replica trajectory.
    
    Returns dict with analysis results.
    """
    print(f"  Loading: {dcd}")
    u = mda.Universe(prmtop, dcd)
    n_frames = len(u.trajectory)
    traj_ns = n_frames * dt_ns
    print(f"  Frames: {n_frames} ({traj_ns:.2f} ns)")
    
    t_start, t_end = trim_range
    c_start, c_end = cgas_range
    
    trim_ca = u.select_atoms(f"resid {t_start}-{t_end} and name CA")
    cgas_ca = u.select_atoms(f"resid {c_start}-{c_end} and name CA")
    protein_ca = u.select_atoms("protein and name CA")
    
    print(f"  TRIM41 CA: {len(trim_ca)}, cGAS CA: {len(cgas_ca)}")
    
    # 1. RMSD (aligned on protein CA)
    print("  [1/5] RMSD...")
    rmsd_data = []
    ref = mda.Universe(prmtop, dcd)
    aligner = align.AlignTraj(u, ref, select="protein and name CA", in_memory=True)
    aligner.run()
    
    for ts in u.trajectory:
        r = rms.rmsd(protein_ca.positions, ref.select_atoms("protein and name CA").positions, superposition=False)
        rmsd_data.append(r)
    rmsd_arr = np.array(rmsd_data)
    n_frames = len(rmsd_arr)  # Use actual count, may differ from len(u.trajectory)
    print(f"    Mean RMSD: {rmsd_arr.mean():.2f} ± {rmsd_arr.std():.2f} Å")
    
    # 2. RMSF (need to reload and align since RMSD modified trajectory)
    print("  [2/5] RMSF...")
    protein_ca2 = protein_ca  # reuse aligned universe
    rmsf = rms.RMSF(protein_ca2).run()
    rmsf_vals = rmsf.results.rmsf
    rmsf_resids = protein_ca2.resids
    
    trim_mask = (rmsf_resids >= t_start) & (rmsf_resids <= t_end)
    cgas_mask = (rmsf_resids >= c_start) & (rmsf_resids <= c_end)
    print(f"    TRIM41 RMSF: {rmsf_vals[trim_mask].mean():.2f} ± {rmsf_vals[trim_mask].std():.2f} Å")
    print(f"    cGAS RMSF:   {rmsf_vals[cgas_mask].mean():.2f} ± {rmsf_vals[cgas_mask].std():.2f} Å")
    
    # 3. Interface contacts (CA-CA < cutoff)
    print(f"  [3/5] Interface contacts (cutoff {contact_cutoff} Å)...")
    contact_counts = Counter()
    
    for ts in u.trajectory:
        dist_mat = distance_array(trim_ca.positions, cgas_ca.positions)
        contacts = np.argwhere(dist_mat < contact_cutoff)
        for i, j in contacts:
            pair = (int(trim_ca.resids[i]), int(cgas_ca.resids[j]))
            contact_counts[pair] += 1
    
    occupancy = {f"{t}_{c}": count / n_frames for (t, c), count in contact_counts.items()}
    top_contacts = sorted(occupancy.items(), key=lambda x: x[1], reverse=True)[:20]
    print(f"    Unique contacts: {len(occupancy)}")
    for key, occ in top_contacts[:5]:
        t, c = key.split('_')
        print(f"      TRIM41-{t} -- cGAS-{c}: {occ:.3f}")
    
    # 4. COM distance
    print("  [4/5] COM distance...")
    com_distances = []
    for ts in u.trajectory:
        trim_com = trim_ca.center_of_mass()
        cgas_com = cgas_ca.center_of_mass()
        com_dist = np.linalg.norm(trim_com - cgas_com)
        com_distances.append(com_dist)
    com_arr = np.array(com_distances)
    print(f"    COM: {com_arr.mean():.2f} ± {com_arr.std():.2f} Å")
    
    # 5. Active site distances (optional)
    active_results = {}
    if active_sites:
        print("  [5/5] Active site distances...")
        sel_parts = []
        for name, resid in active_sites.items():
            sel_parts.append(f"(resid {resid} and name CA)")
        active_sel = " or ".join(sel_parts)
        active_atoms = u.select_atoms(active_sel)
        
        for idx, (name, resid) in enumerate(active_sites.items()):
            dists = []
            for ts in u.trajectory:
                dist_mat = distance_array(active_atoms[idx].position[None, :], trim_ca.positions)
                min_dist = dist_mat[0].min()
                dists.append(min_dist)
            arr = np.array(dists)
            active_results[name] = {
                "mean": float(arr.mean()),
                "std": float(arr.std()),
                "min": float(arr.min()),
                "max": float(arr.max()),
                "values": arr.tolist(),
            }
            print(f"    {name} -> TRIM41: {arr.mean():.2f} ± {arr.std():.2f} Å")
    else:
        print("  [5/5] Active site distances: skipped (not configured)")
    
    return {
        "n_frames": n_frames,
        "traj_ns": n_frames * dt_ns,
        "rmsd": rmsd_arr.tolist(),
        "rmsd_mean": float(rmsd_arr.mean()),
        "rmsd_std": float(rmsd_arr.std()),
        "rmsf_resids": rmsf_resids.tolist(),
        "rmsf_vals": rmsf_vals.tolist(),
        "rmsf_trim_mean": float(rmsf_vals[trim_mask].mean()),
        "rmsf_cgas_mean": float(rmsf_vals[cgas_mask].mean()),
        "occupancy": occupancy,
        "n_unique_contacts": len(occupancy),
        "top_contacts": top_contacts,
        "com_distances": com_arr.tolist(),
        "com_mean": float(com_arr.mean()),
        "com_std": float(com_arr.std()),
        "active_sites": active_results,
    }


def plot_single_system(results, outdir, system_name, dt_ns=0.1):
    """Generate all plots for a single system."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(results)))
    
    # 1. RMSD time series (all replicas)
    fig, ax = plt.subplots(figsize=(10, 3.5))
    for i, (rep_name, data) in enumerate(results.items()):
        time = np.arange(data["n_frames"]) * dt_ns
        rmsd = np.array(data["rmsd"])
        ax.plot(time, rmsd, lw=0.5, alpha=0.7, color=colors[i], label=f"{rep_name} ({rmsd.mean():.2f}Å)")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title(f"{system_name} — Backbone RMSD")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"{system_name}_rmsd.png", dpi=200)
    plt.close(fig)
    print(f"  Saved: {outdir / f'{system_name}_rmsd.png'}")
    
    # 2. RMSF (cross-replica overlay)
    fig, ax = plt.subplots(figsize=(10, 4))
    for i, (rep_name, data) in enumerate(results.items()):
        ax.plot(data["rmsf_resids"], data["rmsf_vals"], lw=0.6, alpha=0.7, color=colors[i], label=rep_name)
    # Draw chain boundary
    cgas_start = min(data["rmsf_resids"] for data in results.values())
    # Find where cGAS starts (first resid > 218)
    sample_resids = np.array(list(results.values())[0]["rmsf_resids"])
    boundary = sample_resids[sample_resids > 218][0] if any(sample_resids > 218) else 219
    ax.axvline(boundary - 0.5, color="black", ls="--", alpha=0.3)
    ax.text(boundary - 0.5, ax.get_ylim()[1] * 0.95, "cGAS →", ha="right", fontsize=8, alpha=0.5)
    ax.text(boundary + 0.5, ax.get_ylim()[1] * 0.95, "← TRIM41", ha="left", fontsize=8, alpha=0.5)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title(f"{system_name} — CA-RMSF per residue")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"{system_name}_rmsf.png", dpi=200)
    plt.close(fig)
    print(f"  Saved: {outdir / f'{system_name}_rmsf.png'}")
    
    # 3. COM distance
    fig, ax = plt.subplots(figsize=(10, 3.5))
    for i, (rep_name, data) in enumerate(results.items()):
        time = np.arange(data["n_frames"]) * dt_ns
        ax.plot(time, data["com_distances"], lw=0.5, alpha=0.7, color=colors[i], label=f"{rep_name}")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("COM distance (Å)")
    ax.set_title(f"{system_name} — TRIM41-cGAS COM distance")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"{system_name}_com.png", dpi=200)
    plt.close(fig)
    print(f"  Saved: {outdir / f'{system_name}_com.png'}")
    
    # 4. Contact heatmap (first replica only, for clarity)
    first_rep = list(results.keys())[0]
    occ = results[first_rep]["occupancy"]
    if occ:
        trim_res = sorted(set([int(k.split('_')[0]) for k in occ.keys()]))
        cgas_res = sorted(set([int(k.split('_')[1]) for k in occ.keys()]))
        
        mat = np.zeros((len(trim_res), len(cgas_res)))
        for i, tr in enumerate(trim_res):
            for j, cr in enumerate(cgas_res):
                key = f"{tr}_{cr}"
                mat[i, j] = occ.get(key, 0)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", origin="lower")
        ax.set_xlabel("cGAS resid")
        ax.set_ylabel("TRIM41 resid")
        ax.set_title(f"{system_name} — Interface contact occupancy ({first_rep})")
        plt.colorbar(im, ax=ax, label="Occupancy")
        fig.tight_layout()
        fig.savefig(outdir / f"{system_name}_contacts.png", dpi=200)
        plt.close(fig)
        print(f"  Saved: {outdir / f'{system_name}_contacts.png'}")
    
    # 5. Active site distances (if available)
    first_data = list(results.values())[0]
    if first_data.get("active_sites"):
        n_sites = len(first_data["active_sites"])
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.flatten()
        for idx, name in enumerate(first_data["active_sites"].keys()):
            ax = axes[idx]
            for i, (rep_name, data) in enumerate(results.items()):
                dists = np.array(data["active_sites"][name]["values"])
                ax.hist(dists, bins=30, alpha=0.4, color=colors[i],
                        label=f"{rep_name} (μ={dists.mean():.1f}Å)", density=True)
            ax.set_xlabel("Distance to nearest TRIM41 CA (Å)")
            ax.set_ylabel("Density")
            ax.set_title(f"{name}")
            ax.legend(fontsize=7)
        fig.suptitle(f"{system_name} — Active site distances", fontsize=12)
        fig.tight_layout()
        fig.savefig(outdir / f"{system_name}_active_sites.png", dpi=200)
        plt.close(fig)
        print(f"  Saved: {outdir / f'{system_name}_active_sites.png'}")


def main():
    parser = argparse.ArgumentParser(description="Analyze MD system")
    parser.add_argument("--system", required=True, help="System name (e.g. Hgal_WT)")
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--trajectories", nargs="+", required=True,
                        help="One or more DCD files (replicas)")
    parser.add_argument("--replica-names", nargs="+", default=None,
                        help="Names for each trajectory (default: rep1, rep2, ...)")
    parser.add_argument("--trim-range", nargs=2, type=int, default=[1, 218],
                        help="TRIM41 resid range (start end)")
    parser.add_argument("--cgas-range", nargs=2, type=int, required=True,
                        help="cGAS resid range (start end)")
    parser.add_argument("--active-sites", default=None,
                        help="JSON dict of active sites, e.g. '{\"S463\": 482, \"E511\": 530}'")
    parser.add_argument("--dt-ns", type=float, default=0.1,
                        help="Time per frame in ns (default 0.1 = 100ps)")
    parser.add_argument("--contact-cutoff", type=float, default=5.0)
    parser.add_argument("--outdir", default="data/analysis")
    args = parser.parse_args()
    
    active_sites = None
    if args.active_sites:
        active_sites = json.loads(args.active_sites)
    
    replica_names = args.replica_names or [f"rep{i+1}" for i in range(len(args.trajectories))]
    if len(replica_names) != len(args.trajectories):
        print("Error: --replica-names must match --trajectories count")
        sys.exit(1)
    
    print(f"[{datetime.now()}] Analyzing system: {args.system}")
    print(f"  TRIM41: resid {args.trim_range[0]}-{args.trim_range[1]}")
    print(f"  cGAS:   resid {args.cgas_range[0]}-{args.cgas_range[1]}")
    print(f"  Replicas: {len(args.trajectories)}")
    
    all_results = {}
    for rep_name, dcd in zip(replica_names, args.trajectories):
        print(f"\n{'='*60}")
        print(f"Replica: {rep_name}")
        print(f"{'='*60}")
        result = analyze_replica(
            args.prmtop, dcd,
            tuple(args.trim_range), tuple(args.cgas_range),
            active_sites=active_sites,
            dt_ns=args.dt_ns,
            contact_cutoff=args.contact_cutoff,
        )
        all_results[rep_name] = result
    
    print(f"\n{'='*60}")
    print("Generating plots...")
    print(f"{'='*60}")
    plot_single_system(all_results, args.outdir, args.system, dt_ns=args.dt_ns)
    
    # Save JSON summary
    summary_path = Path(args.outdir) / f"{args.system}_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved summary: {summary_path}")
    
    print(f"[{datetime.now()}] Done.")


if __name__ == "__main__":
    main()
