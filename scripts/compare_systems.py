#!/usr/bin/env python3
"""
Cross-system comparison for MD analysis results.
Compares WT vs mutant (or any two systems) on:
  - RMSD distributions
  - RMSF differences (ΔRMSF = mutant - WT)
  - Interface contact differences
  - COM distance distributions
  - Active site distances (if available)

Outputs: comparison plots + statistical test results.
"""
import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def load_summary(path):
    with open(path) as f:
        return json.load(f)


def welch_ttest(a, b):
    """Welch's t-test (unequal variances)."""
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return t, p


def compare_rmsd(data_a, data_b, name_a, name_b, outdir):
    """RMSD distribution comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    all_rmsd_a = []
    all_rmsd_b = []
    colors_a = plt.cm.Blues(np.linspace(0.4, 0.9, len(data_a)))
    colors_b = plt.cm.Reds(np.linspace(0.4, 0.9, len(data_b)))
    
    ax = axes[0]
    for i, (rep, d) in enumerate(data_a.items()):
        r = np.array(d["rmsd"])
        all_rmsd_a.extend(r)
        ax.hist(r, bins=50, alpha=0.4, color=colors_a[i], label=f"{name_a} {rep}", density=True)
    for i, (rep, d) in enumerate(data_b.items()):
        r = np.array(d["rmsd"])
        all_rmsd_b.extend(r)
        ax.hist(r, bins=50, alpha=0.4, color=colors_b[i], label=f"{name_b} {rep}", density=True)
    ax.set_xlabel("RMSD (Å)")
    ax.set_ylabel("Density")
    ax.set_title("RMSD distributions")
    ax.legend(fontsize=7)
    
    ax = axes[1]
    ax.boxplot([all_rmsd_a, all_rmsd_b], tick_labels=[name_a, name_b])
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("RMSD comparison")
    
    t, p = welch_ttest(all_rmsd_a, all_rmsd_b)
    fig.suptitle(f"RMSD: {name_a} vs {name_b} (Welch t={t:.2f}, p={p:.2e})", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / f"compare_{name_a}_vs_{name_b}_rmsd.png", dpi=200)
    plt.close(fig)
    print(f"  Saved RMSD comparison (p={p:.2e})")
    return {"rmsd_t": float(t), "rmsd_p": float(p)}


def compare_rmsf(data_a, data_b, name_a, name_b, outdir):
    """RMSF difference map."""
    # Use first replica from each (or average if multiple)
    first_a = list(data_a.values())[0]
    first_b = list(data_b.values())[0]
    
    resids_a = np.array(first_a["rmsf_resids"])
    resids_b = np.array(first_b["rmsf_resids"])
    
    # Interpolate to common grid if needed (should be same for same species)
    if len(resids_a) == len(resids_b) and np.allclose(resids_a, resids_b):
        resids = resids_a
        drmsf = np.array(first_b["rmsf_vals"]) - np.array(first_a["rmsf_vals"])
    else:
        # Interpolate B onto A's grid
        from scipy.interpolate import interp1d
        f = interp1d(resids_b, first_b["rmsf_vals"], kind="linear", bounds_error=False, fill_value=0)
        resids = resids_a
        drmsf = f(resids) - np.array(first_a["rmsf_vals"])
    
    fig, ax = plt.subplots(figsize=(10, 4))
    colors = ["blue" if x < 0 else "red" for x in drmsf]
    ax.bar(resids, drmsf, color=colors, width=1.0, alpha=0.6)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("ΔRMSF (Å)")
    ax.set_title(f"RMSF difference: {name_b} − {name_a}")
    fig.tight_layout()
    fig.savefig(outdir / f"compare_{name_a}_vs_{name_b}_rmsf.png", dpi=200)
    plt.close(fig)
    print(f"  Saved ΔRMSF plot")
    
    # Report significant changes (> 1Å)
    sig_idx = np.where(np.abs(drmsf) > 1.0)[0]
    print(f"    Residues with |ΔRMSF| > 1.0 Å: {len(sig_idx)}")
    return {"drmsf": drmsf.tolist(), "resids": resids.tolist()}


def compare_contacts(data_a, data_b, name_a, name_b, outdir):
    """Compare interface contacts between systems."""
    # Aggregate occupancy across replicas
    def aggregate_occupancy(data):
        all_occ = {}
        total_frames = 0
        for rep_data in data.values():
            n = rep_data["n_frames"]
            total_frames += n
            for key, occ in rep_data["occupancy"].items():
                all_occ[key] = all_occ.get(key, 0) + occ * n
        for key in all_occ:
            all_occ[key] /= total_frames
        return all_occ
    
    occ_a = aggregate_occupancy(data_a)
    occ_b = aggregate_occupancy(data_b)
    
    all_keys = sorted(set(occ_a.keys()) | set(occ_b.keys()))
    
    vals_a = [occ_a.get(k, 0) for k in all_keys]
    vals_b = [occ_b.get(k, 0) for k in all_keys]
    diff = np.array(vals_b) - np.array(vals_a)
    
    # Top lost / gained contacts
    lost = sorted([(all_keys[i], diff[i]) for i in range(len(all_keys)) if diff[i] < -0.1],
                  key=lambda x: x[1])[:10]
    gained = sorted([(all_keys[i], diff[i]) for i in range(len(all_keys)) if diff[i] > 0.1],
                    key=lambda x: x[1], reverse=True)[:10]
    
    print(f"    Contacts lost (Δ < -0.1): {len(lost)}")
    for key, d in lost[:5]:
        t, c = key.split('_')
        print(f"      TRIM41-{t} -- cGAS-{c}: Δ={d:.3f}")
    print(f"    Contacts gained (Δ > +0.1): {len(gained)}")
    for key, d in gained[:5]:
        t, c = key.split('_')
        print(f"      TRIM41-{t} -- cGAS-{c}: Δ={d:.3f}")
    
    # Scatter plot
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(vals_a, vals_b, alpha=0.3, s=10)
    ax.plot([0, 1], [0, 1], "k--", lw=0.5)
    ax.set_xlabel(f"{name_a} occupancy")
    ax.set_ylabel(f"{name_b} occupancy")
    ax.set_title("Interface contact occupancy comparison")
    fig.tight_layout()
    fig.savefig(outdir / f"compare_{name_a}_vs_{name_b}_contacts.png", dpi=200)
    plt.close(fig)
    print(f"  Saved contact comparison")
    
    return {"lost": lost, "gained": gained}


def compare_com(data_a, data_b, name_a, name_b, outdir):
    """Compare COM distance distributions."""
    all_com_a = []
    all_com_b = []
    for d in data_a.values():
        all_com_a.extend(d["com_distances"])
    for d in data_b.values():
        all_com_b.extend(d["com_distances"])
    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(all_com_a, bins=50, alpha=0.5, label=f"{name_a} (μ={np.mean(all_com_a):.1f}Å)", density=True)
    ax.hist(all_com_b, bins=50, alpha=0.5, label=f"{name_b} (μ={np.mean(all_com_b):.1f}Å)", density=True)
    ax.set_xlabel("COM distance (Å)")
    ax.set_ylabel("Density")
    ax.set_title("TRIM41-cGAS COM distance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / f"compare_{name_a}_vs_{name_b}_com.png", dpi=200)
    plt.close(fig)
    
    t, p = welch_ttest(all_com_a, all_com_b)
    print(f"  Saved COM comparison (p={p:.2e})")
    return {"com_t": float(t), "com_p": float(p)}


def main():
    parser = argparse.ArgumentParser(description="Compare two MD systems")
    parser.add_argument("--a", required=True, help="System A JSON summary")
    parser.add_argument("--b", required=True, help="System B JSON summary")
    parser.add_argument("--name-a", default="WT")
    parser.add_argument("--name-b", default="Mutant")
    parser.add_argument("--outdir", default="data/analysis")
    args = parser.parse_args()
    
    print(f"[{datetime.now()}] Comparing {args.name_a} vs {args.name_b}")
    data_a = load_summary(args.a)
    data_b = load_summary(args.b)
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    print("\n[1/4] RMSD comparison...")
    results.update(compare_rmsd(data_a, data_b, args.name_a, args.name_b, outdir))
    
    print("\n[2/4] RMSF comparison...")
    results.update(compare_rmsf(data_a, data_b, args.name_a, args.name_b, outdir))
    
    print("\n[3/4] Contact comparison...")
    results.update(compare_contacts(data_a, data_b, args.name_a, args.name_b, outdir))
    
    print("\n[4/4] COM comparison...")
    results.update(compare_com(data_a, data_b, args.name_a, args.name_b, outdir))
    
    summary_path = outdir / f"compare_{args.name_a}_vs_{args.name_b}_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved comparison summary: {summary_path}")


if __name__ == "__main__":
    main()
