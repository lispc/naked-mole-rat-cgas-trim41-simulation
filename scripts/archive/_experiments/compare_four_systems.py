#!/usr/bin/env python3
"""
Four-system comprehensive MD comparison.
Systems: Hgal_WT, Hgal_4mut_rev, Hsap_WT, Hsap_4mut

Generates:
  - 4-system RMSD boxplot + stats
  - 4-system COM distance distributions
  - Intra-species ΔRMSF (WT vs mutant)
  - Cross-species active site distance comparison (homolog-mapped)
  - Interface contact summary table
  - JSON summary with all statistics

FIXED (2026-04-27): Replaced naive Welch t-test on correlated MD frames
with effective-sample-size t-test using autocorrelation time estimation.
"""
import argparse
import json
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def load_summary(path):
    with open(path) as f:
        return json.load(f)


def effective_sample_size(data, max_lag=None):
    """Estimate effective sample size using integrated autocorrelation time."""
    data = np.asarray(data, dtype=float)
    n = len(data)
    if n < 4:
        return float(n), 1.0
    
    if max_lag is None:
        max_lag = min(n // 3, 2000)
    
    data = data - np.mean(data)
    c0 = np.mean(data ** 2)
    if c0 == 0:
        return float(n), 1.0
    
    tau = 1.0
    for lag in range(0, max_lag - 1, 2):
        c_lag = np.mean(data[:n-lag] * data[lag:]) / c0
        c_lag1 = np.mean(data[:n-(lag+1)] * data[lag+1:]) / c0 if lag + 1 < n else 0.0
        gamma = c_lag + c_lag1
        if gamma < 0:
            break
        tau += 2 * gamma
    
    tau = max(tau, 1.0)
    n_eff = n / tau
    return float(n_eff), float(tau)


def correlated_ttest(a, b):
    """Welch-style t-test for correlated time series."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    
    n_eff_a, tau_a = effective_sample_size(a)
    n_eff_b, tau_b = effective_sample_size(b)
    
    mean_a, mean_b = np.mean(a), np.mean(b)
    var_a = np.var(a, ddof=1)
    var_b = np.var(b, ddof=1)
    
    se = np.sqrt(var_a / n_eff_a + var_b / n_eff_b)
    if se == 0:
        return 0.0, 1.0, n_eff_a, n_eff_b, tau_a, tau_b
    
    t = (mean_a - mean_b) / se
    num = (var_a / n_eff_a + var_b / n_eff_b) ** 2
    den = (var_a / n_eff_a) ** 2 / max(n_eff_a - 1, 1) + (var_b / n_eff_b) ** 2 / max(n_eff_b - 1, 1)
    df = num / den if den > 0 else 1.0
    p = 2 * (1 - stats.t.cdf(abs(t), df))
    return float(t), float(p), n_eff_a, n_eff_b, tau_a, tau_b


def flatten_replicas(data, key):
    """Flatten a metric across all replicas."""
    vals = []
    for rep in data.values():
        vals.extend(rep[key])
    return np.array(vals)


def plot_four_system_rmsd(summaries, names, outdir):
    """Boxplot of RMSD for all four systems."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    all_data = []
    labels = []
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12"]
    
    for name, summary in zip(names, summaries):
        rmsd = flatten_replicas(summary, "rmsd")
        all_data.append(rmsd)
        labels.append(f"{name}\nμ={rmsd.mean():.2f}Å\nσ={rmsd.std():.2f}Å")
    
    bp = ax.boxplot(all_data, labels=labels, patch_artist=True,
                     medianprops=dict(color="black", lw=1.5))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Backbone RMSD — Four Systems (200ns)")
    fig.tight_layout()
    fig.savefig(outdir / "four_systems_rmsd.png", dpi=200)
    plt.close(fig)
    print("  Saved: four_systems_rmsd.png")
    
    # Pairwise stats
    stats_results = {}
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            t, p, n_eff_i, n_eff_j, tau_i, tau_j = correlated_ttest(all_data[i], all_data[j])
            key = f"{names[i]}_vs_{names[j]}"
            stats_results[key] = {
                "t": float(t), "p": float(p),
                "n_eff_i": float(n_eff_i), "n_eff_j": float(n_eff_j),
                "tau_i": float(tau_i), "tau_j": float(tau_j),
                "mean_i": float(all_data[i].mean()),
                "mean_j": float(all_data[j].mean()),
            }
            print(f"    {key}: t={t:.2f}, p={p:.3e}, N_eff={n_eff_i:.0f}/{n_eff_j:.0f}")
    return stats_results


def plot_four_system_com(summaries, names, outdir):
    """COM distance distributions for all four systems."""
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12"]
    
    for name, summary, color in zip(names, summaries, colors):
        com = flatten_replicas(summary, "com_distances")
        ax.hist(com, bins=50, alpha=0.4, color=color,
                label=f"{name} (μ={com.mean():.1f}Å)", density=True)
    
    ax.set_xlabel("COM distance (Å)")
    ax.set_ylabel("Density")
    ax.set_title("TRIM41-cGAS COM Distance — Four Systems")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / "four_systems_com.png", dpi=200)
    plt.close(fig)
    print("  Saved: four_systems_com.png")
    
    stats_results = {}
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            com_i = flatten_replicas(summaries[i], "com_distances")
            com_j = flatten_replicas(summaries[j], "com_distances")
            t, p, n_eff_i, n_eff_j, tau_i, tau_j = correlated_ttest(com_i, com_j)
            key = f"{names[i]}_vs_{names[j]}"
            stats_results[key] = {
                "t": float(t), "p": float(p),
                "n_eff_i": float(n_eff_i), "n_eff_j": float(n_eff_j),
                "tau_i": float(tau_i), "tau_j": float(tau_j),
                "mean_i": float(com_i.mean()),
                "mean_j": float(com_j.mean()),
            }
            print(f"    COM {key}: t={t:.2f}, p={p:.3e}, N_eff={n_eff_i:.0f}/{n_eff_j:.0f}")
    return stats_results


def plot_intrinsic_rmsf_delta(summary_wt, summary_mut, name_wt, name_mut, outdir):
    """ΔRMSF within same species."""
    wt = list(summary_wt.values())[0]
    mut = list(summary_mut.values())[0]
    
    resids = np.array(wt["rmsf_resids"])
    rmsf_wt = np.array(wt["rmsf_vals"])
    rmsf_mut = np.array(mut["rmsf_vals"])
    
    if len(resids) != len(rmsf_mut):
        from scipy.interpolate import interp1d
        mut_resids = np.array(mut["rmsf_resids"])
        f = interp1d(mut_resids, rmsf_mut, kind="linear", bounds_error=False, fill_value=0)
        rmsf_mut = f(resids)
    
    drmsf = rmsf_mut - rmsf_wt
    
    fig, ax = plt.subplots(figsize=(12, 4))
    colors = ["blue" if x < 0 else "red" for x in drmsf]
    ax.bar(resids, drmsf, color=colors, width=1.0, alpha=0.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(218.5, color="black", ls="--", alpha=0.3)
    ax.text(100, ax.get_ylim()[1]*0.9, "TRIM41", ha="center", fontsize=9, alpha=0.5)
    ax.text(400, ax.get_ylim()[1]*0.9, "cGAS", ha="center", fontsize=9, alpha=0.5)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("ΔRMSF (Å)")
    ax.set_title(f"RMSF: {name_mut} − {name_wt}")
    fig.tight_layout()
    fig.savefig(outdir / f"delta_rmsf_{name_wt}_vs_{name_mut}.png", dpi=200)
    plt.close(fig)
    print(f"  Saved: delta_rmsf_{name_wt}_vs_{name_mut}.png")
    
    sig = np.where(np.abs(drmsf) > 1.0)[0]
    print(f"    |ΔRMSF| > 1.0 Å: {len(sig)} residues")
    return {"drmsf": drmsf.tolist(), "resids": resids.tolist(), "n_significant": int(len(sig))}


def compare_active_sites(summaries, names, outdir):
    """
    Compare active site distances across all four systems.
    Homolog mapping: Hgal S463↔Hsap D431, E511↔K479, Y527↔L495, T530↔K498
    """
    hgal_sites = ["S463", "E511", "Y527", "T530"]
    hsap_sites = ["D431", "K479", "L495", "K498"]
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12"]
    
    results = {}
    
    for idx in range(4):
        ax = axes[idx]
        hgal_name = hgal_sites[idx]
        hsap_name = hsap_sites[idx]
        title = f"Hgal {hgal_name} ↔ Hsap {hsap_name}"
        
        means = []
        stds = []
        labels = []
        
        for name, summary, color in zip(names, summaries, colors):
            first_rep = list(summary.values())[0]
            sites = first_rep.get("active_sites", {})
            
            if hgal_name in sites:
                site_name = hgal_name
            elif hsap_name in sites:
                site_name = hsap_name
            else:
                continue
            
            vals = []
            for rep in summary.values():
                vals.extend(rep["active_sites"][site_name]["values"])
            vals = np.array(vals)
            
            ax.hist(vals, bins=30, alpha=0.4, color=color,
                    label=f"{name} (μ={vals.mean():.1f}Å)", density=True)
            means.append(vals.mean())
            stds.append(vals.std())
            labels.append(name)
        
        ax.set_xlabel("Distance to nearest TRIM41 CA (Å)")
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.legend(fontsize=7)
        results[title] = {"means": means, "stds": stds, "labels": labels}
    
    fig.suptitle("Active Site Distances — Four Systems (homolog-mapped)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "four_systems_active_sites.png", dpi=200)
    plt.close(fig)
    print("  Saved: four_systems_active_sites.png")
    return results


def summarize_contacts(summaries, names):
    """Aggregate top contacts per system."""
    def aggregate(summary):
        all_occ = {}
        total = 0
        for rep in summary.values():
            n = rep["n_frames"]
            total += n
            for key, occ in rep["occupancy"].items():
                all_occ[key] = all_occ.get(key, 0) + occ * n
        for k in all_occ:
            all_occ[k] /= total
        return all_occ
    
    results = {}
    for name, summary in zip(names, summaries):
        occ = aggregate(summary)
        top = sorted(occ.items(), key=lambda x: x[1], reverse=True)[:10]
        results[name] = top
        print(f"\n  Top contacts — {name}:")
        for key, v in top[:5]:
            t, c = key.split("_")
            print(f"    TRIM41-{t} ↔ cGAS-{c}: {v:.3f}")
    return results


def main():
    parser = argparse.ArgumentParser(description="Four-system MD comparison")
    parser.add_argument("--summaries", nargs=4, required=True,
                        help="Four JSON summary files in order: Hgal_WT, Hgal_4mut_rev, Hsap_WT, Hsap_4mut")
    parser.add_argument("--names", nargs=4,
                        default=["Hgal_WT", "Hgal_4mut_rev", "Hsap_WT", "Hsap_4mut"])
    parser.add_argument("--outdir", default="data/analysis/final_200ns")
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    summaries = [load_summary(p) for p in args.summaries]
    names = args.names
    
    print(f"[{datetime.now()}] Four-system comparison")
    print(f"  Systems: {', '.join(names)}")
    for name, s in zip(names, summaries):
        n_reps = len(s)
        n_frames = sum(r["n_frames"] for r in s.values())
        print(f"    {name}: {n_reps} replica(s), {n_frames} total frames")
    
    # 1. RMSD
    print("\n[1/5] RMSD comparison...")
    rmsd_stats = plot_four_system_rmsd(summaries, names, outdir)
    
    # 2. COM
    print("\n[2/5] COM distance comparison...")
    com_stats = plot_four_system_com(summaries, names, outdir)
    
    # 3. ΔRMSF (intra-species)
    print("\n[3/5] ΔRMSF (intra-species)...")
    drmsf_hgal = plot_intrinsic_rmsf_delta(
        summaries[0], summaries[1], names[0], names[1], outdir)
    drmsf_hsap = plot_intrinsic_rmsf_delta(
        summaries[2], summaries[3], names[2], names[3], outdir)
    
    # 4. Active sites
    print("\n[4/5] Active site distances...")
    active_results = compare_active_sites(summaries, names, outdir)
    
    # 5. Contacts
    print("\n[5/5] Interface contacts summary...")
    contact_results = summarize_contacts(summaries, names)
    
    # Save master summary
    master = {
        "timestamp": str(datetime.now()),
        "systems": names,
        "rmsd_stats": rmsd_stats,
        "com_stats": com_stats,
        "drmsf_hgal": drmsf_hgal,
        "drmsf_hsap": drmsf_hsap,
        "active_sites": active_results,
        "top_contacts": contact_results,
    }
    
    summary_path = outdir / "four_systems_comparison_summary.json"
    with open(summary_path, "w") as f:
        json.dump(master, f, indent=2)
    print(f"\nSaved master summary: {summary_path}")
    print(f"[{datetime.now()}] Done.")


if __name__ == "__main__":
    main()
