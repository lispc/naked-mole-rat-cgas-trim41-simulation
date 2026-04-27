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

FIXED (2026-04-27): Replaced naive Welch t-test on correlated MD frames
with effective-sample-size t-test using autocorrelation time estimation.
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


def effective_sample_size(data, max_lag=None):
    """Estimate effective sample size using integrated autocorrelation time.
    
    Uses Geyer's initial positive sequence estimator on the normalized
    autocorrelation function.  For MD time series the adjacent frames are
    highly correlated, so the effective N is much smaller than len(data).
    """
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
    """Welch-style t-test for correlated time series.
    
    Computes autocorrelation time for each series, derives effective
    sample size (N_eff = N / tau), and performs Welch's t-test on the
    effective samples.
    
    Returns:
        t, p, n_eff_a, n_eff_b, tau_a, tau_b
    """
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
    
    # Welch-Satterthwaite df with effective N
    num = (var_a / n_eff_a + var_b / n_eff_b) ** 2
    den = (var_a / n_eff_a) ** 2 / max(n_eff_a - 1, 1) + (var_b / n_eff_b) ** 2 / max(n_eff_b - 1, 1)
    df = num / den if den > 0 else 1.0
    
    p = 2 * (1 - stats.t.cdf(abs(t), df))
    return float(t), float(p), n_eff_a, n_eff_b, tau_a, tau_b


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
    
    t, p, n_eff_a, n_eff_b, tau_a, tau_b = correlated_ttest(all_rmsd_a, all_rmsd_b)
    fig.suptitle(
        f"RMSD: {name_a} vs {name_b}\n"
        f"t={t:.2f}, p={p:.3e}  (N_eff={n_eff_a:.0f}/{n_eff_b:.0f}, τ={tau_a:.1f}/{tau_b:.1f} frames)",
        fontsize=10
    )
    fig.tight_layout()
    fig.savefig(outdir / f"compare_{name_a}_vs_{name_b}_rmsd.png", dpi=200)
    plt.close(fig)
    print(f"  Saved RMSD comparison (p={p:.3e}, N_eff={n_eff_a:.0f}/{n_eff_b:.0f})")
    return {
        "rmsd_t": float(t), "rmsd_p": float(p),
        "rmsd_n_eff_a": float(n_eff_a), "rmsd_n_eff_b": float(n_eff_b),
        "rmsd_tau_a": float(tau_a), "rmsd_tau_b": float(tau_b),
    }


def compare_rmsf(data_a, data_b, name_a, name_b, outdir):
    """RMSF difference map."""
    first_a = list(data_a.values())[0]
    first_b = list(data_b.values())[0]
    
    resids_a = np.array(first_a["rmsf_resids"])
    resids_b = np.array(first_b["rmsf_resids"])
    
    if len(resids_a) == len(resids_b) and np.allclose(resids_a, resids_b):
        resids = resids_a
        drmsf = np.array(first_b["rmsf_vals"]) - np.array(first_a["rmsf_vals"])
    else:
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
    
    sig_idx = np.where(np.abs(drmsf) > 1.0)[0]
    print(f"    Residues with |ΔRMSF| > 1.0 Å: {len(sig_idx)}")
    return {"drmsf": drmsf.tolist(), "resids": resids.tolist()}


def compare_contacts(data_a, data_b, name_a, name_b, outdir):
    """Compare interface contacts between systems."""
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
    
    t, p, n_eff_a, n_eff_b, tau_a, tau_b = correlated_ttest(all_com_a, all_com_b)
    print(f"  Saved COM comparison (p={p:.3e}, N_eff={n_eff_a:.0f}/{n_eff_b:.0f})")
    return {
        "com_t": float(t), "com_p": float(p),
        "com_n_eff_a": float(n_eff_a), "com_n_eff_b": float(n_eff_b),
        "com_tau_a": float(tau_a), "com_tau_b": float(tau_b),
    }


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
