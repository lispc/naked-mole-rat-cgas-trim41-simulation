#!/usr/bin/env python3
"""
Route A: Test the "4mut rigidifies SPRY-cGAS interface" hypothesis.

Metrics computed from existing binary MD (3×200ns WT, 3×200ns 4mut):
  1. Cα distance autocorrelation functions for interface residue pairs
     → decorrelation time τ (1/e decay) quantifies breathing timescale
  2. Interface contact persistence lifetimes
  3. Per-residue Cα position autocorrelation (local rigidity)

If the Discussion hypothesis is correct, 4mut interface τ should be
significantly longer than WT.

Output: data/analysis/rigidity_test/
"""
import json
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import distance_array
from scipy.optimize import curve_fit
from scipy import stats

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

# ── Paths ───────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent.parent
OUTDIR = BASE / "data" / "analysis" / "rigidity_test"
OUTDIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "WT": {
        "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        "trajs": [
            BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        ],
    },
    "4mut": {
        "prmtop": BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
        "trajs": [
            BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
        ],
    },
}

# TRIM41 SPRY = 1-218, cGAS = 219-541 in the binary PDB
TRIM_SEL = "resid 1-218 and name CA"
CGAS_SEL = "resid 219-541 and name CA"
TRIM_HEAVY = "resid 1-218 and (name N* or name O* or name C*)"
CGAS_HEAVY = "resid 219-541 and (name N* or name O* or name C*)"
INTERFACE_CUTOFF = 5.0  # Å, heavy-atom contact
STRIDE = 100  # ps per frame in the DCD trajectory
N_LAGS = 1001  # autocorrelation lags (100 ns at 100 ps stride)


def build_universe(prmtop, traj):
    u = mda.Universe(str(prmtop), str(traj))
    return u


def get_interface_residues(u, cutoff=INTERFACE_CUTOFF):
    """Identify interface residues from first frame."""
    cgas = u.select_atoms(CGAS_HEAVY)
    trim = u.select_atoms(TRIM_HEAVY)
    u.trajectory[0]
    dist_mat = distance_array(cgas.positions, trim.positions)
    close = dist_mat < cutoff
    cgas_iface = set(cgas.resids[np.unique(np.where(close)[0])])
    trim_iface = set(trim.resids[np.unique(np.where(close)[1])])
    return sorted(cgas_iface), sorted(trim_iface)


def compute_distance_acf(dist_ts, max_lag):
    """
    Compute autocorrelation function of a distance time series.
    C(t) = <(d(τ)-μ)(d(τ+t)-μ)> / <(d-μ)²>
    """
    n = len(dist_ts)
    d = dist_ts - dist_ts.mean()
    var = np.var(dist_ts)
    if var < 1e-12:
        return np.ones(max_lag)
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        if lag == 0:
            acf[lag] = 1.0
        else:
            acf[lag] = np.mean(d[:n - lag] * d[lag:]) / var
    return acf


def exp_decay(t, tau):
    return np.exp(-t / tau)


def fit_decorrelation_time(acf, stride_ps, max_lag_ns=50):
    """Fit exponential decay to ACF, return τ in ns and R²."""
    max_lag_frames = min(int(max_lag_ns * 1000 / stride_ps), len(acf))
    t_ns = np.arange(max_lag_frames) * stride_ps / 1000.0
    y = acf[:max_lag_frames]
    # Only fit positive values and up to where ACF crosses zero
    pos_mask = y > 0.01
    if np.sum(pos_mask) < 5:
        return np.nan, np.nan
    t_fit = t_ns[pos_mask]
    y_fit = y[pos_mask]
    try:
        popt, pcov = curve_fit(exp_decay, t_fit, y_fit, p0=[1.0], bounds=(0.1, 500))
        tau = popt[0]
        y_pred = exp_decay(t_fit, tau)
        ss_res = np.sum((y_fit - y_pred) ** 2)
        ss_tot = np.sum((y_fit - y_fit.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return tau, r2
    except Exception:
        return np.nan, np.nan


def analyze_system(name, prmtop, traj_paths):
    """Full rigidity analysis on one system across all replicas."""
    all_tau = []
    all_pair_tau = {}
    all_lifetimes = []
    rep_results = []

    for rep_idx, traj_path in enumerate(traj_paths):
        print(f"  [{name}] Analyzing replica {rep_idx + 1}...")
        u = build_universe(prmtop, traj_path)
        cgas_iface, trim_iface = get_interface_residues(u)
        print(f"    cGAS interface residues: {cgas_iface}")
        print(f"    TRIM41 interface residues: {trim_iface}")

        # Align trajectory to first frame (CA only)
        ca = u.select_atoms("name CA")
        u.trajectory[0]
        aligner = align.AlignTraj(u, u, select="name CA", in_memory=True).run()

        n_frames = len(u.trajectory)
        max_lag = min(N_LAGS, n_frames // 2)

        # ── Efficient: extract all Cα positions in one pass ──────────
        print(f"    Extracting {len(cgas_iface)}+{len(trim_iface)} CA positions over {n_frames} frames...")
        # Build position arrays: (n_frames, n_cgas, 3) and (n_frames, n_trim, 3)
        cgas_ca_pos = np.zeros((n_frames, len(cgas_iface), 3), dtype=np.float32)
        trim_ca_pos = np.zeros((n_frames, len(trim_iface), 3), dtype=np.float32)

        cgas_sel_str = f"resid {' '.join(map(str, cgas_iface))} and name CA"
        trim_sel_str = f"resid {' '.join(map(str, trim_iface))} and name CA"
        cgas_atoms = u.select_atoms(cgas_sel_str)
        trim_atoms = u.select_atoms(trim_sel_str)

        # Also for heavy-atom contact tracking
        cgas_heavy_sel = f"resid {' '.join(map(str, cgas_iface))} and (name N* or name O* or name C*)"
        trim_heavy_sel = f"resid {' '.join(map(str, trim_iface))} and (name N* or name O* or name C*)"
        cgas_heavy = u.select_atoms(cgas_heavy_sel)
        trim_heavy = u.select_atoms(trim_heavy_sel)

        for i, ts in enumerate(u.trajectory):
            cgas_ca_pos[i] = cgas_atoms.positions
            trim_ca_pos[i] = trim_atoms.positions

        # ── 1. Pairwise Cα distance ACF ────────────────────────────
        print(f"    Computing distance ACF for {len(cgas_iface) * len(trim_iface)} pairs...")
        pair_count = 0
        rep_tau = []

        for ci, cg in enumerate(cgas_iface):
            for tj, tr in enumerate(trim_iface):
                # Distance time series for this pair
                delta = cgas_ca_pos[:, ci, :] - trim_ca_pos[:, tj, :]
                dist_ts = np.sqrt(np.sum(delta ** 2, axis=1))

                # Compute ACF and fit
                acf = compute_distance_acf(dist_ts, max_lag)
                tau, r2 = fit_decorrelation_time(acf, STRIDE)
                if not np.isnan(tau) and r2 > 0.3:
                    rep_tau.append(tau)
                    key = (cg, tr)
                    if key not in all_pair_tau:
                        all_pair_tau[key] = []
                    all_pair_tau[key].append(tau)
                pair_count += 1
                if pair_count % 50 == 0:
                    print(f"      {pair_count} pairs processed...")

        all_tau.extend(rep_tau)
        print(f"    Valid τ fits: {len(rep_tau)} / {pair_count}")

        # ── 2. Contact lifetimes (vectorized over frames) ──────────
        print(f"    Computing contact lifetimes...")
        pair_lifetimes = compute_contact_lifetimes_vectorized(
            u, cgas_iface, trim_iface, cgas_heavy, trim_heavy
        )
        all_lifetimes.extend(pair_lifetimes.values())

        rep_results.append({
            "rep": rep_idx + 1,
            "n_interface_pairs": len(rep_tau),
            "tau_mean_ns": np.mean(rep_tau) if rep_tau else np.nan,
            "tau_median_ns": np.median(rep_tau) if rep_tau else np.nan,
            "tau_std_ns": np.std(rep_tau) if rep_tau else np.nan,
            "lifetime_mean_ns": np.mean(list(pair_lifetimes.values())) if pair_lifetimes else np.nan,
        })

    return {
        "all_tau": all_tau,
        "pair_tau": all_pair_tau,
        "all_lifetimes": all_lifetimes,
        "rep_results": rep_results,
    }


def compute_contact_lifetimes_vectorized(u, cgas_iface, trim_iface,
                                         cgas_heavy, trim_heavy):
    """Vectorized contact detection and lifetime computation."""
    n_frames = len(u.trajectory)
    n_cgas_h = len(cgas_heavy) // len(cgas_iface)  # heavy atoms per residue
    n_trim_h = len(trim_heavy) // len(trim_iface)

    lifetimes = {}
    for ci, cg in enumerate(cgas_iface):
        c_sel = cgas_heavy.select_atoms(f"resid {cg}")
        for tj, tr in enumerate(trim_iface):
            t_sel = trim_heavy.select_atoms(f"resid {tr}")
            contact = np.zeros(n_frames, dtype=bool)
            for i, ts in enumerate(u.trajectory):
                d = distance_array(c_sel.positions, t_sel.positions)
                contact[i] = np.min(d) < INTERFACE_CUTOFF

            # Compute average contact event duration
            events = []
            in_event = False
            event_start = 0
            for i, c in enumerate(contact):
                if c and not in_event:
                    in_event = True
                    event_start = i
                elif not c and in_event:
                    in_event = False
                    events.append(i - event_start)
            if in_event:
                events.append(n_frames - event_start)

            lifetimes[(cg, tr)] = np.mean(events) * STRIDE / 1000.0 if events else 0.0

    return lifetimes


def main():
    print("=" * 60)
    print("Route A: Interface Rigidity Hypothesis Test")
    print("=" * 60)

    results = {}
    for name, cfg in SYSTEMS.items():
        print(f"\n{'─' * 40}")
        print(f"System: {name}")
        print(f"{'─' * 40}")
        results[name] = analyze_system(name, cfg["prmtop"], cfg["trajs"])

    # ── Statistical comparison ──────────────────────────────────────────
    wt_tau = results["WT"]["all_tau"]
    mut_tau = results["4mut"]["all_tau"]
    wt_life = results["WT"]["all_lifetimes"]
    mut_life = results["4mut"]["all_lifetimes"]

    print(f"\n{'=' * 60}")
    print("RESULTS")
    print(f"{'=' * 60}")

    print(f"\nWT  interface pairs with valid τ: {len(wt_tau)}")
    print(f"4mut interface pairs with valid τ: {len(mut_tau)}")
    print(f"WT  τ mean/median/std: {np.mean(wt_tau):.2f} / {np.median(wt_tau):.2f} / {np.std(wt_tau):.2f} ns")
    print(f"4mut τ mean/median/std: {np.mean(mut_tau):.2f} / {np.median(mut_tau):.2f} / {np.std(mut_tau):.2f} ns")

    # Mann-Whitney U test (τ distribution is typically non-normal)
    u_stat, p_value = stats.mannwhitneyu(wt_tau, mut_tau, alternative="two-sided")
    print(f"\nMann-Whitney U test on τ: U={u_stat:.1f}, p={p_value:.4f}")
    # One-sided: 4mut τ > WT τ (rigidification)
    u_stat_1s, p_value_1s = stats.mannwhitneyu(mut_tau, wt_tau, alternative="greater")
    print(f"One-sided (4mut > WT): U={u_stat_1s:.1f}, p={p_value_1s:.4f}")

    # Contact lifetime comparison
    u_life, p_life = stats.mannwhitneyu(wt_life, mut_life, alternative="two-sided")
    print(f"\nContact lifetime Mann-Whitney: U={u_life:.1f}, p={p_life:.4f}")
    print(f"WT  contact lifetime mean: {np.mean(wt_life):.2f} ns")
    print(f"4mut contact lifetime mean: {np.mean(mut_life):.2f} ns")

    # ── Plots ───────────────────────────────────────────────────────────
    # 1. τ distribution histogram
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    max_tau = max(np.percentile(wt_tau, 99), np.percentile(mut_tau, 99))
    bins = np.linspace(0, max_tau, 50)

    axes[0].hist(wt_tau, bins=bins, alpha=0.6, label=f"WT (n={len(wt_tau)})", color="steelblue", density=True)
    axes[0].hist(mut_tau, bins=bins, alpha=0.6, label=f"4mut (n={len(mut_tau)})", color="coral", density=True)
    axes[0].axvline(np.median(wt_tau), color="steelblue", ls="--", lw=1.5)
    axes[0].axvline(np.median(mut_tau), color="coral", ls="--", lw=1.5)
    axes[0].set_xlabel("Decorrelation time τ (ns)")
    axes[0].set_ylabel("Density")
    axes[0].set_title(f"Cα distance ACF decay time\np={p_value:.4f} (two-sided MW)")
    axes[0].legend(fontsize=9)

    # 2. Contact lifetime distribution
    max_life = max(np.percentile(wt_life, 99), np.percentile(mut_life, 99))
    bins_l = np.linspace(0, max_life, 50)
    axes[1].hist(wt_life, bins=bins_l, alpha=0.6, label=f"WT", color="steelblue", density=True)
    axes[1].hist(mut_life, bins=bins_l, alpha=0.6, label=f"4mut", color="coral", density=True)
    axes[1].set_xlabel("Contact lifetime (ns)")
    axes[1].set_ylabel("Density")
    axes[1].set_title(f"Interface contact persistence\np={p_life:.4f}")
    axes[1].legend(fontsize=9)

    # 3. Per-replica summary
    for name in ["WT", "4mut"]:
        reps = results[name]["rep_results"]
        x = [r["rep"] for r in reps]
        y = [r["tau_mean_ns"] for r in reps]
        yerr = [r["tau_std_ns"] for r in reps]
        color = "steelblue" if name == "WT" else "coral"
        axes[2].errorbar([xi + (-0.1 if name == "WT" else 0.1) for xi in x], y, yerr=yerr,
                         fmt="o", capsize=5, color=color, label=name, markersize=8)
    axes[2].set_xlabel("Replica")
    axes[2].set_ylabel("Mean τ (ns)")
    axes[2].set_title("Per-replica decorrelation time")
    axes[2].set_xticks([1, 2, 3])
    axes[2].legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTDIR / "rigidity_summary.png", dpi=300)
    plt.close(fig)
    print(f"\nSaved: {OUTDIR / 'rigidity_summary.png'}")

    # 4. Average ACF curves (pooled)
    fig, ax = plt.subplots(figsize=(7, 4))
    t_ns = np.arange(N_LAGS) * STRIDE / 1000.0

    for name, color in [("WT", "steelblue"), ("4mut", "coral")]:
        pair_tau = results[name]["pair_tau"]
        # Average ACF across all pairs for each system is complex since
        # each pair has different trajectory. Instead, show τ CDF.
        all_t = sorted(results[name]["all_tau"])
        cdf = np.arange(1, len(all_t) + 1) / len(all_t)
        ax.plot(all_t, cdf, lw=2, color=color, label=f"{name} (n={len(all_t)})")

    ax.set_xlabel("Decorrelation time τ (ns)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("CDF of interface residue pair decorrelation times")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "tau_cdf.png", dpi=300)
    plt.close(fig)
    print(f"Saved: {OUTDIR / 'tau_cdf.png'}")

    # ── Export data ─────────────────────────────────────────────────────
    export = {
        "analysis": "Route A: Interface rigidity hypothesis test",
        "method": "Cα distance autocorrelation + exponential decay fit",
        "interface_cutoff_A": INTERFACE_CUTOFF,
        "stride_ps": STRIDE,
        "WT": {
            "n_pairs": len(wt_tau),
            "tau_mean_ns": float(np.mean(wt_tau)),
            "tau_median_ns": float(np.median(wt_tau)),
            "tau_std_ns": float(np.std(wt_tau)),
            "tau_sem_ns": float(np.std(wt_tau) / np.sqrt(len(wt_tau))),
            "lifetime_mean_ns": float(np.mean(wt_life)),
            "lifetime_median_ns": float(np.median(wt_life)),
            "rep_results": results["WT"]["rep_results"],
        },
        "4mut": {
            "n_pairs": len(mut_tau),
            "tau_mean_ns": float(np.mean(mut_tau)),
            "tau_median_ns": float(np.median(mut_tau)),
            "tau_std_ns": float(np.std(mut_tau)),
            "tau_sem_ns": float(np.std(mut_tau) / np.sqrt(len(mut_tau))),
            "lifetime_mean_ns": float(np.mean(mut_life)),
            "lifetime_median_ns": float(np.median(mut_life)),
            "rep_results": results["4mut"]["rep_results"],
        },
        "statistics": {
            "test": "Mann-Whitney U (two-sided)",
            "tau_u_statistic": float(u_stat),
            "tau_p_value": float(p_value),
            "tau_one_sided_4mut_greater_p": float(p_value_1s),
            "lifetime_u_statistic": float(u_life),
            "lifetime_p_value": float(p_life),
        },
    }

    with open(OUTDIR / "rigidity_results.json", "w") as f:
        json.dump(export, f, indent=2, default=str)
    print(f"Saved: {OUTDIR / 'rigidity_results.json'}")

    # ── Scientific interpretation ───────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("INTERPRETATION")
    print(f"{'=' * 60}")
    if p_value < 0.05:
        if np.median(mut_tau) > np.median(wt_tau):
            direction = "LONGER → supports rigidification hypothesis"
        else:
            direction = "SHORTER → contradicts rigidification hypothesis"
        print(f"Significant difference in τ (p={p_value:.4f}): 4mut τ is {direction}")
    else:
        ratio = np.median(mut_tau) / np.median(wt_tau) if np.median(wt_tau) > 0 else float("inf")
        print(f"No significant difference in τ (p={p_value:.4f})")
        print(f"4mut/WT τ ratio: {ratio:.2f}")
        print("→ The 'rigidification' hypothesis is NOT supported by quantitative analysis")
        print("→ Discussion §3's claim that 4mut 'rigidifies the SPRY-cGAS interface'")
        print("  should be revised or removed unless other evidence emerges.")

    print(f"\n{'=' * 60}")
    print("Done.")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
