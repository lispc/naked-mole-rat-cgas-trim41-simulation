#!/usr/bin/env python3
"""
US Convergence Check (fast version — no full WHAM).

Three checks:
  1. CV mean stability: split each window into 10 blocks, check if mean drifts
  2. Adjacent window histogram overlap
  3. Half-split KS test: first-half vs second-half CV distribution

Output: data/analysis/us_convergence/
"""
import json
import warnings
from pathlib import Path
from glob import glob
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

BASE = Path(__file__).resolve().parent.parent.parent
OUTDIR = BASE / "data" / "analysis" / "us_convergence"
OUTDIR.mkdir(parents=True, exist_ok=True)

US_DATA = {
    "WT": BASE / "data/analysis/final/us_Hsap_WT_Lys315",
    "4mut": BASE / "data/analysis/final/us_Hsap_4mut_Lys315",
}


def load_cv_data(us_dir):
    """Load CV data, prefer latest version (v2/v3 over v1)."""
    cv_files = sorted(glob(str(us_dir / "us_w*_cv.dat")))
    window_data = {}
    for f in cv_files:
        name = Path(f).stem
        m = re.match(r"us_w(\d+A)(?:_v(\d+))?_cv", name)
        if m:
            wid = m.group(1)
            ver = int(m.group(2)) if m.group(2) else 1
            if wid not in window_data or ver > window_data[wid]["version"]:
                data = np.loadtxt(f)
                window_data[wid] = {"version": ver, "data": data,
                                     "n_frames": len(data)}
    sorted_wids = sorted(window_data.keys(),
                         key=lambda x: int(x.replace("A", "")))
    return {w: window_data[w] for w in sorted_wids}


def check_cv_mean_stability(window_data, n_blocks=10):
    """Check if per-window CV mean stabilizes across blocks."""
    results = {}
    for wid, info in window_data.items():
        cv = info["data"][:, 2]
        n = len(cv)
        block_size = n // n_blocks
        block_means = []
        for b in range(n_blocks):
            start = b * block_size
            end = start + block_size if b < n_blocks - 1 else n
            block_means.append(np.mean(cv[start:end]))
        block_means = np.array(block_means)
        # Check: is the trend slope significantly non-zero?
        x = np.arange(n_blocks)
        slope, intercept, r, p, std_err = stats.linregress(x, block_means)
        # Also check: drift between first and last third
        first_third = np.mean(block_means[:3])
        last_third = np.mean(block_means[-3:])
        drift = last_third - first_third
        results[wid] = {
            "n_blocks": n_blocks,
            "block_means": block_means.tolist(),
            "mean": float(np.mean(cv)),
            "std": float(np.std(cv)),
            "slope_per_block": float(slope),
            "drift_first_to_last": float(drift),
            "drift_significant": abs(drift) > 0.5,  # >0.5 Å drift
        }
    return results


def check_histogram_overlap(window_data):
    """Check CV histogram overlap between adjacent windows."""
    wids = sorted(window_data.keys(),
                   key=lambda x: int(x.replace("A", "")))
    bins_global = np.linspace(0, 35, 80)
    overlaps = []
    for i in range(len(wids) - 1):
        cv1 = window_data[wids[i]]["data"][:, 2]
        cv2 = window_data[wids[i + 1]]["data"][:, 2]
        h1, _ = np.histogram(cv1, bins=bins_global, density=True)
        h2, _ = np.histogram(cv2, bins=bins_global, density=True)
        min_sum = np.sum(np.minimum(h1, h2))
        max_sum = np.sum(np.maximum(h1, h2))
        overlap = min_sum / max_sum if max_sum > 0 else 0
        overlaps.append((wids[i], wids[i + 1], float(overlap)))
    return overlaps


def half_split_ks_test(window_data):
    """KS test: first-half vs second-half CV distribution per window."""
    results = []
    for wid in sorted(window_data.keys(),
                       key=lambda x: int(x.replace("A", ""))):
        cv = window_data[wid]["data"][:, 2]
        mid = len(cv) // 2
        ks_stat, ks_p = stats.ks_2samp(cv[:mid], cv[mid:])
        results.append({
            "window": wid,
            "ks_statistic": float(ks_stat),
            "ks_p_value": float(ks_p),
            "significant_drift": ks_p < 0.01,
        })
    return results


def main():
    print("=" * 60)
    print("US Convergence Check (fast)")
    print("=" * 60)

    all_results = {}

    for name, us_dir in US_DATA.items():
        print(f"\n{'─' * 40}")
        print(f"System: {name}")
        print(f"{'─' * 40}")

        window_data = load_cv_data(us_dir)
        wids = sorted(window_data.keys(),
                       key=lambda x: int(x.replace("A", "")))
        total_frames = sum(w["n_frames"] for w in window_data.values())
        print(f"Windows: {wids}")
        print(f"Total frames: {total_frames}")
        for w in wids:
            print(f"  {w}: {window_data[w]['n_frames']} frames, "
                  f"mean={np.mean(window_data[w]['data'][:, 2]):.1f} Å, "
                  f"std={np.std(window_data[w]['data'][:, 2]):.2f} Å")

        # 1. CV mean stability
        print("\n[1] CV mean stability (block analysis):")
        stability = check_cv_mean_stability(window_data)
        unstable = []
        for wid, r in stability.items():
            flag = ""
            if r["drift_significant"]:
                flag = f" *** DRIFT {r['drift_first_to_last']:+.2f} Å"
                unstable.append(wid)
            print(f"  {wid}: mean={r['mean']:.2f} Å, slope={r['slope_per_block']:.4f}/block, "
                  f"drift={r['drift_first_to_last']:+.3f} Å{flag}")

        # 2. Histogram overlap
        print("\n[2] Adjacent window histogram overlap:")
        overlaps = check_histogram_overlap(window_data)
        low_overlap = []
        for w1, w2, ov in overlaps:
            status = "OK" if ov > 0.05 else "LOW"
            print(f"  {w1} ↔ {w2}: {ov:.3f} [{status}]")
            if ov < 0.05:
                low_overlap.append((w1, w2))

        # 3. KS test
        print("\n[3] Half-split KS test:")
        ks_results = half_split_ks_test(window_data)
        drifted = []
        for r in ks_results:
            flag = " *** DRIFT" if r["significant_drift"] else ""
            print(f"  {r['window']}: KS={r['ks_statistic']:.3f}, p={r['ks_p_value']:.4f}{flag}")
            if r["significant_drift"]:
                drifted.append(r["window"])

        # ── Plots ──────────────────────────────────────────────────────
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Panel 1: All window CV distributions
        ax = axes[0, 0]
        cv_list = [window_data[w]["data"][:, 2] for w in wids]
        bins_g = np.linspace(0, 35, 80)
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(wids)))
        for wid, cv, c in zip(wids, cv_list, colors):
            ax.hist(cv, bins=bins_g, alpha=0.4, color=c, label=wid, density=True)
        ax.set_xlabel("CV (Å)")
        ax.set_ylabel("Density")
        ax.set_title(f"{name}: Window CV distributions")
        ax.legend(fontsize=7, ncol=2)

        # Panel 2: Overlap between adjacent windows
        ax = axes[0, 1]
        ov_vals = [o[2] for o in overlaps]
        pair_labels = [f"{w1[:4]}\n{w2[:4]}" for w1, w2, _ in overlaps]
        x = np.arange(len(ov_vals))
        bars = ax.bar(x, ov_vals,
                      color=["green" if o > 0.05 else "red" for o in ov_vals])
        ax.axhline(0.05, color="gray", ls="--", alpha=0.7, label="0.05 threshold")
        ax.set_xticks(x)
        ax.set_xticklabels(pair_labels, fontsize=7)
        ax.set_ylabel("Histogram overlap")
        ax.set_title("Adjacent window overlap")
        ax.legend(fontsize=8)

        # Panel 3: Block mean stability (per window)
        ax = axes[1, 0]
        for wid, r in stability.items():
            ax.plot(np.arange(r["n_blocks"]) + 1, r["block_means"],
                    "o-", ms=3, lw=0.8, label=wid)
        ax.set_xlabel("Block index")
        ax.set_ylabel("CV mean (Å)")
        ax.set_title(f"{name}: Per-window CV mean stability")
        ax.legend(fontsize=5, ncol=3, loc="best")

        # Panel 4: KS statistics
        ax = axes[1, 1]
        ks_stats = [r["ks_statistic"] for r in ks_results]
        ks_ps = [r["ks_p_value"] for r in ks_results]
        x_w = np.arange(len(wids))
        ax.bar(x_w, ks_stats,
               color=["red" if p < 0.01 else "steelblue" for p in ks_ps])
        ax.set_xticks(x_w)
        ax.set_xticklabels(wids, rotation=45, fontsize=7)
        ax.set_ylabel("KS statistic")
        ax.set_title(f"{name}: Half-split KS test\n(red = p<0.01, significant drift)")
        ax.axhline(0.1, color="gray", ls="--", alpha=0.5)

        fig.tight_layout()
        fig.savefig(OUTDIR / f"us_convergence_{name}.png", dpi=300)
        plt.close(fig)
        print(f"\nSaved: {OUTDIR / f'us_convergence_{name}.png'}")

        all_results[name] = {
            "n_windows": len(wids),
            "windows": wids,
            "total_frames": total_frames,
            "unstable_windows": unstable,
            "low_overlap_pairs": low_overlap,
            "drifted_windows": drifted,
        }

    # ── Summary ────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    for name in ["WT", "4mut"]:
        r = all_results[name]
        issues = []
        if r["unstable_windows"]:
            issues.append(f"CV mean drift: {r['unstable_windows']}")
        if r["low_overlap_pairs"]:
            issues.append(f"Low overlap: {r['low_overlap_pairs']}")
        if r["drifted_windows"]:
            issues.append(f"KS drift: {r['drifted_windows']}")
        if issues:
            print(f"\n{name}: ISSUES FOUND")
            for issue in issues:
                print(f"  - {issue}")
            print(f"  → Consider extending US to 20 ns/window")
        else:
            print(f"\n{name}: ALL CHECKS PASSED")
            print(f"  → 10 ns/window appears sufficient")

    with open(OUTDIR / "us_convergence_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved: {OUTDIR / 'us_convergence_results.json'}")
    print("Done.")


if __name__ == "__main__":
    main()
