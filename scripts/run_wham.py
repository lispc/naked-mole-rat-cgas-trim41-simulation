#!/usr/bin/env python3
"""
WHAM (Weighted Histogram Analysis Method) for Umbrella Sampling PMF reconstruction.
Reads CV trajectories from multiple US windows and computes the Potential of Mean Force.
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


# Energy units: kJ/mol throughout (matches OpenMM)
KB = 0.008314462618  # kJ/mol/K
T = 310.0  # K
BETA = 1.0 / (KB * T)


def read_cv(path, skip=0):
    """Read CV values from file, skipping header and initial frames."""
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                vals.append(float(parts[2]))
    vals = np.array(vals)
    if skip > 0:
        vals = vals[skip:]
    return vals


def run_wham(window_data, bin_edges, max_iter=10000, tol=1e-6):
    """
    window_data: list of dicts with keys 'cv' (array), 'center', 'k'
    bin_edges: array of bin edges for histogramming
    Returns: bin_centers, pmf (kcal/mol), convergence info
    """
    n_windows = len(window_data)
    n_bins = len(bin_edges) - 1
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Histogram each window
    N = np.zeros((n_windows, n_bins))
    for i, wd in enumerate(window_data):
        counts, _ = np.histogram(wd["cv"], bins=bin_edges)
        N[i, :] = counts
    
    # Bias potential matrix U[i, j] = 0.5 * k_i * (bin_center_j - center_i)^2
    U = np.zeros((n_windows, n_bins))
    for i, wd in enumerate(window_data):
        U[i, :] = 0.5 * wd["k"] * (bin_centers - wd["center"])**2
    
    # Initialize free energy shifts f_i = 0
    f = np.zeros(n_windows)
    
    # Iterative WHAM
    for iteration in range(max_iter):
        # Denominator for each bin: sum over windows of N_i * exp(-beta*(U_i - f_i))
        denom = np.zeros(n_bins)
        for i in range(n_windows):
            denom += N[i, :] * np.exp(-BETA * (U[i, :] - f[i]))
        denom = np.maximum(denom, 1e-300)  # avoid division by zero
        
        # New f_i
        f_new = np.zeros(n_windows)
        for i in range(n_windows):
            # sum over bins of N_i,j / denom * exp(-beta*(U_i,j - f_i))
            # Actually: exp(f_i) = sum_j [N_total,j / denom * exp(-beta*U_i,j)]
            integrand = np.zeros(n_bins)
            for j in range(n_bins):
                if N[:, j].sum() > 0:
                    integrand[j] = (N[:, j].sum() / denom[j]) * np.exp(-BETA * U[i, j])
            f_new[i] = -KB * T * np.log(np.sum(integrand))
        
        # Normalize f so that f[0] = 0
        f_new = f_new - f_new[0]
        
        delta = np.max(np.abs(f_new - f))
        f = f_new
        
        if delta < tol:
            print(f"  WHAM converged in {iteration+1} iterations (Δ={delta:.2e})")
            break
    else:
        print(f"  WHAM did NOT converge in {max_iter} iterations (Δ={delta:.2e})")
    
    # Compute PMF
    pmf = np.zeros(n_bins)
    for j in range(n_bins):
        if denom[j] > 0:
            pmf[j] = -KB * T * np.log(denom[j])
        else:
            pmf[j] = np.nan
    
    # Shift PMF so that minimum = 0
    pmf = pmf - np.nanmin(pmf)
    
    # Convert to kcal/mol for output
    pmf = pmf / 4.184
    
    return bin_centers, pmf, f, iteration


def estimate_error_by_bootstrap(window_data, bin_edges, n_bootstrap=100):
    """Bootstrap error estimation for PMF."""
    pmf_boot = []
    rng = np.random.default_rng(42)
    
    for b in range(n_bootstrap):
        boot_data = []
        for wd in window_data:
            n = len(wd["cv"])
            idx = rng.integers(0, n, size=n)
            boot_data.append({
                "cv": wd["cv"][idx],
                "center": wd["center"],
                "k": wd["k"],
            })
        bc, pmf, _, _ = run_wham(boot_data, bin_edges, max_iter=5000, tol=1e-5)
        pmf_boot.append(pmf)
    
    pmf_boot = np.array(pmf_boot)
    pmf_mean = np.nanmean(pmf_boot, axis=0)
    pmf_std = np.nanstd(pmf_boot, axis=0)
    return pmf_mean, pmf_std


def plot_pmf(bin_centers, pmf, pmf_std, outpath, title="PMF"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bin_centers, pmf, lw=2, color="#2c3e50")
    ax.fill_between(bin_centers, pmf - pmf_std, pmf + pmf_std, alpha=0.3, color="#3498db")
    ax.set_xlabel("RING → Lys-334 distance (Å)")
    ax.set_ylabel("PMF (kcal/mol)")
    ax.set_title(title)
    ax.axhline(0, color="black", lw=0.5, ls="--")
    # Mark window centers
    ax.set_xlim(bin_centers.min(), bin_centers.max())
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def main():
    parser = argparse.ArgumentParser(description="WHAM PMF from US windows")
    parser.add_argument("--windows", nargs="+", required=True,
                        help="Window specs: name:center:k e.g. w04A:4.0:1000")
    parser.add_argument("--cv-dir", required=True,
                        help="Directory containing CV files (name_cv.dat)")
    parser.add_argument("--bin-min", type=float, default=3.0)
    parser.add_argument("--bin-max", type=float, default=22.0)
    parser.add_argument("--bin-width", type=float, default=0.2)
    parser.add_argument("--skip", type=int, default=0,
                        help="Skip first N frames from each window")
    parser.add_argument("--bootstrap", type=int, default=100)
    parser.add_argument("--outdir", default="data/analysis/final/us_Hsap_WT_Lys334")
    parser.add_argument("--title", default="Hsap_WT Lys-334 US PMF")
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Parse window specs
    window_data = []
    for spec in args.windows:
        parts = spec.split(":")
        name = parts[0]
        center = float(parts[1])
        k = float(parts[2])
        cv_path = Path(args.cv_dir) / f"us_{name}_cv.dat"
        if not cv_path.exists():
            print(f"ERROR: CV file not found: {cv_path}")
            sys.exit(1)
        cv = read_cv(cv_path, skip=args.skip)
        print(f"  {name}: center={center}Å, k={k}, N={len(cv)}, CV range=[{cv.min():.2f}, {cv.max():.2f}]")
        window_data.append({"name": name, "cv": cv, "center": center, "k": k})
    
    # Setup bins
    bin_edges = np.arange(args.bin_min, args.bin_max + args.bin_width, args.bin_width)
    print(f"\nBins: {len(bin_edges)-1} bins, [{args.bin_min}, {args.bin_max}] Å, width={args.bin_width} Å")
    
    # Run WHAM
    print("\n[1/3] Running WHAM...")
    bin_centers, pmf, f, n_iter = run_wham(window_data, bin_edges)
    
    # Bootstrap error
    print(f"\n[2/3] Bootstrap error estimation (n={args.bootstrap})...")
    _, pmf_std = estimate_error_by_bootstrap(window_data, bin_edges, n_bootstrap=args.bootstrap)
    
    # Plot
    print(f"\n[3/3] Plotting...")
    plot_pmf(bin_centers, pmf, pmf_std, outdir / "wham_pmf.png", title=args.title)
    
    # Save data
    data_path = outdir / "wham_pmf.dat"
    with open(data_path, "w") as f_out:
        f_out.write("#distance_A PMF_kcal_mol PMF_err_kcal_mol\n")
        for x, y, e in zip(bin_centers, pmf, pmf_std):
            if not np.isnan(y):
                f_out.write(f"{x:.3f} {y:.4f} {e:.4f}\n")
    print(f"  Saved: {data_path}")
    
    # Summary stats
    min_idx = np.nanargmin(pmf)
    min_dist = bin_centers[min_idx]
    min_pmf = pmf[min_idx]
    
    # Find FWHM-like width or approximate basin
    threshold = 2.0  # kcal/mol
    below = np.where(pmf < threshold)[0]
    if len(below) > 0:
        basin_min = bin_centers[below[0]]
        basin_max = bin_centers[below[-1]]
    else:
        basin_min = basin_max = min_dist
    
    summary = {
        "timestamp": str(datetime.now()),
        "n_windows": len(window_data),
        "total_frames": sum(len(w["cv"]) for w in window_data),
        "bin_width": args.bin_width,
        "pmf_minimum": {"distance_A": float(min_dist), "pmf_kcal_mol": float(min_pmf)},
        "basin_2kcal": {"min_A": float(basin_min), "max_A": float(basin_max), "width_A": float(basin_max - basin_min)},
        "free_energy_barrier": float(np.nanmax(pmf) - min_pmf),
    }
    
    with open(outdir / "wham_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  PMF minimum: {min_dist:.2f} Å, {min_pmf:.2f} kcal/mol")
    print(f"  2 kcal/mol basin: [{basin_min:.2f}, {basin_max:.2f}] Å (width={basin_max-basin_min:.2f} Å)")
    print(f"  Free energy barrier: {summary['free_energy_barrier']:.2f} kcal/mol")
    print(f"\nDone.")


if __name__ == "__main__":
    main()
