#!/usr/bin/env python3
"""
WHAM analysis: 1D and 2D PMF from umbrella sampling trajectories.

1D PMF: K315 NZ → Ub G76 C distance
2D PMF: K315 distance + attack angle (NZ-CE vs CE→UbG76)

Usage:
  python wham_pmf.py --system WT --windows 12,13,14,15,16,17,18,19,20,21,22
  python wham_pmf.py --system 4mut --windows 13,14,15,16,17,18,19,20,21,22
"""
import sys, argparse, numpy as np, MDAnalysis as mda, warnings; warnings.filterwarnings("ignore")
from scipy import optimize
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
USDIR = BASE / "data/md_runs/quaternary_full"

# Atom indices (FULL quaternary)
K315_NZ = 11306
K315_CE = 11304
UB_G76_C = 5964


def wham_1d(rc_values, window_centers, k_restraint, T=300):
    """
    Simple 1D WHAM for harmonic umbrella restraints.

    Parameters:
      rc_values: list of arrays, each array is the reaction coordinate for one window
      window_centers: list of target distances (Å) for each window
      k_restraint: force constant in kJ/mol/nm²
      T: temperature in K

    Returns:
      bins: bin centers
      pmf: free energy in kcal/mol
      pmf_err: bootstrap error estimate
    """
    kT = 0.008314 * T / 4.184  # kcal/mol

    # Convert k from kJ/mol/nm² to kcal/mol/Å²
    k_kcal = k_restraint * 0.239 / 100  # kJ→kcal, nm→Å

    # Combine all data
    all_rc = np.concatenate(rc_values)
    rc_min, rc_max = all_rc.min(), all_rc.max()
    bins = np.linspace(rc_min - 0.5, rc_max + 0.5, 50)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins = len(bin_centers)

    # Build histogram counts
    N_w = np.array([len(rc) for rc in rc_values])
    n_iw = np.zeros((n_bins, len(rc_values)))
    for w in range(len(rc_values)):
        n_iw[:, w], _ = np.histogram(rc_values[w], bins=bins)

    # Iterative WHAM
    F = np.zeros(len(rc_values))  # window free energies
    for iteration in range(1000):
        F_old = F.copy()

        # Unbiased probability
        denominator = np.zeros(n_bins)
        for w in range(len(rc_values)):
            bias = 0.5 * k_kcal * (bin_centers - window_centers[w]) ** 2
            denominator += N_w[w] * np.exp(-(F[w] + bias) / kT)

        # New F
        for w in range(len(rc_values)):
            bias_w = 0.5 * k_kcal * (bin_centers - window_centers[w]) ** 2
            numerator = np.sum(n_iw[:, w] / np.maximum(denominator, 1e-30) * np.exp(-bias_w / kT))
            F[w] = -kT * np.log(max(numerator, 1e-30))

        if np.max(np.abs(F - F_old)) < 1e-7:
            print(f"  WHAM converged in {iteration+1} iterations")
            break

    # Compute PMF
    denominator_final = np.zeros(n_bins)
    for w in range(len(rc_values)):
        bias = 0.5 * k_kcal * (bin_centers - window_centers[w]) ** 2
        denominator_final += N_w[w] * np.exp(-(F[w] + bias) / kT)

    pmf = -kT * np.log(np.maximum(denominator_final, 1e-30))
    pmf -= pmf.min()

    # Bootstrap error
    n_bootstrap = 50
    pmf_samples = np.zeros((n_bootstrap, n_bins))
    for b in range(n_bootstrap):
        rc_boot = []
        for w in range(len(rc_values)):
            idx = np.random.choice(len(rc_values[w]), len(rc_values[w]), replace=True)
            rc_boot.append(rc_values[w][idx])
        N_w_b = [len(r) for r in rc_boot]
        n_iw_b = np.zeros((n_bins, len(rc_values)))
        for w in range(len(rc_values)):
            n_iw_b[:, w], _ = np.histogram(rc_boot[w], bins=bins)
        F_b = F.copy()
        for _ in range(50):
            F_b_old = F_b.copy()
            denom_b = np.zeros(n_bins)
            for w in range(len(rc_values)):
                b_bias = 0.5 * k_kcal * (bin_centers - window_centers[w]) ** 2
                denom_b += N_w_b[w] * np.exp(-(F_b[w] + b_bias) / kT)
            for w in range(len(rc_values)):
                b_bias_w = 0.5 * k_kcal * (bin_centers - window_centers[w]) ** 2
                num = np.sum(n_iw_b[:, w] / np.maximum(denom_b, 1e-30) * np.exp(-b_bias_w / kT))
                F_b[w] = -kT * np.log(max(num, 1e-30))
            if np.max(np.abs(F_b - F_b_old)) < 1e-6:
                break
        denom_bf = np.zeros(n_bins)
        for w in range(len(rc_values)):
            b_bias = 0.5 * k_kcal * (bin_centers - window_centers[w]) ** 2
            denom_bf += N_w_b[w] * np.exp(-(F_b[w] + b_bias) / kT)
        pmf_b = -kT * np.log(np.maximum(denom_bf, 1e-30))
        pmf_samples[b] = pmf_b - pmf_b.min()

    pmf_err = pmf_samples.std(axis=0)
    return bin_centers, pmf, pmf_err


def wham_2d(rc1_values, rc2_values, window_centers, k_restraint, T=300,
             n_bins1=25, n_bins2=25):
    """2D WHAM: first dimension has umbrella restraint, second is unbiased."""
    kT = 0.008314 * T / 4.184
    k_kcal = k_restraint * 0.239 / 100

    all_rc1 = np.concatenate(rc1_values)
    all_rc2 = np.concatenate(rc2_values)

    bins1 = np.linspace(all_rc1.min() - 0.5, all_rc1.max() + 0.5, n_bins1 + 1)
    bins2 = np.linspace(0, 180, n_bins2 + 1)
    cents1 = (bins1[:-1] + bins1[1:]) / 2
    cents2 = (bins2[:-1] + bins2[1:]) / 2

    # Build 2D histograms
    N_w = np.array([len(r) for r in rc1_values])
    n_iw = np.zeros((n_bins1, n_bins2, len(rc1_values)))
    for w in range(len(rc1_values)):
        h, _, _ = np.histogram2d(rc1_values[w], rc2_values[w], bins=[bins1, bins2])
        n_iw[:, :, w] = h

    # Iterative WHAM
    F = np.zeros(len(rc1_values))
    for iteration in range(500):
        F_old = F.copy()
        denominator = np.zeros((n_bins1, n_bins2))
        for w in range(len(rc1_values)):
            bias_1d = 0.5 * k_kcal * (cents1[:, None] - window_centers[w]) ** 2
            denominator += N_w[w] * np.exp(-(F[w] + bias_1d) / kT)

        for w in range(len(rc1_values)):
            bias_w = 0.5 * k_kcal * (cents1[:, None] - window_centers[w]) ** 2
            ratio = n_iw[:, :, w] / np.maximum(denominator, 1e-30) * np.exp(-bias_w / kT)
            numerator = np.sum(ratio)
            F[w] = -kT * np.log(max(numerator, 1e-30))

        if np.max(np.abs(F - F_old)) < 1e-6:
            print(f"  2D WHAM converged in {iteration+1} iterations")
            break

    # PMF
    denominator_final = np.zeros((n_bins1, n_bins2))
    for w in range(len(rc1_values)):
        bias = 0.5 * k_kcal * (cents1[:, None] - window_centers[w]) ** 2
        denominator_final += N_w[w] * np.exp(-(F[w] + bias) / kT)

    pmf_2d = -kT * np.log(np.maximum(denominator_final, 1e-30))
    pmf_2d -= pmf_2d.min()
    return cents1, cents2, pmf_2d


def load_window(label, w, stride=5, equil_fraction=0.3):
    """Load equilibrated RC values from a US window."""
    dcd = USDIR / f"us_{label}/window_{w}A/us_{label}_win{w}.dcd"
    if not dcd.exists():
        return None, None, None

    u = mda.Universe(str(dcd), str(dcd))
    n = len(u.trajectory)
    start = int(n * equil_fraction)

    rc1, rc2 = [], []
    for i, ts in enumerate(u.trajectory):
        if i < start or i % stride != 0:
            continue
        p = ts.positions
        # RC1: K315 NZ → Ub G76 C
        rc1.append(np.linalg.norm(p[K315_NZ] - p[UB_G76_C]))
        # RC2: attack angle
        nz, ce, ub = p[K315_NZ], p[K315_CE], p[UB_G76_C]
        v1, v2 = nz - ce, ub - ce
        cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        rc2.append(np.degrees(np.arccos(np.clip(cos_a, -1, 1))))

    return np.array(rc1), np.array(rc2), n


def plot_pmf(bin_centers, pmf, pmf_err, label, outpath):
    """Plot 1D PMF with error bars."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bin_centers, pmf, 'b-', lw=2, label=label)
    ax.fill_between(bin_centers, pmf - pmf_err, pmf + pmf_err, alpha=0.3, color='b')
    ax.set_xlabel("K315 NZ → Ub G76 C (Å)", fontsize=12)
    ax.set_ylabel("Free Energy (kcal/mol)", fontsize=12)
    ax.set_title(f"PMF: {label}", fontsize=14)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_2d_pmf(cents1, cents2, pmf_2d, label, outpath):
    """Plot 2D PMF as contour + heatmap."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 7))
    pmf_clipped = np.clip(pmf_2d, 0, 15)
    cs = ax.contourf(cents1, cents2, pmf_clipped.T, levels=30, cmap='viridis')
    ax.contour(cents1, cents2, pmf_clipped.T, levels=[1, 2, 3, 5, 8], colors='white', linewidths=0.5, alpha=0.5)
    plt.colorbar(cs, ax=ax, label='Free Energy (kcal/mol)')
    ax.set_xlabel("K315 NZ → Ub G76 C (Å)", fontsize=12)
    ax.set_ylabel("Attack Angle (°)", fontsize=12)
    ax.set_title(f"2D PMF: {label}", fontsize=14)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--system', required=True, choices=['WT', '4mut'])
    p.add_argument('--windows', required=True, help='Comma-separated window list, e.g. 12,13,14,...,22')
    p.add_argument('--k', type=float, default=500, help='Restraint force constant (kJ/mol/nm²)')
    p.add_argument('--stride', type=int, default=5, help='Frame stride')
    p.add_argument('--equil', type=float, default=0.3, help='Equilibration fraction to discard')
    p.add_argument('--plot', action='store_true', default=True)
    a = p.parse_args()

    system = a.system
    windows = [int(w.strip()) for w in a.windows.split(',')]
    label = f"{system} (US, {len(windows)} windows)"

    print(f"WHAM Analysis: {label}")
    print(f"  k = {a.k} kJ/mol/nm², stride = {a.stride}, equil = {a.equil*100:.0f}%")

    # Load all windows
    rc1_all, rc2_all, win_centers = [], [], []
    for w in windows:
        rc1, rc2, n_frames = load_window(system, w, a.stride, a.equil)
        if rc1 is None:
            print(f"  Window {w}Å: NOT FOUND, skipping")
            continue
        rc1_all.append(rc1)
        rc2_all.append(rc2)
        win_centers.append(w)
        print(f"  Window {w:2d}Å: {len(rc1)} frames (from {n_frames} total)")

    if len(rc1_all) < 5:
        print("ERROR: Too few windows for WHAM")
        return

    # 1D WHAM
    print(f"\n  Running 1D WHAM...")
    bins, pmf, pmf_err = wham_1d(rc1_all, win_centers, a.k)

    # Find minimum
    min_idx = np.argmin(pmf)
    print(f"\n  PMF minimum: {bins[min_idx]:.1f} Å (F = {pmf[min_idx]:.2f} kcal/mol)")

    # Report free energy at key distances
    for dist in [12, 15, 18, 20, 22]:
        idx = np.argmin(np.abs(bins - dist))
        print(f"  F({dist:2d}Å) = {pmf[idx]:5.2f} ± {pmf_err[idx]:.2f} kcal/mol")

    # 2D WHAM
    print(f"\n  Running 2D WHAM (distance × angle)...")
    cents1, cents2, pmf_2d = wham_2d(rc1_all, rc2_all, win_centers, a.k)

    # 2D PMF minimum
    min2d_idx = np.unravel_index(np.argmin(pmf_2d), pmf_2d.shape)
    print(f"  2D PMF minimum: K315={cents1[min2d_idx[0]]:.1f} Å, angle={cents2[min2d_idx[1]]:.0f}°")

    # 1D projection: angle at K315 distance minimum
    dist_at_min = cents1[min2d_idx[0]]
    angle_slice = pmf_2d[min2d_idx[0], :]
    print(f"  Angle profile at K315={dist_at_min:.1f}Å: min angle={cents2[np.argmin(angle_slice)]:.0f}°")

    # Plot
    if a.plot:
        outdir = BASE / "data/analysis/pfm"
        outdir.mkdir(parents=True, exist_ok=True)
        plot_pmf(bins, pmf, pmf_err, label, outdir / f"pmf_1d_{system}.png")
        plot_2d_pmf(cents1, cents2, pmf_2d, label, outdir / f"pmf_2d_{system}.png")

    # Save data
    np.savez(outdir / f"pmf_data_{system}.npz",
             bins=bins, pmf=pmf, pmf_err=pmf_err,
             cents1=cents1, cents2=cents2, pmf_2d=pmf_2d)
    print(f"  Data saved: {outdir / f'pmf_data_{system}.npz'}")


if __name__ == '__main__':
    main()
