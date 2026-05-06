#!/usr/bin/env python3
"""
Optimized four-system MD comparison (Hgal_WT vs Hgal_4mut_rev vs Hsap_WT vs Hsap_4mut).
Metrics: COM distance, CA-based interface contacts, per-domain Rg, RMSF.
Avoids heavy all-atom distance arrays; uses CA-only for contacts.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from pathlib import Path
import sys
import json
import warnings
warnings.filterwarnings("ignore")

# Ensure scripts/lib/ is importable from any working directory
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lib.paths import BASE, FOUR_SYSTEM_DIR as OUTDIR  # noqa: E402
OUTDIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "Hgal_WT": {
        "prmtop": BASE / "data/md_runs/Hgal_WT/rep1/Hgal_WT.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hgal_WT/rep1/Hgal_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hgal_WT/rep2/Hgal_WT_rep2_prod_fixed.dcd",
            BASE / "data/md_runs/Hgal_WT/rep2/Hgal_WT_rep2_restart.dcd",
            BASE / "data/md_runs/Hgal_WT/rep3/Hgal_WT_rep3_prod.dcd",
        ],
        "cgas_last": 218,
    },
    "Hgal_4mut_rev": {
        "prmtop": BASE / "data/md_runs/Hgal_4mut_rev/rep1/Hgal_4mut_rev.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hgal_4mut_rev/rep1/Hgal_4mut_rev_rep1_prod.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_prod.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_restart.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep3/Hgal_4mut_rev_rep3_prod.dcd",
        ],
        "cgas_last": 218,
    },
    "Hsap_WT": {
        "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        ],
        "cgas_last": 198,
    },
    "Hsap_4mut": {
        "prmtop": BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
        ],
        "cgas_last": 198,
    },
}

# Analysis parameters
CONTACT_CUTOFF = 8.0  # Angstrom, CA-CA distance
DT_NS = 0.1  # frame spacing in ns


def compute_system(name, info):
    """Compute metrics for one system across all reps."""
    prmtop = info["prmtop"]
    dcd_paths = [str(p) for p in info["dcds"] if p.exists()]
    cgas_last = info["cgas_last"]

    if not dcd_paths:
        raise ValueError(f"No DCD files found for {name}")

    # Load all reps into one Universe
    u = mda.Universe(str(prmtop), *dcd_paths)
    print(f"  [{name}] Loaded {len(dcd_paths)} DCDs, {len(u.trajectory)} frames")

    # Selections
    cgas = u.select_atoms(f"protein and resid 1:{cgas_last} and name CA")
    trim41 = u.select_atoms(f"protein and resid {cgas_last+1}:9999 and name CA")
    cgas_all = u.select_atoms(f"protein and resid 1:{cgas_last}")
    trim41_all = u.select_atoms(f"protein and resid {cgas_last+1}:9999")
    all_ca = u.select_atoms("protein and name CA")

    n_cgas_ca = len(cgas)
    n_trim41_ca = len(trim41)
    n_all_ca = len(all_ca)

    print(f"  [{name}] cGAS CA={n_cgas_ca}, TRIM41 CA={n_trim41_ca}, total CA={n_all_ca}")

    # Pre-allocate arrays
    n_frames = len(u.trajectory)
    com_dist = np.empty(n_frames, dtype=np.float32)
    contacts = np.empty(n_frames, dtype=np.int32)
    rg_cgas = np.empty(n_frames, dtype=np.float32)
    rg_trim41 = np.empty(n_frames, dtype=np.float32)
    rg_total = np.empty(n_frames, dtype=np.float32)
    rmsd_all = np.empty(n_frames, dtype=np.float32)

    # Reference for RMSD: first frame, all CA
    u.trajectory[0]
    ref_pos = all_ca.positions.copy()

    # Process frame by frame
    for i, ts in enumerate(u.trajectory):
        if i % 200 == 0:
            print(f"  [{name}] frame {i}/{n_frames} ...", flush=True)

        # COM distance
        com_cgas = cgas.center_of_mass()
        com_trim41 = trim41.center_of_mass()
        com_dist[i] = np.linalg.norm(com_cgas - com_trim41)

        # CA-CA contacts (< 8A)
        # Use broadcasting for CA-CA distances, much smaller than all-atom
        d = np.sqrt(np.sum((cgas.positions[:, None, :] - trim41.positions[None, :, :]) ** 2, axis=2))
        contacts[i] = np.count_nonzero(d < CONTACT_CUTOFF)

        # Rg
        rg_cgas[i] = cgas_all.radius_of_gyration()
        rg_trim41[i] = trim41_all.radius_of_gyration()
        rg_total[i] = (cgas_all + trim41_all).radius_of_gyration()

        # RMSD vs first frame (no alignment — raw drift)
        rmsd_all[i] = rms.rmsd(all_ca.positions, ref_pos, superposition=False)

    # Compute per-residue RMSF (all CA)
    print(f"  [{name}] Computing RMSF ...")
    coords = np.empty((n_frames, n_all_ca, 3), dtype=np.float32)
    for i, ts in enumerate(u.trajectory):
        coords[i] = all_ca.positions
    rmsf = np.sqrt(np.mean((coords - coords.mean(axis=0)) ** 2, axis=(0, 2)))

    resids = all_ca.resids
    resnames = all_ca.resnames

    return {
        "com_dist": com_dist,
        "contacts": contacts,
        "rg_cgas": rg_cgas,
        "rg_trim41": rg_trim41,
        "rg_total": rg_total,
        "rmsd": rmsd_all,
        "rmsf": rmsf,
        "resids": resids,
        "resnames": resnames,
        "n_frames": n_frames,
        "cgas_last": cgas_last,
    }


def plot_comparison(results):
    """Generate comparison plots for all four systems."""
    systems = list(results.keys())
    colors = {"Hgal_WT": "#1f77b4", "Hgal_4mut_rev": "#ff7f0e",
              "Hsap_WT": "#2ca02c", "Hsap_4mut": "#d62728"}

    time_ns = {s: np.arange(results[s]["n_frames"]) * DT_NS for s in systems}

    # === Figure 1: Time series ===
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 2, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, 0])
    ax6 = fig.add_subplot(gs[2, 1])

    for s in systems:
        t = time_ns[s]
        ax1.plot(t, results[s]["com_dist"], label=s, color=colors[s], alpha=0.7, lw=0.5)
        ax2.plot(t, results[s]["contacts"], label=s, color=colors[s], alpha=0.7, lw=0.5)
        ax3.plot(t, results[s]["rg_total"], label=s, color=colors[s], alpha=0.7, lw=0.5)
        ax4.plot(t, results[s]["rg_cgas"], label=s, color=colors[s], alpha=0.7, lw=0.5)
        ax5.plot(t, results[s]["rg_trim41"], label=s, color=colors[s], alpha=0.7, lw=0.5)
        ax6.plot(t, results[s]["rmsd"], label=s, color=colors[s], alpha=0.7, lw=0.5)

    ax1.set_ylabel("COM distance (Å)")
    ax1.set_title("cGAS-TRIM41 COM Distance")
    ax1.legend(loc="upper right", fontsize=7)

    ax2.set_ylabel("CA-CA contacts (<8Å)")
    ax2.set_title("Interface Contacts (CA only)")
    ax2.legend(loc="upper right", fontsize=7)

    ax3.set_ylabel("Rg (Å)")
    ax3.set_title("Total Complex Rg")
    ax3.legend(loc="upper right", fontsize=7)

    ax4.set_ylabel("Rg (Å)")
    ax4.set_title("cGAS Domain Rg")
    ax4.legend(loc="upper right", fontsize=7)

    ax5.set_ylabel("Rg (Å)")
    ax5.set_title("TRIM41 Domain Rg")
    ax5.legend(loc="upper right", fontsize=7)

    ax6.set_ylabel("RMSD (Å)")
    ax6.set_title("All-CA RMSD vs Frame 0 (no alignment)")
    ax6.legend(loc="upper right", fontsize=7)
    ax6.set_xlabel("Time (ns)")

    fig.suptitle("Four-System MD Comparison (200 ns × 3 reps)", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(OUTDIR / "four_system_time_series.png", dpi=300)
    plt.close(fig)
    print(f"  Saved: {OUTDIR / 'four_system_time_series.png'}")

    # === Figure 2: RMSF comparison ===
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    axes = axes.flatten()

    for idx, pair in enumerate([("Hgal_WT", "Hgal_4mut_rev"), ("Hsap_WT", "Hsap_4mut")]):
        s1, s2 = pair
        ax = axes[idx]
        resids = results[s1]["resids"]
        cgas_last = results[s1]["cgas_last"]

        ax.plot(resids, results[s1]["rmsf"], label=s1, color=colors[s1], lw=1.2)
        ax.plot(resids, results[s2]["rmsf"], label=s2, color=colors[s2], lw=1.2)
        ax.axvline(cgas_last + 0.5, color="gray", ls="--", alpha=0.5)
        ax.text(cgas_last / 2, ax.get_ylim()[1] * 0.9, "cGAS", ha="center", fontsize=9)
        ax.text(cgas_last + (resids[-1] - cgas_last) / 2, ax.get_ylim()[1] * 0.9, "TRIM41", ha="center", fontsize=9)
        ax.set_xlabel("Residue Number")
        ax.set_ylabel("RMSF (Å)")
        ax.set_title(f"{s1} vs {s2}")
        ax.legend(loc="upper right")

    # Delta RMSF plots
    for idx, pair in enumerate([("Hgal_WT", "Hgal_4mut_rev"), ("Hsap_WT", "Hsap_4mut")], start=2):
        s1, s2 = pair
        ax = axes[idx]
        resids = results[s1]["resids"]
        cgas_last = results[s1]["cgas_last"]
        delta = results[s2]["rmsf"] - results[s1]["rmsf"]

        ax.plot(resids, delta, color="black", lw=1.0)
        ax.axhline(0, color="gray", ls="--", alpha=0.5)
        ax.axvline(cgas_last + 0.5, color="gray", ls="--", alpha=0.5)
        ax.fill_between(resids, delta, 0, where=(delta > 0), color="red", alpha=0.3, label=f"{s2} more flexible")
        ax.fill_between(resids, delta, 0, where=(delta < 0), color="blue", alpha=0.3, label=f"{s2} more rigid")
        ax.set_xlabel("Residue Number")
        ax.set_ylabel("ΔRMSF (Å)")
        ax.set_title(f"ΔRMSF: {s2} − {s1}")
        ax.legend(loc="upper right")

    fig.suptitle("Per-Residue RMSF Comparison", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(OUTDIR / "four_system_rmsf.png", dpi=300)
    plt.close(fig)
    print(f"  Saved: {OUTDIR / 'four_system_rmsf.png'}")

    # === Figure 3: Distributions ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    metrics = [("com_dist", "COM Distance (Å)"), ("contacts", "CA-CA Contacts"),
               ("rg_total", "Total Rg (Å)"), ("rmsd", "RMSD (Å)")]

    for ax, (metric, label) in zip(axes.flatten(), metrics):
        for s in systems:
            data = results[s][metric]
            ax.hist(data, bins=50, alpha=0.5, color=colors[s], label=s, density=True)
            ax.axvline(np.mean(data), color=colors[s], ls="--", lw=1.5)
        ax.set_xlabel(label)
        ax.set_ylabel("Density")
        ax.legend(loc="upper right", fontsize=7)

    fig.suptitle("Distribution Comparison", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(OUTDIR / "four_system_distributions.png", dpi=300)
    plt.close(fig)
    print(f"  Saved: {OUTDIR / 'four_system_distributions.png'}")


def summarize(results):
    """Print summary statistics and save to JSON."""
    summary = {}
    print("\n" + "=" * 60)
    print("FOUR-SYSTEM COMPARISON SUMMARY")
    print("=" * 60)

    for s in results:
        d = results[s]
        summary[s] = {
            "n_frames": int(d["n_frames"]),
            "com_dist_mean": float(np.mean(d["com_dist"])),
            "com_dist_std": float(np.std(d["com_dist"])),
            "contacts_mean": float(np.mean(d["contacts"])),
            "contacts_std": float(np.std(d["contacts"])),
            "rg_total_mean": float(np.mean(d["rg_total"])),
            "rg_total_std": float(np.std(d["rg_total"])),
            "rg_cgas_mean": float(np.mean(d["rg_cgas"])),
            "rg_trim41_mean": float(np.mean(d["rg_trim41"])),
            "rmsd_mean": float(np.mean(d["rmsd"])),
            "rmsd_std": float(np.std(d["rmsd"])),
        }
        print(f"\n{s}:")
        print(f"  Frames: {d['n_frames']}")
        print(f"  COM distance: {summary[s]['com_dist_mean']:.2f} ± {summary[s]['com_dist_std']:.2f} Å")
        print(f"  CA contacts:  {summary[s]['contacts_mean']:.1f} ± {summary[s]['contacts_std']:.1f}")
        print(f"  Total Rg:     {summary[s]['rg_total_mean']:.2f} ± {summary[s]['rg_total_std']:.2f} Å")
        print(f"  cGAS Rg:      {summary[s]['rg_cgas_mean']:.2f} Å")
        print(f"  TRIM41 Rg:    {summary[s]['rg_trim41_mean']:.2f} Å")
        print(f"  RMSD:         {summary[s]['rmsd_mean']:.2f} ± {summary[s]['rmsd_std']:.2f} Å")

    with open(OUTDIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved: {OUTDIR / 'summary.json'}")

    # Cross-system deltas
    print("\n" + "-" * 60)
    print("CROSS-SYSTEM DELTAS (mut − WT)")
    print("-" * 60)
    for pair in [("Hgal_WT", "Hgal_4mut_rev"), ("Hsap_WT", "Hsap_4mut")]:
        wt, mut = pair
        print(f"\n{mut} vs {wt}:")
        for metric in ["com_dist", "contacts", "rg_total", "rmsd"]:
            d_mut = np.mean(results[mut][metric])
            d_wt = np.mean(results[wt][metric])
            print(f"  Δ{metric}: {d_mut - d_wt:+.3f}")

    # Save full data as NPZ
    npz_data = {}
    for s in results:
        for key in ["com_dist", "contacts", "rg_cgas", "rg_trim41", "rg_total", "rmsd", "rmsf", "resids"]:
            npz_data[f"{s}_{key}"] = results[s][key]
    np.savez(OUTDIR / "four_system_data.npz", **npz_data)
    print(f"\n  Saved: {OUTDIR / 'four_system_data.npz'}")


def main():
    print("=" * 60)
    print("FOUR-SYSTEM MD COMPARISON (Optimized)")
    print("=" * 60)

    results = {}
    for name, info in SYSTEMS.items():
        print(f"\n>>> Processing {name} ...")
        try:
            results[name] = compute_system(name, info)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    print("\n>>> Plotting comparisons ...")
    plot_comparison(results)

    print("\n>>> Generating summary ...")
    summarize(results)

    print("\n>>> All done!")


if __name__ == "__main__":
    main()
