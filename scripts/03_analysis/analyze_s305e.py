#!/usr/bin/env python3
"""
S305E vs WT analysis (200ns, 3 reps each).
Metrics: COM distance, H-bonds, RMSD, Rg.
Outputs: data/analysis/s305e_vs_wt/
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.distances import distance_array

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/s305e_vs_wt"
OUTDIR.mkdir(parents=True, exist_ok=True)


def analyze_system(name, prmtop, dcd_paths, dt_ns=0.1):
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")

    results = {
        'name': name,
        'com': [],
        'hbonds': [],
        'rg_cgas': [],
        'rg_trim': [],
        'rmsd': [],
        'time_ns': None,
    }

    for rep_idx, dcd in enumerate(dcd_paths):
        print(f"\n  Replica {rep_idx+1}: {dcd}")
        u = mda.Universe(str(prmtop), str(dcd))
        n_frames = len(u.trajectory)
        time_ns = np.arange(n_frames) * dt_ns
        if results['time_ns'] is None:
            results['time_ns'] = time_ns

        cgas = u.select_atoms("resid 219-541 and name CA")
        trim = u.select_atoms("resid 1-218 and name CA")
        protein_ca = u.select_atoms("protein and name CA")

        # H-bond analysis between cGAS and TRIM41
        hb = HydrogenBondAnalysis(
            universe=u,
            between=["resid 219-541", "resid 1-218"],
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150,
        )

        ref = mda.Universe(str(prmtop), str(dcd))
        ref_protein_ca = ref.select_atoms("protein and name CA")

        com_rep = []
        hb_rep = []
        rg_cgas_rep = []
        rg_trim_rep = []
        rmsd_rep = []

        # H-bond analysis (frame by frame to avoid memory issues)
        print("    Computing H-bonds...")
        hb.run(verbose=False)
        n_hb_per_frame = np.zeros(n_frames, dtype=int)
        for frame_idx, _ in hb.results.hbonds:
            n_hb_per_frame[int(frame_idx)] += 1
        # H-bonds already counted both directions by between parameter

        print("    Computing COM/Rg/RMSD...")
        for ts in u.trajectory:
            com_cgas = cgas.center_of_mass()
            com_trim = trim.center_of_mass()
            com_dist = np.linalg.norm(com_cgas - com_trim)
            com_rep.append(com_dist)

            rg_cgas_rep.append(cgas.radius_of_gyration())
            rg_trim_rep.append(trim.radius_of_gyration())

            r = rms.rmsd(protein_ca.positions, ref_protein_ca.positions, superposition=False)
            rmsd_rep.append(r)

            hb_rep.append(n_hb_per_frame[ts.frame])

        results['com'].append(np.array(com_rep))
        results['hbonds'].append(np.array(hb_rep))
        results['rg_cgas'].append(np.array(rg_cgas_rep))
        results['rg_trim'].append(np.array(rg_trim_rep))
        results['rmsd'].append(np.array(rmsd_rep))

        print(f"    COM: {np.mean(com_rep):.1f} ± {np.std(com_rep):.1f} Å")
        print(f"    H-bonds: {np.mean(hb_rep):.1f} ± {np.std(hb_rep):.1f}")
        print(f"    Rg cGAS: {np.mean(rg_cgas_rep):.1f} ± {np.std(rg_cgas_rep):.1f} Å")
        print(f"    RMSD: {np.mean(rmsd_rep):.1f} ± {np.std(rmsd_rep):.1f} Å")

    return results


def plot_comparison(wt, mutant, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    time = wt['time_ns']
    colors_wt = ['#1f77b4', '#4a90d9', '#87ceeb']
    colors_mut = ['#d62728', '#ff7f0e', '#ff9999']

    # 1. COM distance
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for i in range(3):
        axes[0].plot(time, wt['com'][i], color=colors_wt[i], alpha=0.7, label=f'WT rep{i+1}')
    axes[0].set_ylabel('COM Distance (Å)')
    axes[0].set_title('WT: cGAS-TRIM41 COM Distance')
    axes[0].legend(loc='upper right', fontsize=8)
    axes[0].set_ylim(30, 120)

    for i in range(3):
        axes[1].plot(time, mutant['com'][i], color=colors_mut[i], alpha=0.7, label=f'S305E rep{i+1}')
    axes[1].set_ylabel('COM Distance (Å)')
    axes[1].set_xlabel('Time (ns)')
    axes[1].set_title('S305E: cGAS-TRIM41 COM Distance')
    axes[1].legend(loc='upper right', fontsize=8)
    axes[1].set_ylim(30, 120)

    plt.tight_layout()
    fig.savefig(outdir / 'com_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'com_comparison.png'}")

    # 2. H-bonds
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for i in range(3):
        axes[0].plot(time, wt['hbonds'][i], color=colors_wt[i], alpha=0.7, label=f'WT rep{i+1}')
    axes[0].set_ylabel('H-bond Count')
    axes[0].set_title('WT: Interface H-bonds')
    axes[0].legend(loc='upper right', fontsize=8)

    for i in range(3):
        axes[1].plot(time, mutant['hbonds'][i], color=colors_mut[i], alpha=0.7, label=f'S305E rep{i+1}')
    axes[1].set_ylabel('H-bond Count')
    axes[1].set_xlabel('Time (ns)')
    axes[1].set_title('S305E: Interface H-bonds')
    axes[1].legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    fig.savefig(outdir / 'hbonds_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'hbonds_comparison.png'}")

    # 3. Rg cGAS
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(3):
        ax.plot(time, wt['rg_cgas'][i], color=colors_wt[i], alpha=0.5, linewidth=0.8)
        ax.plot(time, mutant['rg_cgas'][i], color=colors_mut[i], alpha=0.5, linewidth=0.8)

    wt_mean = np.mean([wt['rg_cgas'][i] for i in range(3)], axis=0)
    mut_mean = np.mean([mutant['rg_cgas'][i] for i in range(3)], axis=0)
    ax.plot(time, wt_mean, color='#1f77b4', linewidth=2, label='WT mean')
    ax.plot(time, mut_mean, color='#d62728', linewidth=2, label='S305E mean')

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('cGAS Rg (Å)')
    ax.set_title('cGAS Radius of Gyration')
    ax.legend()
    plt.tight_layout()
    fig.savefig(outdir / 'rg_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'rg_comparison.png'}")

    # 4. RMSD
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(3):
        ax.plot(time, wt['rmsd'][i], color=colors_wt[i], alpha=0.5, linewidth=0.8)
        ax.plot(time, mutant['rmsd'][i], color=colors_mut[i], alpha=0.5, linewidth=0.8)

    wt_mean_rmsd = np.mean([wt['rmsd'][i] for i in range(3)], axis=0)
    mut_mean_rmsd = np.mean([mutant['rmsd'][i] for i in range(3)], axis=0)
    ax.plot(time, wt_mean_rmsd, color='#1f77b4', linewidth=2, label='WT mean')
    ax.plot(time, mut_mean_rmsd, color='#d62728', linewidth=2, label='S305E mean')

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('Protein CA RMSD')
    ax.legend()
    plt.tight_layout()
    fig.savefig(outdir / 'rmsd_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'rmsd_comparison.png'}")

    # 5. Distributions
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    wt_com_all = np.concatenate(wt['com'])
    mut_com_all = np.concatenate(mutant['com'])
    axes[0, 0].hist(wt_com_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_com_all):.1f})', density=True)
    axes[0, 0].hist(mut_com_all, bins=50, alpha=0.6, color='#d62728', label=f'S305E (μ={np.mean(mut_com_all):.1f})', density=True)
    axes[0, 0].set_xlabel('COM Distance (Å)')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].set_title('COM Distance Distribution')
    axes[0, 0].legend()

    wt_hb_all = np.concatenate(wt['hbonds'])
    mut_hb_all = np.concatenate(mutant['hbonds'])
    axes[0, 1].hist(wt_hb_all, bins=range(0, max(wt_hb_all.max(), mut_hb_all.max())+2), alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_hb_all):.1f})', density=True)
    axes[0, 1].hist(mut_hb_all, bins=range(0, max(wt_hb_all.max(), mut_hb_all.max())+2), alpha=0.6, color='#d62728', label=f'S305E (μ={np.mean(mut_hb_all):.1f})', density=True)
    axes[0, 1].set_xlabel('H-bond Count')
    axes[0, 1].set_ylabel('Density')
    axes[0, 1].set_title('H-bond Distribution')
    axes[0, 1].legend()

    wt_rg_all = np.concatenate(wt['rg_cgas'])
    mut_rg_all = np.concatenate(mutant['rg_cgas'])
    axes[1, 0].hist(wt_rg_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_rg_all):.1f})', density=True)
    axes[1, 0].hist(mut_rg_all, bins=50, alpha=0.6, color='#d62728', label=f'S305E (μ={np.mean(mut_rg_all):.1f})', density=True)
    axes[1, 0].set_xlabel('cGAS Rg (Å)')
    axes[1, 0].set_ylabel('Density')
    axes[1, 0].set_title('cGAS Rg Distribution')
    axes[1, 0].legend()

    wt_rmsd_all = np.concatenate(wt['rmsd'])
    mut_rmsd_all = np.concatenate(mutant['rmsd'])
    axes[1, 1].hist(wt_rmsd_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_rmsd_all):.1f})', density=True)
    axes[1, 1].hist(mut_rmsd_all, bins=50, alpha=0.6, color='#d62728', label=f'S305E (μ={np.mean(mut_rmsd_all):.1f})', density=True)
    axes[1, 1].set_xlabel('RMSD (Å)')
    axes[1, 1].set_ylabel('Density')
    axes[1, 1].set_title('RMSD Distribution')
    axes[1, 1].legend()

    plt.tight_layout()
    fig.savefig(outdir / 'distributions_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'distributions_comparison.png'}")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS (all 200ns)")
    print(f"{'='*60}")
    print(f"{'Metric':<25} {'WT (3 reps)':<25} {'S305E (3 reps)':<25}")
    print("-" * 75)

    def fmt(data_list):
        arr = np.concatenate(data_list)
        return f"{np.mean(arr):.1f} ± {np.std(arr):.1f}"

    print(f"{'COM (Å)':<25} {fmt(wt['com']):<25} {fmt(mutant['com']):<25}")
    print(f"{'H-bonds':<25} {fmt(wt['hbonds']):<25} {fmt(mutant['hbonds']):<25}")
    print(f"{'Rg cGAS (Å)':<25} {fmt(wt['rg_cgas']):<25} {fmt(mutant['rg_cgas']):<25}")
    print(f"{'Rg TRIM41 (Å)':<25} {fmt(wt['rg_trim']):<25} {fmt(mutant['rg_trim']):<25}")
    print(f"{'RMSD (Å)':<25} {fmt(wt['rmsd']):<25} {fmt(mutant['rmsd']):<25}")

    # Final 50ns
    print(f"\n{'='*60}")
    print("FINAL 50ns AVERAGE (150-200ns)")
    print(f"{'='*60}")
    idx_150 = int(150 / 0.1)

    def fmt_final(data_list):
        arr = np.concatenate([d[idx_150:] for d in data_list])
        return f"{np.mean(arr):.1f} ± {np.std(arr):.1f}"

    print(f"{'COM (Å)':<25} {fmt_final(wt['com']):<25} {fmt_final(mutant['com']):<25}")
    print(f"{'H-bonds':<25} {fmt_final(wt['hbonds']):<25} {fmt_final(mutant['hbonds']):<25}")


def main():
    outdir = BASE / 'data/analysis/s305e_vs_wt'
    outdir.mkdir(parents=True, exist_ok=True)

    wt_dcds = [BASE / f'data/md_runs/Hsap_WT/rep{i}/Hsap_WT_rep{i}_prod.dcd' for i in range(1, 4)]
    mut_dcds = [BASE / f'data/md_runs/Hsap_WT_S305E/rep{i}/Hsap_WT_S305E_rep{i}_prod.dcd' for i in range(1, 4)]

    wt = analyze_system('WT', BASE / 'data/md_runs/Hsap_WT/Hsap_WT.prmtop', wt_dcds)
    mutant = analyze_system('S305E', BASE / 'data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E.prmtop', mut_dcds)

    plot_comparison(wt, mutant, outdir)

    print(f"\n{'='*60}")
    print(f"Analysis complete. Outputs: {outdir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
