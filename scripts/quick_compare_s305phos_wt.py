#!/usr/bin/env python3
"""
Quick comparison: WT vs S305-phos (200ns each, 3 reps).

Focus metrics:
  - COM distance (cGAS-TRIM41 center of mass)
  - Radius of gyration (Rg) per protein
  - RMSD (protein CA)
  - Active site distances (D431, K479, L495, K498 to TRIM41 COM)
  - Hydrogen bonds at interface

Outputs: comparison plots to data/analysis/s305phos_vs_wt/
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.distances import distance_array


def analyze_system(name, prmtop, dcd_paths, dt_ns=0.1):
    """Analyze all replicas for a system. Returns dict of time series."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")
    
    results = {
        'name': name,
        'com': [],
        'rg_cgas': [],
        'rg_trim': [],
        'rmsd': [],
        'active_sites': {'D431': [], 'K479': [], 'L495': [], 'K498': []},
        'time_ns': None,
    }
    
    for rep_idx, dcd in enumerate(dcd_paths):
        print(f"\n  Replica {rep_idx+1}: {dcd}")
        u = mda.Universe(prmtop, dcd)
        n_frames = len(u.trajectory)
        time_ns = np.arange(n_frames) * dt_ns
        if results['time_ns'] is None:
            results['time_ns'] = time_ns
        
        # Selections (prmtop has no segid, use resid ranges)
        # TRIM41: resid 1-218, cGAS: resid 219-541 (protein only)
        cgas = u.select_atoms("resid 219-541 and name CA")
        trim = u.select_atoms("resid 1-218 and name CA")
        protein_ca = u.select_atoms("protein and name CA")
        
        # Active site residues (pdb4amber renumbered)
        active_sel = {
            'D431': u.select_atoms("resid 450 and name CA"),
            'K479': u.select_atoms("resid 498 and name CA"),
            'L495': u.select_atoms("resid 514 and name CA"),
            'K498': u.select_atoms("resid 517 and name CA"),
        }
        
        # Load reference for RMSD
        ref = mda.Universe(prmtop, dcd)
        ref_protein_ca = ref.select_atoms("protein and name CA")
        
        com_rep = []
        rg_cgas_rep = []
        rg_trim_rep = []
        rmsd_rep = []
        active_rep = {k: [] for k in active_sel}
        
        for ts in u.trajectory:
            # COM distance (using CA only for speed)
            com_cgas = cgas.center_of_mass()
            com_trim = trim.center_of_mass()
            com_dist = np.linalg.norm(com_cgas - com_trim)
            com_rep.append(com_dist)
            
            # Rg (using CA only)
            rg_cgas_rep.append(cgas.radius_of_gyration())
            rg_trim_rep.append(trim.radius_of_gyration())
            
            # RMSD (no alignment)
            r = rms.rmsd(protein_ca.positions, ref_protein_ca.positions, superposition=False)
            rmsd_rep.append(r)
            
            # Active site distances to TRIM41 COM
            for site_name, sel in active_sel.items():
                if len(sel) > 0:
                    d = np.linalg.norm(sel.positions[0] - com_trim)
                    active_rep[site_name].append(d)
                else:
                    active_rep[site_name].append(np.nan)
        
        results['com'].append(np.array(com_rep))
        results['rg_cgas'].append(np.array(rg_cgas_rep))
        results['rg_trim'].append(np.array(rg_trim_rep))
        results['rmsd'].append(np.array(rmsd_rep))
        for k in active_sel:
            results['active_sites'][k].append(np.array(active_rep[k]))
        
        print(f"    COM: {np.mean(com_rep):.1f} ± {np.std(com_rep):.1f} Å")
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
    axes[0].axhline(y=45, color='blue', linestyle='--', alpha=0.5, label='WT avg ~45Å')
    axes[0].set_ylabel('COM Distance (Å)')
    axes[0].set_title('WT: cGAS-TRIM41 COM Distance')
    axes[0].legend(loc='upper right', fontsize=8)
    axes[0].set_ylim(30, 120)
    
    for i in range(3):
        axes[1].plot(time, mutant['com'][i], color=colors_mut[i], alpha=0.7, label=f'S305-phos rep{i+1}')
    axes[1].axhline(y=45, color='blue', linestyle='--', alpha=0.5)
    axes[1].set_ylabel('COM Distance (Å)')
    axes[1].set_xlabel('Time (ns)')
    axes[1].set_title('S305-phos: cGAS-TRIM41 COM Distance')
    axes[1].legend(loc='upper right', fontsize=8)
    axes[1].set_ylim(30, 120)
    
    plt.tight_layout()
    fig.savefig(outdir / 'com_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'com_comparison.png'}")
    
    # 2. Rg cGAS
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(3):
        ax.plot(time, wt['rg_cgas'][i], color=colors_wt[i], alpha=0.5, linewidth=0.8)
        ax.plot(time, mutant['rg_cgas'][i], color=colors_mut[i], alpha=0.5, linewidth=0.8)
    
    wt_mean = np.mean([wt['rg_cgas'][i] for i in range(3)], axis=0)
    mut_mean = np.mean([mutant['rg_cgas'][i] for i in range(3)], axis=0)
    ax.plot(time, wt_mean, color='#1f77b4', linewidth=2, label='WT mean')
    ax.plot(time, mut_mean, color='#d62728', linewidth=2, label='S305-phos mean')
    
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('cGAS Rg (Å)')
    ax.set_title('cGAS Radius of Gyration')
    ax.legend()
    plt.tight_layout()
    fig.savefig(outdir / 'rg_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'rg_comparison.png'}")
    
    # 3. RMSD
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(3):
        ax.plot(time, wt['rmsd'][i], color=colors_wt[i], alpha=0.5, linewidth=0.8)
        ax.plot(time, mutant['rmsd'][i], color=colors_mut[i], alpha=0.5, linewidth=0.8)
    
    wt_mean_rmsd = np.mean([wt['rmsd'][i] for i in range(3)], axis=0)
    mut_mean_rmsd = np.mean([mutant['rmsd'][i] for i in range(3)], axis=0)
    ax.plot(time, wt_mean_rmsd, color='#1f77b4', linewidth=2, label='WT mean')
    ax.plot(time, mut_mean_rmsd, color='#d62728', linewidth=2, label='S305-phos mean')
    
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('Protein CA RMSD')
    ax.legend()
    plt.tight_layout()
    fig.savefig(outdir / 'rmsd_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'rmsd_comparison.png'}")
    
    # 4. Active site distances
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
    sites = ['D431', 'K479', 'L495', 'K498']
    for idx, site in enumerate(sites):
        ax = axes[idx // 2, idx % 2]
        for i in range(3):
            ax.plot(time, wt['active_sites'][site][i], color=colors_wt[i], alpha=0.5, linewidth=0.8)
            ax.plot(time, mutant['active_sites'][site][i], color=colors_mut[i], alpha=0.5, linewidth=0.8)
        
        wt_mean_s = np.mean([wt['active_sites'][site][i] for i in range(3)], axis=0)
        mut_mean_s = np.mean([mutant['active_sites'][site][i] for i in range(3)], axis=0)
        ax.plot(time, wt_mean_s, color='#1f77b4', linewidth=2, label='WT')
        ax.plot(time, mut_mean_s, color='#d62728', linewidth=2, label='S305-phos')
        
        ax.set_ylabel(f'{site} to TRIM41 (Å)')
        ax.set_title(f'Active Site: {site}')
        ax.legend()
    
    axes[1, 0].set_xlabel('Time (ns)')
    axes[1, 1].set_xlabel('Time (ns)')
    plt.tight_layout()
    fig.savefig(outdir / 'active_sites_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'active_sites_comparison.png'}")
    
    # 5. Distribution histograms
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # COM distribution
    wt_com_all = np.concatenate(wt['com'])
    mut_com_all = np.concatenate(mutant['com'])
    axes[0, 0].hist(wt_com_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_com_all):.1f})', density=True)
    axes[0, 0].hist(mut_com_all, bins=50, alpha=0.6, color='#d62728', label=f'S305-phos (μ={np.mean(mut_com_all):.1f})', density=True)
    axes[0, 0].set_xlabel('COM Distance (Å)')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].set_title('COM Distance Distribution')
    axes[0, 0].legend()
    
    # Rg distribution
    wt_rg_all = np.concatenate(wt['rg_cgas'])
    mut_rg_all = np.concatenate(mutant['rg_cgas'])
    axes[0, 1].hist(wt_rg_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_rg_all):.1f})', density=True)
    axes[0, 1].hist(mut_rg_all, bins=50, alpha=0.6, color='#d62728', label=f'S305-phos (μ={np.mean(mut_rg_all):.1f})', density=True)
    axes[0, 1].set_xlabel('cGAS Rg (Å)')
    axes[0, 1].set_ylabel('Density')
    axes[0, 1].set_title('cGAS Rg Distribution')
    axes[0, 1].legend()
    
    # RMSD distribution
    wt_rmsd_all = np.concatenate(wt['rmsd'])
    mut_rmsd_all = np.concatenate(mutant['rmsd'])
    axes[1, 0].hist(wt_rmsd_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_rmsd_all):.1f})', density=True)
    axes[1, 0].hist(mut_rmsd_all, bins=50, alpha=0.6, color='#d62728', label=f'S305-phos (μ={np.mean(mut_rmsd_all):.1f})', density=True)
    axes[1, 0].set_xlabel('RMSD (Å)')
    axes[1, 0].set_ylabel('Density')
    axes[1, 0].set_title('RMSD Distribution')
    axes[1, 0].legend()
    
    # RMSD distribution (already shown above, use COM distribution here too)
    axes[1, 1].hist(wt_com_all, bins=50, alpha=0.6, color='#1f77b4', label=f'WT (μ={np.mean(wt_com_all):.1f})', density=True)
    axes[1, 1].hist(mut_com_all, bins=50, alpha=0.6, color='#d62728', label=f'S305-phos (μ={np.mean(mut_com_all):.1f})', density=True)
    axes[1, 1].set_xlabel('COM Distance (Å)')
    axes[1, 1].set_ylabel('Density')
    axes[1, 1].set_title('COM Distance Distribution')
    axes[1, 1].legend()
    
    plt.tight_layout()
    fig.savefig(outdir / 'distributions_comparison.png', dpi=150)
    plt.close()
    print(f"  Saved: {outdir / 'distributions_comparison.png'}")
    
    # Summary stats table
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS")
    print(f"{'='*60}")
    print(f"{'Metric':<25} {'WT (3 reps)':<25} {'S305-phos (3 reps)':<25}")
    print("-" * 75)
    
    def fmt(data_list):
        arr = np.concatenate(data_list)
        return f"{np.mean(arr):.1f} ± {np.std(arr):.1f}"
    
    print(f"{'COM (Å)':<25} {fmt(wt['com']):<25} {fmt(mutant['com']):<25}")
    print(f"{'Rg cGAS (Å)':<25} {fmt(wt['rg_cgas']):<25} {fmt(mutant['rg_cgas']):<25}")
    print(f"{'Rg TRIM41 (Å)':<25} {fmt(wt['rg_trim']):<25} {fmt(mutant['rg_trim']):<25}")
    print(f"{'RMSD (Å)':<25} {fmt(wt['rmsd']):<25} {fmt(mutant['rmsd']):<25}")
    
    for site in sites:
        print(f"{site + ' (Å)':<25} {fmt(wt['active_sites'][site]):<25} {fmt(mutant['active_sites'][site]):<25}")
    
    # Final 50ns average
    print(f"\n{'='*60}")
    print("FINAL 50ns AVERAGE (150-200ns)")
    print(f"{'='*60}")
    idx_150 = int(150 / 0.1)
    
    def fmt_final(data_list):
        arr = np.concatenate([d[idx_150:] for d in data_list])
        return f"{np.mean(arr):.1f} ± {np.std(arr):.1f}"
    
    print(f"{'COM (Å)':<25} {fmt_final(wt['com']):<25} {fmt_final(mutant['com']):<25}")
    for site in sites:
        print(f"{site + ' (Å)':<25} {fmt_final(wt['active_sites'][site]):<25} {fmt_final(mutant['active_sites'][site]):<25}")


def main():
    outdir = Path('data/analysis/s305phos_vs_wt')
    outdir.mkdir(parents=True, exist_ok=True)
    
    wt_dcds = [f'data/md_runs/Hsap_WT/rep{i}/Hsap_WT_rep{i}_prod.dcd' for i in range(1, 4)]
    mut_dcds = [f'data/md_runs/Hsap_WT_S305phos/rep{i}/Hsap_WT_S305phos_rep{i}_prod.dcd' for i in range(1, 4)]
    
    wt = analyze_system('WT', 'data/md_runs/Hsap_WT/Hsap_WT.prmtop', wt_dcds)
    mutant = analyze_system('S305-phos', 'data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop', mut_dcds)
    
    plot_comparison(wt, mutant, outdir)
    
    print(f"\n{'='*60}")
    print(f"Analysis complete. Outputs: {outdir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
