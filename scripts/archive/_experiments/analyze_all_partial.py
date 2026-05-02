#!/usr/bin/env python3
"""Analyze all partial results: GROMACS 2026 + S305-phos."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def analyze_system(top, traj, label, dt_ps, max_frames=None, unwrap=False):
    u = mda.Universe(top, traj)
    protein = u.select_atoms('protein')
    bb = u.select_atoms('protein and backbone')
    cgas = u.select_atoms('protein and resid 1-218')
    trim41 = u.select_atoms('protein and resid 219-541')
    
    n_frames = len(u.trajectory) if max_frames is None else min(max_frames, len(u.trajectory))
    
    # Reference = first frame
    u.trajectory[0]
    ref_bb = bb.positions.copy()
    ref_prot = protein.positions.copy()
    
    times = []
    rmsds = []
    coms = []
    rgs = []
    
    for i, ts in enumerate(u.trajectory[:n_frames]):
        R, _ = align.rotation_matrix(bb.positions, ref_bb)
        mobile_cog = bb.center_of_geometry()
        ref_cog = ref_bb.mean(axis=0)
        aligned = np.dot(protein.positions - mobile_cog, R.T) + ref_cog
        rmsd = np.sqrt(np.mean(np.sum((aligned - ref_prot)**2, axis=1)))
        
        times.append(i * dt_ps / 1000)
        rmsds.append(rmsd)
        coms.append(np.linalg.norm(cgas.center_of_geometry() - trim41.center_of_geometry()))
        rgs.append(protein.radius_of_gyration())
    
    return (np.array(times), np.array(rmsds), np.array(coms), np.array(rgs),
            label, len(u.trajectory))


def plot_comparison(datasets, outpath):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    colors = {'GROMACS 2026': '#1f77b4', 'OpenMM WT': '#ff7f0e',
              'S305-phos rep1': '#2ca02c', 'S305-phos rep2': '#d62728', 'S305-phos rep3': '#9467bd'}
    
    for t, rmsd, com, rg, label, _ in datasets:
        c = colors.get(label, '#333333')
        axes[0,0].plot(t, rmsd, label=label, color=c, alpha=0.8)
        axes[0,1].plot(t, com, label=label, color=c, alpha=0.8)
        axes[1,0].plot(t, rg, label=label, color=c, alpha=0.8)
    
    axes[0,0].set_xlabel('Time (ns)')
    axes[0,0].set_ylabel('Self-RMSD (Å)')
    axes[0,0].set_title('Backbone RMSD')
    axes[0,0].legend(fontsize=8)
    axes[0,0].grid(True, alpha=0.3)
    
    axes[0,1].set_xlabel('Time (ns)')
    axes[0,1].set_ylabel('COM distance (Å)')
    axes[0,1].set_title('cGAS-TRIM41 COM Distance')
    axes[0,1].legend(fontsize=8)
    axes[0,1].grid(True, alpha=0.3)
    
    axes[1,0].set_xlabel('Time (ns)')
    axes[1,0].set_ylabel('Rg (Å)')
    axes[1,0].set_title('Radius of Gyration')
    axes[1,0].legend(fontsize=8)
    axes[1,0].grid(True, alpha=0.3)
    
    # Summary table
    axes[1,1].axis('off')
    rows = []
    for t, rmsd, com, rg, label, total in datasets:
        rows.append([label, f"{t[-1]:.1f}ns", f"{rmsd.mean():.1f}±{rmsd.std():.1f}",
                     f"{com.mean():.1f}±{com.std():.1f}", f"{rg.mean():.1f}±{rg.std():.1f}"])
    
    table = axes[1,1].table(cellText=rows,
                            colLabels=['System', 'Length', 'RMSD (Å)', 'COM (Å)', 'Rg (Å)'],
                            loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)
    axes[1,1].set_title('Summary Statistics', pad=20)
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    print(f"Plot saved to {outpath}")


def main():
    project = Path('/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation')
    datasets = []
    
    # 1. GROMACS 2026 (77ns)
    print("Loading GROMACS 2026...")
    d = analyze_system(
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/Hsap_WT_rep1_ionized.gro',
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/prod_unwrapped_latest.xtc',
        'GROMACS 2026', dt_ps=10.0)
    datasets.append(d)
    print(f"  {d[4]}: {d[5]} frames, {d[0][-1]:.1f}ns, RMSD={d[1].mean():.1f}Å")
    
    # 2. OpenMM WT (first 80ns for fair comparison)
    print("Loading OpenMM WT...")
    d = analyze_system(
        project / 'data/md_runs/Hsap_WT/Hsap_WT.prmtop',
        project / 'data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd',
        'OpenMM WT', dt_ps=100.0, max_frames=800)
    datasets.append(d)
    print(f"  {d[4]}: {d[5]} frames, {d[0][-1]:.1f}ns, RMSD={d[1].mean():.1f}Å")
    
    # 3-5. S305-phos 3 replicas
    for rep in [1, 2, 3]:
        print(f"Loading S305-phos rep{rep}...")
        dcd = project / f'data/md_runs/Hsap_WT_S305phos/rep{rep}/Hsap_WT_S305phos_rep{rep}_prod.dcd'
        if dcd.exists():
            d = analyze_system(
                project / 'data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop',
                dcd, f'S305-phos rep{rep}', dt_ps=100.0)
            datasets.append(d)
            print(f"  {d[4]}: {d[5]} frames, {d[0][-1]:.1f}ns, RMSD={d[1].mean():.1f}Å")
    
    # Plot
    out = project / 'figures/partial_results_comparison.png'
    out.parent.mkdir(parents=True, exist_ok=True)
    plot_comparison(datasets, out)
    
    # Print comparison
    print("\n" + "="*70)
    print("GROMACS 2026 vs OpenMM WT (first ~80ns)")
    print("="*70)
    gmx = datasets[0]
    omm = datasets[1]
    print(f"{'Metric':<20} {'GROMACS':<15} {'OpenMM':<15} {'Ratio':<10}")
    print("-"*60)
    print(f"{'RMSD mean (Å)':<20} {gmx[1].mean():<15.2f} {omm[1].mean():<15.2f} {gmx[1].mean()/omm[1].mean():<10.2f}")
    print(f"{'COM mean (Å)':<20} {gmx[2].mean():<15.2f} {omm[2].mean():<15.2f} {gmx[2].mean()/omm[2].mean():<10.2f}")
    print(f"{'Rg mean (Å)':<20} {gmx[3].mean():<15.2f} {omm[3].mean():<15.2f} {gmx[3].mean()/omm[3].mean():<10.2f}")
    
    print("\n" + "="*70)
    print("S305-phos vs WT (first ~80ns each)")
    print("="*70)
    wt = datasets[1]
    for d in datasets[2:]:
        print(f"\n{d[4]}:")
        print(f"  RMSD: {d[1].mean():.2f} Å (WT: {wt[1].mean():.2f} Å, ratio: {d[1].mean()/wt[1].mean():.2f})")
        print(f"  COM:  {d[2].mean():.2f} Å (WT: {wt[2].mean():.2f} Å)")
        print(f"  Rg:   {d[3].mean():.2f} Å (WT: {wt[3].mean():.2f} Å)")
    
    print(f"\nPlot saved to: {out}")

if __name__ == '__main__':
    main()
