#!/usr/bin/env python3
"""Compare GROMACS 2026 vs OpenMM using self-referenced RMSD (fair comparison)."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from pathlib import Path

def analyze_self_rmsd(u, label, dt_ps=10.0):
    """Compute RMSD referenced to the trajectory's own first frame."""
    protein = u.select_atoms('protein')
    bb = u.select_atoms('protein and backbone')
    cgas = u.select_atoms('protein and resid 1-218')
    trim41 = u.select_atoms('protein and resid 219-541')
    
    # Reference = first frame
    u.trajectory[0]
    ref_bb = bb.positions.copy()
    ref_prot = protein.positions.copy()
    
    times = []
    rmsds = []
    com_dists = []
    rgs = []
    
    for i, ts in enumerate(u.trajectory):
        # Kabsch align to reference backbone
        R, _ = align.rotation_matrix(bb.positions, ref_bb)
        aligned = np.dot(protein.positions - bb.center_of_geometry(), R.T) + ref_bb.mean(axis=0)
        # Actually align.alignto modifies in place, let's use it properly
        
    # Re-do with proper alignment
    u.trajectory[0]
    ref_bb = bb.positions.copy()
    ref_cog = bb.center_of_geometry()
    
    for i, ts in enumerate(u.trajectory):
        mobile_bb = bb.positions
        R, _ = align.rotation_matrix(mobile_bb, ref_bb)
        # Center mobile, rotate, then add ref center
        mobile_cog = bb.center_of_geometry()
        aligned = np.dot(protein.positions - mobile_cog, R.T) + ref_cog
        rmsd = np.sqrt(np.mean(np.sum((aligned - ref_prot)**2, axis=1)))
        
        com_dist = np.linalg.norm(cgas.center_of_geometry() - trim41.center_of_geometry())
        rg = protein.radius_of_gyration()
        
        times.append(i * dt_ps / 1000)
        rmsds.append(rmsd)
        com_dists.append(com_dist)
        rgs.append(rg)
    
    return np.array(times), np.array(rmsds), np.array(com_dists), np.array(rgs)


def main():
    project = Path('/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation')
    
    # GROMACS 2026 (self-referenced: first prod frame)
    print("="*60)
    print("GROMACS 2026 (native amber19sb.ff)")
    u_gmx = mda.Universe(
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/Hsap_WT_rep1_ionized.gro',
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/prod.xtc'
    )
    # GROMACS xtc dt: check
    dt_gmx = getattr(u_gmx.trajectory, 'dt', 10.0)
    t_gmx, rmsd_gmx, com_gmx, rg_gmx = analyze_self_rmsd(u_gmx, 'GROMACS', dt_ps=dt_gmx)
    print(f"  Frames: {len(t_gmx)}, Time: {t_gmx[0]:.2f}-{t_gmx[-1]:.2f} ns")
    print(f"  Self-RMSD: mean={rmsd_gmx.mean():.2f} Å, max={rmsd_gmx.max():.2f} Å")
    print(f"  COM: mean={com_gmx.mean():.2f} Å, std={com_gmx.std():.2f} Å")
    print(f"  Rg: mean={rg_gmx.mean():.2f} Å, std={rg_gmx.std():.2f} Å")
    
    # OpenMM (self-referenced: first prod frame, first 240 frames = 24ns)
    print("\n" + "="*60)
    print("OpenMM Hsap_WT rep1 (first 24ns)")
    u_omm = mda.Universe(
        project / 'data/md_runs/Hsap_WT/Hsap_WT.prmtop',
        project / 'data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd'
    )
    t_omm, rmsd_omm, com_omm, rg_omm = analyze_self_rmsd(u_omm, 'OpenMM', dt_ps=100.0)
    # Only first 24ns
    mask = t_omm <= 24.5
    t_omm = t_omm[mask]
    rmsd_omm = rmsd_omm[mask]
    com_omm = com_omm[mask]
    rg_omm = rg_omm[mask]
    print(f"  Frames: {len(t_omm)}, Time: {t_omm[0]:.2f}-{t_omm[-1]:.2f} ns")
    print(f"  Self-RMSD: mean={rmsd_omm.mean():.2f} Å, max={rmsd_omm.max():.2f} Å")
    print(f"  COM: mean={com_omm.mean():.2f} Å, std={com_omm.std():.2f} Å")
    print(f"  Rg: mean={rg_omm.mean():.2f} Å, std={rg_omm.std():.2f} Å")
    
    # Also full OpenMM for context
    print("\n" + "="*60)
    print("OpenMM Hsap_WT rep1 (full 200ns)")
    t_ommf, rmsd_ommf, com_ommf, rg_ommf = analyze_self_rmsd(u_omm, 'OpenMM', dt_ps=100.0)
    print(f"  Self-RMSD: mean={rmsd_ommf.mean():.2f} Å, max={rmsd_ommf.max():.2f} Å")
    print(f"  COM: mean={com_ommf.mean():.2f} Å, std={com_ommf.std():.2f} Å")
    
    # Comparison
    print("\n" + "="*60)
    print("COMPARISON (self-referenced RMSD)")
    print("="*60)
    print(f"{'Metric':<25} {'GROMACS 24ns':<18} {'OpenMM 24ns':<18} {'Ratio':<10}")
    print("-"*65)
    print(f"{'Self-RMSD mean (Å)':<25} {rmsd_gmx.mean():<18.2f} {rmsd_omm.mean():<18.2f} {rmsd_gmx.mean()/rmsd_omm.mean():<10.2f}")
    print(f"{'Self-RMSD max (Å)':<25} {rmsd_gmx.max():<18.2f} {rmsd_omm.max():<18.2f} {rmsd_gmx.max()/rmsd_omm.max():<10.2f}")
    print(f"{'COM mean (Å)':<25} {com_gmx.mean():<18.2f} {com_omm.mean():<18.2f} {com_gmx.mean()/com_omm.mean():<10.2f}")
    print(f"{'COM std (Å)':<25} {com_gmx.std():<18.2f} {com_omm.std():<18.2f} {'-':<10}")
    print(f"{'Rg mean (Å)':<25} {rg_gmx.mean():<18.2f} {rg_omm.mean():<18.2f} {rg_gmx.mean()/rg_omm.mean():<10.2f}")
    
    # Verdict
    ratio = rmsd_gmx.mean() / rmsd_omm.mean()
    print("\n" + "="*60)
    print("VERDICT (self-referenced)")
    print("="*60)
    if ratio < 1.5:
        print(f"✅ CMAP FIX SUCCESSFUL: Self-RMSD ratio = {ratio:.2f} (< 1.5)")
    elif ratio < 2.5:
        print(f"🟡 PARTIAL IMPROVEMENT: Self-RMSD ratio = {ratio:.2f} (1.5-2.5)")
    else:
        print(f"🔴 STILL DIVERGENT: Self-RMSD ratio = {ratio:.2f} (> 2.5)")
        print("   GROMACS 2026 prod dynamics still deviate from OpenMM.")
        print("   CMAP fix may not be sufficient; other factors at play.")
    
    # Save data
    outdir = project / 'data/md_runs_gmx2026/Hsap_WT/rep1'
    np.savetxt(outdir / 'gmx2026_self_rmsd.dat', np.column_stack([t_gmx, rmsd_gmx, com_gmx, rg_gmx]),
               header='time_ns self_RMSD_A COM_A Rg_A', fmt='%.4f')
    np.savetxt(outdir / 'openmm_self_rmsd_24ns.dat', np.column_stack([t_omm, rmsd_omm, com_omm, rg_omm]),
               header='time_ns self_RMSD_A COM_A Rg_A', fmt='%.4f')

if __name__ == '__main__':
    main()
