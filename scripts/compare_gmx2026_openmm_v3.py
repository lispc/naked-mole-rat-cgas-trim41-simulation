#!/usr/bin/env python3
"""Compare GROMACS 2026 vs OpenMM with PBC unwrapping."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.transformations import unwrap
from pathlib import Path

def analyze_trajectory(u, label, dt_ps=10.0, max_frames=None):
    """Analyze with PBC unwrapping."""
    protein = u.select_atoms('protein')
    bb = u.select_atoms('protein and backbone')
    cgas = u.select_atoms('protein and resid 1-218')
    trim41 = u.select_atoms('protein and resid 219-541')
    
    # Using pre-unwrapped trajectory (gmx trjconv -pbc mol)
    
    # Reference = first frame (after unwrap)
    u.trajectory[0]
    ref_bb = bb.positions.copy()
    ref_prot = protein.positions.copy()
    
    times = []
    rmsds = []
    com_dists = []
    rgs = []
    
    n_frames = len(u.trajectory) if max_frames is None else min(max_frames, len(u.trajectory))
    
    for i, ts in enumerate(u.trajectory[:n_frames]):
        # Unwrap is applied automatically by the transformation
        
        # Kabsch align to reference backbone
        R, _ = align.rotation_matrix(bb.positions, ref_bb)
        mobile_cog = bb.center_of_geometry()
        ref_cog = ref_bb.mean(axis=0)
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
    
    # GROMACS 2026 prod (with PBC unwrap)
    print("="*60)
    print("GROMACS 2026 (native amber19sb.ff) — PBC unwrapped")
    u_gmx = mda.Universe(
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/Hsap_WT_rep1_ionized.gro',
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/prod_unwrapped.xtc'
    )
    t_gmx, rmsd_gmx, com_gmx, rg_gmx = analyze_trajectory(u_gmx, 'GROMACS', dt_ps=10.0)
    print(f"  Frames: {len(t_gmx)}, Time: {t_gmx[0]:.2f}-{t_gmx[-1]:.2f} ns")
    print(f"  Self-RMSD: mean={rmsd_gmx.mean():.2f} Å, max={rmsd_gmx.max():.2f} Å")
    print(f"  COM: mean={com_gmx.mean():.2f} Å, std={com_gmx.std():.2f} Å")
    print(f"  Rg: mean={rg_gmx.mean():.2f} Å, std={rg_gmx.std():.2f} Å")
    
    # OpenMM prod (first 25ns)
    print("\n" + "="*60)
    print("OpenMM Hsap_WT rep1 — first 25ns")
    u_omm = mda.Universe(
        project / 'data/md_runs/Hsap_WT/Hsap_WT.prmtop',
        project / 'data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd'
    )
    t_omm, rmsd_omm, com_omm, rg_omm = analyze_trajectory(u_omm, 'OpenMM', dt_ps=100.0, max_frames=250)
    print(f"  Frames: {len(t_omm)}, Time: {t_omm[0]:.2f}-{t_omm[-1]:.2f} ns")
    print(f"  Self-RMSD: mean={rmsd_omm.mean():.2f} Å, max={rmsd_omm.max():.2f} Å")
    print(f"  COM: mean={com_omm.mean():.2f} Å, std={com_omm.std():.2f} Å")
    print(f"  Rg: mean={rg_omm.mean():.2f} Å, std={rg_omm.std():.2f} Å")
    
    # Also compare GROMACS first frame to OpenMM reference
    print("\n" + "="*60)
    print("Absolute comparison (GROMACS vs OpenMM NPT reference)")
    u_ref = mda.Universe(
        project / 'data/md_runs/Hsap_WT/Hsap_WT.prmtop',
        project / 'data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_npt.dcd'
    )
    u_ref.trajectory[-1]
    ref_prot = u_ref.select_atoms('protein')
    ref_bb = ref_prot.select_atoms('backbone')
    ref_pos = ref_prot.positions.copy()
    
    u_gmx.trajectory[0]
    R, rmsd = align.rotation_matrix(u_gmx.select_atoms('protein and backbone').positions, ref_bb.positions)
    aligned = np.dot(u_gmx.select_atoms('protein').positions - u_gmx.select_atoms('protein and backbone').center_of_geometry(), R.T) + ref_bb.center_of_geometry()
    abs_rmsd = np.sqrt(np.mean(np.sum((aligned - ref_pos)**2, axis=1)))
    
    print(f"  GROMACS prod frame 0 vs OpenMM NPT: {abs_rmsd:.2f} Å")
    
    # Summary
    print("\n" + "="*60)
    print("COMPARISON (first ~25 ns, PBC-unwrapped)")
    print("="*60)
    print(f"{'Metric':<25} {'GROMACS':<18} {'OpenMM':<18} {'Ratio':<10}")
    print("-"*65)
    print(f"{'Self-RMSD mean (Å)':<25} {rmsd_gmx.mean():<18.2f} {rmsd_omm.mean():<18.2f} {rmsd_gmx.mean()/rmsd_omm.mean():<10.2f}")
    print(f"{'Self-RMSD max (Å)':<25} {rmsd_gmx.max():<18.2f} {rmsd_omm.max():<18.2f} {rmsd_gmx.max()/rmsd_omm.max():<10.2f}")
    print(f"{'COM mean (Å)':<25} {com_gmx.mean():<18.2f} {com_omm.mean():<18.2f} {com_gmx.mean()/com_omm.mean():<10.2f}")
    print(f"{'Rg mean (Å)':<25} {rg_gmx.mean():<18.2f} {rg_omm.mean():<18.2f} {rg_gmx.mean()/rg_omm.mean():<10.2f}")
    
    ratio = rmsd_gmx.mean() / rmsd_omm.mean()
    print("\n" + "="*60)
    print("VERDICT")
    print("="*60)
    if ratio < 1.5:
        print(f"✅ CMAP FIX SUCCESSFUL: Self-RMSD ratio = {ratio:.2f} (< 1.5)")
        print("   GROMACS 2026 native ff produces consistent dynamics with OpenMM.")
    elif ratio < 2.5:
        print(f"🟡 PARTIAL IMPROVEMENT: Self-RMSD ratio = {ratio:.2f} (1.5-2.5)")
        print("   Better than old GROMACS, but some differences remain.")
    else:
        print(f"🔴 STILL DIVERGENT: Self-RMSD ratio = {ratio:.2f} (> 2.5)")
        print("   Even after PBC correction, dynamics differ significantly.")
        print("   May need additional fixes (LINCS, ensemble, or other parameters).")
    
    # Save
    outdir = project / 'data/md_runs_gmx2026/Hsap_WT/rep1'
    np.savetxt(outdir / 'gmx2026_unwrapped.dat', np.column_stack([t_gmx, rmsd_gmx, com_gmx, rg_gmx]),
               header='time_ns self_RMSD_A COM_A Rg_A', fmt='%.4f')
    np.savetxt(outdir / 'openmm_unwrapped_25ns.dat', np.column_stack([t_omm, rmsd_omm, com_omm, rg_omm]),
               header='time_ns self_RMSD_A COM_A Rg_A', fmt='%.4f')

if __name__ == '__main__':
    main()
