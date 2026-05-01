#!/usr/bin/env python3
"""Compare GROMACS 2026 (native ff) partial results with OpenMM."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from pathlib import Path
import sys

def analyze_trajectory(u, label, max_frames=None):
    """Analyze a trajectory: RMSD, COM distance, Rg."""
    protein = u.select_atoms('protein')
    
    # Determine cGAS and TRIM41 by residue range
    # Both GROMACS (pdb2gmx) and OpenMM (AMBER) use same resid numbering:
    # resid 1-218 = cGAS (NLEU...CPRO), resid 219-541 = TRIM41 (NASP...CPHE)
    cgas = u.select_atoms('protein and resid 1-218')
    trim41 = u.select_atoms('protein and resid 219-541')
    
    print(f"\n[{label}]")
    print(f"  Total protein atoms: {len(protein)}")
    print(f"  cGAS: {len(cgas.residues)} residues, {len(cgas)} atoms")
    print(f"  TRIM41: {len(trim41.residues)} residues, {len(trim41)} atoms")
    print(f"  Total frames: {len(u.trajectory)}")
    
    # Reference frame for RMSD (first frame)
    u.trajectory[0]
    ref_pos = protein.positions.copy()
    
    # Align on backbone for RMSD calculation
    bb = u.select_atoms('protein and backbone')
    ref_bb = bb.positions.copy()
    
    times = []
    rmsds = []
    com_dists = []
    rgs = []
    
    n_frames = len(u.trajectory) if max_frames is None else min(max_frames, len(u.trajectory))
    
    # Time per frame (ps)
    dt = 100.0  # OpenMM writes every 100ps (50000 steps @ 2fs)
    if 'xtc' in str(u.trajectory.filename).lower():
        # Try to get dt from trajectory, fallback to 10ps (GROMACS default nstxout=5000 @ 2fs)
        dt = getattr(u.trajectory, 'dt', 10.0)
    
    for i, ts in enumerate(u.trajectory[:n_frames]):
        # Align to reference on backbone
        mobile_bb = bb.positions
        R, rmsd_align = align.rotation_matrix(mobile_bb, ref_bb)
        
        # Apply rotation to full protein
        aligned = np.dot(protein.positions - bb.center_of_geometry(), R.T) + bb.center_of_geometry()
        
        # RMSD of full protein after alignment
        rmsd = np.sqrt(np.mean(np.sum((aligned - ref_pos)**2, axis=1)))
        
        # COM distance (cGAS vs TRIM41)
        com_cgas = cgas.center_of_geometry()
        com_trim41 = trim41.center_of_geometry()
        com_dist = np.linalg.norm(com_cgas - com_trim41)
        
        # Radius of gyration of full protein
        rg = protein.radius_of_gyration()
        
        # Use frame index × dt for consistent timing (DCD sometimes has wrong timestamps)
        t_ps = i * dt
        times.append(t_ps / 1000)  # ns
        rmsds.append(rmsd)
        com_dists.append(com_dist)
        rgs.append(rg)
    
    times = np.array(times)
    rmsds = np.array(rmsds)
    com_dists = np.array(com_dists)
    rgs = np.array(rgs)
    
    print(f"  Time range: {times[0]:.2f} - {times[-1]:.2f} ns")
    print(f"  RMSD:  mean={rmsds.mean():.2f} Å, std={rmsds.std():.2f} Å, max={rmsds.max():.2f} Å")
    print(f"  COM:   mean={com_dists.mean():.2f} Å, std={com_dists.std():.2f} Å, range=[{com_dists.min():.1f}, {com_dists.max():.1f}] Å")
    print(f"  Rg:    mean={rgs.mean():.2f} Å, std={rgs.std():.2f} Å")
    
    return times, rmsds, com_dists, rgs


def main():
    project = Path('/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation')
    
    # GROMACS 2026 (native ff)
    print("="*60)
    print("Loading GROMACS 2026 (native amber19sb.ff)...")
    u_gmx = mda.Universe(
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/Hsap_WT_rep1_ionized.gro',
        project / 'data/md_runs_gmx2026/Hsap_WT/rep1/prod.xtc'
    )
    t_gmx, rmsd_gmx, com_gmx, rg_gmx = analyze_trajectory(u_gmx, 'GROMACS 2026 native')
    
    # OpenMM Hsap_WT rep1 (same time range for fair comparison)
    print("\n" + "="*60)
    print("Loading OpenMM Hsap_WT rep1...")
    u_omm = mda.Universe(
        project / 'data/md_runs/Hsap_WT/Hsap_WT.prmtop',
        project / 'data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd'
    )
    # GROMACS: 2417 frames @ 10ps = 24.17ns
    # OpenMM: 2000 frames @ 100ps = 200ns
    # For fair comparison, take first ~24ns of OpenMM = 240 frames
    # But also analyze full 200ns for context
    print(f"  OpenMM total frames: {len(u_omm.trajectory)} (200ns)")
    print("  Analyzing first 240 frames (~24ns) for fair comparison with GROMACS")
    t_omm, rmsd_omm, com_omm, rg_omm = analyze_trajectory(u_omm, 'OpenMM (first 24ns)', max_frames=240)
    
    # Also analyze full OpenMM for reference
    t_omm_full, rmsd_omm_full, com_omm_full, rg_omm_full = analyze_trajectory(u_omm, 'OpenMM (full 200ns)', max_frames=None)
    
    # Summary comparison
    print("\n" + "="*60)
    print("COMPARISON (first ~24 ns)")
    print("="*60)
    print(f"{'Metric':<20} {'GROMACS 2026':<18} {'OpenMM 24ns':<18} {'Ratio':<10}")
    print("-"*60)
    print(f"{'RMSD mean (Å)':<20} {rmsd_gmx.mean():<18.2f} {rmsd_omm.mean():<18.2f} {rmsd_gmx.mean()/rmsd_omm.mean():<10.2f}")
    print(f"{'RMSD max (Å)':<20} {rmsd_gmx.max():<18.2f} {rmsd_omm.max():<18.2f} {rmsd_gmx.max()/rmsd_omm.max():<10.2f}")
    print(f"{'COM mean (Å)':<20} {com_gmx.mean():<18.2f} {com_omm.mean():<18.2f} {com_gmx.mean()/com_omm.mean():<10.2f}")
    print(f"{'COM std (Å)':<20} {com_gmx.std():<18.2f} {com_omm.std():<18.2f} {'-':<10}")
    print(f"{'Rg mean (Å)':<20} {rg_gmx.mean():<18.2f} {rg_omm.mean():<18.2f} {rg_gmx.mean()/rg_omm.mean():<10.2f}")
    
    print("\n" + "="*60)
    print("OpenMM full 200ns (for context)")
    print("="*60)
    print(f"{'RMSD mean (Å)':<20} {rmsd_omm_full.mean():<18.2f}")
    print(f"{'RMSD max (Å)':<20} {rmsd_omm_full.max():<18.2f}")
    print(f"{'COM mean (Å)':<20} {com_omm_full.mean():<18.2f}")
    print(f"{'Rg mean (Å)':<20} {rg_omm_full.mean():<18.2f}")
    
    # Verdict
    print("\n" + "="*60)
    print("VERDICT")
    print("="*60)
    rmsd_ratio = rmsd_gmx.mean() / rmsd_omm.mean()
    if rmsd_ratio < 1.5:
        print(f"✅ CMAP FIX LIKELY SUCCESSFUL: RMSD ratio = {rmsd_ratio:.2f} (< 1.5)")
        print("   GROMACS 2026 native ff produces similar RMSD to OpenMM.")
    elif rmsd_ratio < 2.5:
        print(f"🟡 PARTIAL IMPROVEMENT: RMSD ratio = {rmsd_ratio:.2f} (1.5-2.5)")
        print("   Better than old GROMACS (~4×), but still not fully consistent.")
        print("   May need more simulation time or additional parameter tuning.")
    else:
        print(f"🔴 STILL DIVERGENT: RMSD ratio = {rmsd_ratio:.2f} (> 2.5)")
        print("   RMSD still significantly higher than OpenMM.")
        print("   CMAP fix alone may not be sufficient.")
    
    # Write data files for plotting
    outdir = project / 'data/md_runs_gmx2026/Hsap_WT/rep1'
    np.savetxt(outdir / 'gmx2026_rmsd.dat', np.column_stack([t_gmx, rmsd_gmx, com_gmx, rg_gmx]),
               header='time_ns RMSD_A COM_A Rg_A', fmt='%.4f')
    np.savetxt(outdir / 'openmm_rmsd_24ns.dat', np.column_stack([t_omm, rmsd_omm, com_omm, rg_omm]),
               header='time_ns RMSD_A COM_A Rg_A', fmt='%.4f')
    np.savetxt(outdir / 'openmm_rmsd_200ns.dat', np.column_stack([t_omm_full, rmsd_omm_full, com_omm_full, rg_omm_full]),
               header='time_ns RMSD_A COM_A Rg_A', fmt='%.4f')
    print(f"\nData saved to {outdir}/gmx2026_rmsd.dat, openmm_rmsd_24ns.dat, openmm_rmsd_200ns.dat")


if __name__ == '__main__':
    main()
