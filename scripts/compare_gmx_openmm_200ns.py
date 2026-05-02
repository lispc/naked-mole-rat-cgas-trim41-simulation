"""Compare GROMACS 2026 vs OpenMM 200ns trajectories for Hsap_WT."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

# Paths
GRO   = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs_gmx2026/Hsap_WT/rep1/prod.gro"
GMX_XTC = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs_gmx2026/Hsap_WT/rep1/prod.xtc"
OMM_TOP = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/Hsap_WT.prmtop"
OMM_DCD = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd"
OUTDIR  = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/analysis/gmx_openmm_comparison"
os.makedirs(OUTDIR, exist_ok=True)

print("Loading GROMACS trajectory...")
u_gmx = mda.Universe(GRO, GMX_XTC)
print(f"  GROMACS: {len(u_gmx.trajectory)} frames")

print("Loading OpenMM trajectory...")
u_omm = mda.Universe(OMM_TOP, OMM_DCD)
print(f"  OpenMM:  {len(u_omm.trajectory)} frames")

# Selections
# GROMACS: resid 1-218 = TRIM41, 219-541 = cGAS
# OpenMM:  resid 1-218 = TRIM41, 219-541 = cGAS
trim41_gmx = u_gmx.select_atoms("protein and resid 1-218")
cgas_gmx   = u_gmx.select_atoms("protein and resid 219-541")
prot_gmx   = u_gmx.select_atoms("protein")

trim41_omm = u_omm.select_atoms("protein and resid 1-218")
cgas_omm   = u_omm.select_atoms("protein and resid 219-541")
prot_omm   = u_omm.select_atoms("protein")

print(f"  TRIM41 atoms: GMX={len(trim41_gmx)}, OMM={len(trim41_omm)}")
print(f"  cGAS atoms:   GMX={len(cgas_gmx)}, OMM={len(cgas_omm)}")

# Reference structures for RMSD
ref_gmx = mda.Universe(GRO)
ref_omm = mda.Universe(OMM_TOP, OMM_DCD)

# Align reference to itself
align.alignto(ref_gmx, ref_gmx, select="protein and name CA", weights="mass")
align.alignto(ref_omm, ref_omm, select="protein and name CA", weights="mass")

# Analysis functions
def analyze(u, trim41, cgas, prot, ref, label, step=10):
    """Analyze trajectory. step=10 means sample every 10 frames for GROMACS (100 ps)."""
    n_frames = len(u.trajectory)
    indices = range(0, n_frames, step)
    n_sample = len(indices)
    
    times = np.zeros(n_sample)
    com_dist = np.zeros(n_sample)
    rg_trim = np.zeros(n_sample)
    rg_cgas = np.zeros(n_sample)
    rmsd_all = np.zeros(n_sample)
    rmsd_trim = np.zeros(n_sample)
    rmsd_cgas = np.zeros(n_sample)
    
    # Reference CA coords
    ref_ca = ref.select_atoms("protein and name CA")
    ref_ca_coords = ref_ca.positions.copy()
    ref_trim_ca = ref.select_atoms("protein and resid 1-218 and name CA")
    ref_cgas_ca = ref.select_atoms("protein and resid 219-541 and name CA")
    
    prev_com_trim = None
    box_prev = None
    
    for i, idx in enumerate(indices):
        ts = u.trajectory[idx]
        times[i] = ts.time / 1000.0  # ps -> ns
        
        # Align to reference CA
        mobile_ca = prot.select_atoms("name CA")
        R, rmsd_val = align.rotation_matrix(mobile_ca.positions, ref_ca_coords, weights=mobile_ca.masses)
        prot.atoms.translate(-prot.center_of_mass())
        prot.atoms.rotate(R)
        prot.atoms.translate(ref.select_atoms("protein").center_of_mass())
        
        # COM distance
        com_trim = trim41.center_of_mass()
        com_cgas = cgas.center_of_mass()
        
        # PBC correction for COM distance (manual)
        if prev_com_trim is not None and ts.dimensions is not None:
            box = ts.dimensions[:3]
            for dim in range(3):
                if box[dim] > 0:
                    delta = com_trim[dim] - prev_com_trim[dim]
                    if delta > box[dim] / 2:
                        com_trim[dim] -= box[dim]
                    elif delta < -box[dim] / 2:
                        com_trim[dim] += box[dim]
        
        dist = np.linalg.norm(com_trim - com_cgas)
        com_dist[i] = dist
        prev_com_trim = com_trim.copy()
        
        # Rg
        rg_trim[i] = trim41.radius_of_gyration()
        rg_cgas[i] = cgas.radius_of_gyration()
        
        # RMSD (after alignment)
        rmsd_all[i] = rms.rmsd(mobile_ca.positions, ref_ca_coords, weights=mobile_ca.masses)
        rmsd_trim[i] = rms.rmsd(trim41.select_atoms("name CA").positions, 
                                 ref_trim_ca.positions, weights=trim41.select_atoms("name CA").masses)
        rmsd_cgas[i] = rms.rmsd(cgas.select_atoms("name CA").positions,
                                 ref_cgas_ca.positions, weights=cgas.select_atoms("name CA").masses)
        
        if i % 500 == 0:
            print(f"  {label}: {i}/{n_sample} frames processed")
    
    return {
        'time': times,
        'com': com_dist,
        'rg_trim': rg_trim,
        'rg_cgas': rg_cgas,
        'rmsd_all': rmsd_all,
        'rmsd_trim': rmsd_trim,
        'rmsd_cgas': rmsd_cgas,
    }

# Run analysis
print("\nAnalyzing GROMACS...")
gmx_data = analyze(u_gmx, trim41_gmx, cgas_gmx, prot_gmx, ref_gmx, "GROMACS", step=10)

print("\nAnalyzing OpenMM...")
omm_data = analyze(u_omm, trim41_omm, cgas_omm, prot_omm, ref_omm, "OpenMM", step=1)

# Save data
for name, data in [("gmx", gmx_data), ("omm", omm_data)]:
    np.savez(os.path.join(OUTDIR, f"{name}_200ns.npz"), **data)

print("\nData saved.")

# Plotting
fig = plt.figure(figsize=(14, 12))
gs = GridSpec(3, 2, figure=fig)

# COM distance
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(gmx_data['time'], gmx_data['com'], alpha=0.7, label='GROMACS 2026', color='tab:blue')
ax1.plot(omm_data['time'], omm_data['com'], alpha=0.7, label='OpenMM', color='tab:orange')
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('COM Distance (Å)')
ax1.set_title('TRIM41-cGAS Center of Mass Distance')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Rg TRIM41
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(gmx_data['time'], gmx_data['rg_trim'], alpha=0.7, label='GROMACS', color='tab:blue')
ax2.plot(omm_data['time'], omm_data['rg_trim'], alpha=0.7, label='OpenMM', color='tab:orange')
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('Rg (Å)')
ax2.set_title('TRIM41 Radius of Gyration')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Rg cGAS
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(gmx_data['time'], gmx_data['rg_cgas'], alpha=0.7, label='GROMACS', color='tab:blue')
ax3.plot(omm_data['time'], omm_data['rg_cgas'], alpha=0.7, label='OpenMM', color='tab:orange')
ax3.set_xlabel('Time (ns)')
ax3.set_ylabel('Rg (Å)')
ax3.set_title('cGAS Radius of Gyration')
ax3.legend()
ax3.grid(True, alpha=0.3)

# RMSD all
ax4 = fig.add_subplot(gs[2, 0])
ax4.plot(gmx_data['time'], gmx_data['rmsd_all'], alpha=0.7, label='GROMACS', color='tab:blue')
ax4.plot(omm_data['time'], omm_data['rmsd_all'], alpha=0.7, label='OpenMM', color='tab:orange')
ax4.set_xlabel('Time (ns)')
ax4.set_ylabel('RMSD (Å)')
ax4.set_title('Protein CA RMSD (aligned to initial)')
ax4.legend()
ax4.grid(True, alpha=0.3)

# RMSD cGAS
ax5 = fig.add_subplot(gs[2, 1])
ax5.plot(gmx_data['time'], gmx_data['rmsd_cgas'], alpha=0.7, label='GROMACS', color='tab:blue')
ax5.plot(omm_data['time'], omm_data['rmsd_cgas'], alpha=0.7, label='OpenMM', color='tab:orange')
ax5.set_xlabel('Time (ns)')
ax5.set_ylabel('RMSD (Å)')
ax5.set_title('cGAS CA RMSD')
ax5.legend()
ax5.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'gmx_vs_openmm_200ns.png'), dpi=150)
print(f"Plot saved to {OUTDIR}/gmx_vs_openmm_200ns.png")

# Statistics
print("\n=== 200ns Statistics ===")
for label, data in [("GROMACS", gmx_data), ("OpenMM", omm_data)]:
    print(f"\n{label}:")
    print(f"  COM distance:    {data['com'].mean():.2f} ± {data['com'].std():.2f} Å  (min={data['com'].min():.2f}, max={data['com'].max():.2f})")
    print(f"  Rg TRIM41:       {data['rg_trim'].mean():.2f} ± {data['rg_trim'].std():.2f} Å")
    print(f"  Rg cGAS:         {data['rg_cgas'].mean():.2f} ± {data['rg_cgas'].std():.2f} Å")
    print(f"  RMSD all:        {data['rmsd_all'].mean():.2f} ± {data['rmsd_all'].std():.2f} Å")
    print(f"  RMSD cGAS:       {data['rmsd_cgas'].mean():.2f} ± {data['rmsd_cgas'].std():.2f} Å")
