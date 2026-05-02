"""Final corrected comparison plot."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

OUTDIR = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/analysis/gmx_openmm_comparison"

# GROMACS RMSD from gmx rms (whole xtc, CA fit, CA calc)
gmx_rmsd = []
with open("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs_gmx2026/Hsap_WT/rep1/gmx_rmsd_whole_ca.xvg") as f:
    for line in f:
        if line.startswith('#') or line.startswith('@'):
            continue
        parts = line.strip().split()
        if len(parts) >= 2:
            gmx_rmsd.append([float(parts[0]), float(parts[1])])
gmx_rmsd = np.array(gmx_rmsd)

# COM distances (no alignment, already saved)
com_data = np.load(f"{OUTDIR}/com_distance_noalign.npz")
t_gmx = com_data['t_gmx']
d_gmx = com_data['d_gmx']
t_omm = com_data['t_omm']
d_omm = com_data['d_omm']

# OpenMM RMSD (corrected)
print("Computing OpenMM correct RMSD...")
u = mda.Universe("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/Hsap_WT.prmtop",
                  "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd")
ref = mda.Universe("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/Hsap_WT.prmtop",
                    "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd")
ref.trajectory[0]
ref_ca = ref.select_atoms("protein and name CA")
omm_rmsd = np.zeros(len(u.trajectory))
for i in range(len(u.trajectory)):
    u.trajectory[i]
    align.alignto(u, ref, select="protein and name CA", weights="mass")
    mobile_ca = u.select_atoms("protein and name CA")
    omm_rmsd[i] = rms.rmsd(mobile_ca.positions, ref_ca.positions, weights=mobile_ca.masses)
    if i % 500 == 0:
        print(f"  RMSD {i}/{len(u.trajectory)}")
omm_time = np.arange(len(u.trajectory)) * 0.1

# OpenMM Rg
print("Computing OpenMM Rg...")
u2 = mda.Universe("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/Hsap_WT.prmtop",
                   "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd")
rg_t_omm = np.zeros(len(u2.trajectory))
rg_c_omm = np.zeros(len(u2.trajectory))
for i in range(len(u2.trajectory)):
    u2.trajectory[i]
    rg_t_omm[i] = u2.select_atoms("protein and resid 1-218").radius_of_gyration()
    rg_c_omm[i] = u2.select_atoms("protein and resid 219-541").radius_of_gyration()
    if i % 500 == 0:
        print(f"  Rg {i}/{len(u2.trajectory)}")

# GROMACS Rg from v3 npz (no alignment, correct)
gmx_data = np.load(f"{OUTDIR}/gmx_whole_200ns.npz")
rg_t_gmx = gmx_data['rg_trim']
rg_c_gmx = gmx_data['rg_cgas']

# Plot
fig = plt.figure(figsize=(14, 14))
gs = GridSpec(4, 2, figure=fig)

ax1 = fig.add_subplot(gs[0, :])
ax1.plot(t_gmx, d_gmx, alpha=0.7, label='GROMACS 2026', color='tab:blue', lw=0.8)
ax1.plot(t_omm, d_omm, alpha=0.7, label='OpenMM', color='tab:orange', lw=0.8)
ax1.axhline(y=60, color='red', linestyle='--', alpha=0.3, label='Dissociation threshold (~60 Å)')
ax1.set_xlabel('Time (ns)'); ax1.set_ylabel('COM Distance (Å)')
ax1.set_title('TRIM41-cGAS Center of Mass Distance')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2 = fig.add_subplot(gs[1, :])
ax2.plot(gmx_rmsd[:,0], gmx_rmsd[:,1]*10, alpha=0.7, label='GROMACS 2026', color='tab:blue', lw=0.8)
ax2.plot(omm_time, omm_rmsd, alpha=0.7, label='OpenMM', color='tab:orange', lw=0.8)
ax2.set_xlabel('Time (ns)'); ax2.set_ylabel('RMSD (Å)')
ax2.set_title('Protein CA RMSD (corrected alignment)')
ax2.legend(); ax2.grid(True, alpha=0.3)

ax3 = fig.add_subplot(gs[2, 0])
ax3.plot(t_gmx, rg_t_gmx, alpha=0.7, label='GROMACS', color='tab:blue', lw=0.8)
ax3.plot(omm_time, rg_t_omm, alpha=0.7, label='OpenMM', color='tab:orange', lw=0.8)
ax3.set_xlabel('Time (ns)'); ax3.set_ylabel('Rg (Å)')
ax3.set_title('TRIM41 Radius of Gyration')
ax3.legend(); ax3.grid(True, alpha=0.3)

ax4 = fig.add_subplot(gs[2, 1])
ax4.plot(t_gmx, rg_c_gmx, alpha=0.7, label='GROMACS', color='tab:blue', lw=0.8)
ax4.plot(omm_time, rg_c_omm, alpha=0.7, label='OpenMM', color='tab:orange', lw=0.8)
ax4.set_xlabel('Time (ns)'); ax4.set_ylabel('Rg (Å)')
ax4.set_title('cGAS Radius of Gyration')
ax4.legend(); ax4.grid(True, alpha=0.3)

ax5 = fig.add_subplot(gs[3, 0])
mask_gmx = t_gmx <= 150
mask_omm = t_omm <= 150
ax5.plot(t_gmx[mask_gmx], d_gmx[mask_gmx], alpha=0.7, label='GROMACS', color='tab:blue', lw=0.8)
ax5.plot(t_omm[mask_omm], d_omm[mask_omm], alpha=0.7, label='OpenMM', color='tab:orange', lw=0.8)
ax5.set_xlabel('Time (ns)'); ax5.set_ylabel('COM Distance (Å)')
ax5.set_title('COM Distance (0-150 ns)')
ax5.legend(); ax5.grid(True, alpha=0.3)

ax6 = fig.add_subplot(gs[3, 1])
mask_gmx_rms = gmx_rmsd[:,0] <= 150
ax6.plot(gmx_rmsd[mask_gmx_rms,0], gmx_rmsd[mask_gmx_rms,1]*10, alpha=0.7, label='GROMACS', color='tab:blue', lw=0.8)
ax6.plot(omm_time[:1500], omm_rmsd[:1500], alpha=0.7, label='OpenMM', color='tab:orange', lw=0.8)
ax6.set_xlabel('Time (ns)'); ax6.set_ylabel('RMSD (Å)')
ax6.set_title('RMSD (0-150 ns)')
ax6.legend(); ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/gmx_vs_openmm_200ns_FINAL.png", dpi=150)
print(f"\nFinal plot saved.")

print("\n" + "="*70)
print("FINAL CORRECTED STATISTICS")
print("="*70)
print("\n【COM Distance (no alignment)】")
print(f"  GROMACS:  {d_gmx.mean():.2f} ± {d_gmx.std():.2f} Å  (min={d_gmx.min():.2f}, max={d_gmx.max():.2f})")
print(f"  OpenMM:   {d_omm.mean():.2f} ± {d_omm.std():.2f} Å  (min={d_omm.min():.2f}, max={d_omm.max():.2f})")
print(f"  GROMACS 0-150ns: {d_gmx[t_gmx<=150].mean():.2f} ± {d_gmx[t_gmx<=150].std():.2f} Å")
print(f"  OpenMM  0-150ns: {d_omm[t_omm<=150].mean():.2f} ± {d_omm[t_omm<=150].std():.2f} Å")

print("\n【RMSD (CA, correct alignment)】")
print(f"  GROMACS:  {gmx_rmsd[:,1].mean()*10:.2f} ± {gmx_rmsd[:,1].std()*10:.2f} Å  (max={gmx_rmsd[:,1].max()*10:.2f})")
print(f"  OpenMM:   {omm_rmsd.mean():.2f} ± {omm_rmsd.std():.2f} Å  (max={omm_rmsd.max():.2f})")
gmx_rmsd_150 = gmx_rmsd[gmx_rmsd[:,0] <= 150]
print(f"  GROMACS 0-150ns: {gmx_rmsd_150[:,1].mean()*10:.2f} ± {gmx_rmsd_150[:,1].std()*10:.2f} Å")
print(f"  OpenMM  0-150ns: {omm_rmsd[:1500].mean():.2f} ± {omm_rmsd[:1500].std():.2f} Å")

print("\n【Rg TRIM41】")
print(f"  GROMACS:  {rg_t_gmx.mean():.2f} ± {rg_t_gmx.std():.2f} Å")
print(f"  OpenMM:   {rg_t_omm.mean():.2f} ± {rg_t_omm.std():.2f} Å")
print("\n【Rg cGAS】")
print(f"  GROMACS:  {rg_c_gmx.mean():.2f} ± {rg_c_gmx.std():.2f} Å")
print(f"  OpenMM:   {rg_c_omm.mean():.2f} ± {rg_c_omm.std():.2f} Å")
