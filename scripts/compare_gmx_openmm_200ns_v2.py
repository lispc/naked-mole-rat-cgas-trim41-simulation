"""Compare GROMACS 2026 vs OpenMM 200ns trajectories for Hsap_WT (with manual PBC fix)."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

# Paths
GRO     = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs_gmx2026/Hsap_WT/rep1/prod.gro"
GMX_XTC = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs_gmx2026/Hsap_WT/rep1/prod.xtc"
OMM_TOP = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/Hsap_WT.prmtop"
OMM_DCD = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd"
OUTDIR  = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/analysis/gmx_openmm_comparison"
os.makedirs(OUTDIR, exist_ok=True)

print("Loading trajectories...")
u_gmx = mda.Universe(GRO, GMX_XTC)
u_omm = mda.Universe(OMM_TOP, OMM_DCD)
print(f"  GROMACS: {len(u_gmx.trajectory)} frames")
print(f"  OpenMM:  {len(u_omm.trajectory)} frames")

trim41_gmx = u_gmx.select_atoms("protein and resid 1-218")
cgas_gmx   = u_gmx.select_atoms("protein and resid 219-541")
prot_gmx   = u_gmx.select_atoms("protein")

trim41_omm = u_omm.select_atoms("protein and resid 1-218")
cgas_omm   = u_omm.select_atoms("protein and resid 219-541")
prot_omm   = u_omm.select_atoms("protein")

print(f"  TRIM41: GMX={len(trim41_gmx)}, OMM={len(trim41_omm)}")
print(f"  cGAS:   GMX={len(cgas_gmx)}, OMM={len(cgas_omm)}")

def unwrap_chain(atoms, ca_selection, box):
    """Unwrap a protein chain by correcting PBC jumps between CA atoms."""
    ca = ca_selection
    coords = atoms.positions.copy()
    ca_coords = ca.positions.copy()
    
    # Build residue -> atom indices mapping (local indices within the chain)
    res_to_local = {}
    for i, atom in enumerate(atoms):
        res_idx = atom.resindex
        if res_idx not in res_to_local:
            res_to_local[res_idx] = []
        res_to_local[res_idx].append(i)
    
    # Unwrap CA coordinates
    shifts = np.zeros_like(ca_coords)
    for i in range(1, len(ca_coords)):
        delta = ca_coords[i] + shifts[i-1] - (ca_coords[i-1] + shifts[i-1])
        if box is not None:
            for dim in range(3):
                if box[dim] > 0:
                    while delta[dim] > box[dim] / 2:
                        delta[dim] -= box[dim]
                        shifts[i, dim] -= box[dim]
                    while delta[dim] < -box[dim] / 2:
                        delta[dim] += box[dim]
                        shifts[i, dim] += box[dim]
    
    # Apply shifts to all atoms in each residue
    for i, ca_atom in enumerate(ca):
        res_idx = ca_atom.resindex
        if res_idx in res_to_local:
            for local_idx in res_to_local[res_idx]:
                coords[local_idx] += shifts[i]
    
    return coords

def analyze_gmx(u, trim41, cgas, prot, step=10):
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
    
    # Reference (first frame, unwrapped)
    u.trajectory[0]
    box0 = u.trajectory[0].dimensions[:3] if u.trajectory[0].dimensions is not None else None
    ref_trim_coords = unwrap_chain(trim41, trim41.select_atoms("name CA"), box0)
    ref_cgas_coords = unwrap_chain(cgas, cgas.select_atoms("name CA"), box0)
    
    # Reference CA coords
    ref_trim_ca = ref_trim_coords[trim41.select_atoms("name CA").indices - trim41.indices[0]]
    ref_cgas_ca = ref_cgas_coords[cgas.select_atoms("name CA").indices - cgas.indices[0]]
    ref_all_ca = np.vstack([ref_trim_ca, ref_cgas_ca])
    
    for i, idx in enumerate(indices):
        ts = u.trajectory[idx]
        times[i] = ts.time / 1000.0
        box = ts.dimensions[:3] if ts.dimensions is not None else None
        
        trim_coords = unwrap_chain(trim41, trim41.select_atoms("name CA"), box)
        cgas_coords = unwrap_chain(cgas, cgas.select_atoms("name CA"), box)
        
        # COM distance
        com_trim = np.average(trim_coords, axis=0, weights=trim41.masses)
        com_cgas = np.average(cgas_coords, axis=0, weights=cgas.masses)
        com_dist[i] = np.linalg.norm(com_trim - com_cgas)
        
        # Rg
        rg_trim[i] = np.sqrt(np.average(np.sum((trim_coords - com_trim)**2, axis=1), weights=trim41.masses))
        rg_cgas[i] = np.sqrt(np.average(np.sum((cgas_coords - com_cgas)**2, axis=1), weights=cgas.masses))
        
        # RMSD
        curr_trim_ca = trim_coords[trim41.select_atoms("name CA").indices - trim41.indices[0]]
        curr_cgas_ca = cgas_coords[cgas.select_atoms("name CA").indices - cgas.indices[0]]
        curr_all_ca = np.vstack([curr_trim_ca, curr_cgas_ca])
        
        R, _ = align.rotation_matrix(curr_all_ca, ref_all_ca, weights=np.ones(len(curr_all_ca)))
        curr_all_aligned = curr_all_ca @ R.T
        curr_trim_aligned = curr_trim_ca @ R.T
        curr_cgas_aligned = curr_cgas_ca @ R.T
        
        rmsd_all[i] = rms.rmsd(curr_all_aligned, ref_all_ca, weights=np.ones(len(curr_all_ca)))
        rmsd_trim[i] = rms.rmsd(curr_trim_aligned, ref_trim_ca, weights=np.ones(len(curr_trim_ca)))
        rmsd_cgas[i] = rms.rmsd(curr_cgas_aligned, ref_cgas_ca, weights=np.ones(len(curr_cgas_ca)))
        
        if i % 500 == 0:
            print(f"  GROMACS: {i}/{n_sample} frames")
    
    return {
        'time': times, 'com': com_dist, 'rg_trim': rg_trim, 'rg_cgas': rg_cgas,
        'rmsd_all': rmsd_all, 'rmsd_trim': rmsd_trim, 'rmsd_cgas': rmsd_cgas,
    }

def analyze_omm(u, trim41, cgas, prot):
    n_frames = len(u.trajectory)
    n_sample = n_frames
    
    times = np.zeros(n_sample)
    com_dist = np.zeros(n_sample)
    rg_trim = np.zeros(n_sample)
    rg_cgas = np.zeros(n_sample)
    rmsd_all = np.zeros(n_sample)
    rmsd_trim = np.zeros(n_sample)
    rmsd_cgas = np.zeros(n_sample)
    
    ref = mda.Universe(OMM_TOP, OMM_DCD)
    ref.trajectory[0]
    ref_ca = ref.select_atoms("protein and name CA")
    ref_trim_ca = ref.select_atoms("protein and resid 1-218 and name CA")
    ref_cgas_ca = ref.select_atoms("protein and resid 219-541 and name CA")
    
    for i in range(n_frames):
        ts = u.trajectory[i]
        times[i] = i * 0.1  # OpenMM DCD has 100 ps step
        
        mobile_ca = prot.select_atoms("name CA")
        R, _ = align.rotation_matrix(mobile_ca.positions, ref_ca.positions, weights=mobile_ca.masses)
        prot.atoms.translate(-prot.center_of_mass())
        prot.atoms.rotate(R)
        prot.atoms.translate(ref.select_atoms("protein").center_of_mass())
        
        com_dist[i] = np.linalg.norm(trim41.center_of_mass() - cgas.center_of_mass())
        rg_trim[i] = trim41.radius_of_gyration()
        rg_cgas[i] = cgas.radius_of_gyration()
        
        rmsd_all[i] = rms.rmsd(mobile_ca.positions, ref_ca.positions, weights=mobile_ca.masses)
        rmsd_trim[i] = rms.rmsd(trim41.select_atoms("name CA").positions, ref_trim_ca.positions, weights=trim41.select_atoms("name CA").masses)
        rmsd_cgas[i] = rms.rmsd(cgas.select_atoms("name CA").positions, ref_cgas_ca.positions, weights=cgas.select_atoms("name CA").masses)
        
        if i % 500 == 0:
            print(f"  OpenMM: {i}/{n_sample} frames")
    
    return {
        'time': times, 'com': com_dist, 'rg_trim': rg_trim, 'rg_cgas': rg_cgas,
        'rmsd_all': rmsd_all, 'rmsd_trim': rmsd_trim, 'rmsd_cgas': rmsd_cgas,
    }

print("\nAnalyzing GROMACS (with PBC fix)...")
gmx_data = analyze_gmx(u_gmx, trim41_gmx, cgas_gmx, prot_gmx, step=10)

print("\nAnalyzing OpenMM...")
omm_data = analyze_omm(u_omm, trim41_omm, cgas_omm, prot_omm)

for name, data in [("gmx_v2", gmx_data), ("omm_v2", omm_data)]:
    np.savez(os.path.join(OUTDIR, f"{name}_200ns.npz"), **data)

# Plot
fig = plt.figure(figsize=(14, 12))
gs = GridSpec(3, 2, figure=fig)

ax1 = fig.add_subplot(gs[0, :])
ax1.plot(gmx_data['time'], gmx_data['com'], alpha=0.7, label='GROMACS 2026', color='tab:blue')
ax1.plot(omm_data['time'], omm_data['com'], alpha=0.7, label='OpenMM', color='tab:orange')
ax1.set_xlabel('Time (ns)'); ax1.set_ylabel('COM Distance (Å)')
ax1.set_title('TRIM41-cGAS Center of Mass Distance')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(gmx_data['time'], gmx_data['rg_trim'], alpha=0.7, label='GROMACS', color='tab:blue')
ax2.plot(omm_data['time'], omm_data['rg_trim'], alpha=0.7, label='OpenMM', color='tab:orange')
ax2.set_xlabel('Time (ns)'); ax2.set_ylabel('Rg (Å)')
ax2.set_title('TRIM41 Radius of Gyration')
ax2.legend(); ax2.grid(True, alpha=0.3)

ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(gmx_data['time'], gmx_data['rg_cgas'], alpha=0.7, label='GROMACS', color='tab:blue')
ax3.plot(omm_data['time'], omm_data['rg_cgas'], alpha=0.7, label='OpenMM', color='tab:orange')
ax3.set_xlabel('Time (ns)'); ax3.set_ylabel('Rg (Å)')
ax3.set_title('cGAS Radius of Gyration')
ax3.legend(); ax3.grid(True, alpha=0.3)

ax4 = fig.add_subplot(gs[2, 0])
ax4.plot(gmx_data['time'], gmx_data['rmsd_all'], alpha=0.7, label='GROMACS', color='tab:blue')
ax4.plot(omm_data['time'], omm_data['rmsd_all'], alpha=0.7, label='OpenMM', color='tab:orange')
ax4.set_xlabel('Time (ns)'); ax4.set_ylabel('RMSD (Å)')
ax4.set_title('Protein CA RMSD')
ax4.legend(); ax4.grid(True, alpha=0.3)

ax5 = fig.add_subplot(gs[2, 1])
ax5.plot(gmx_data['time'], gmx_data['rmsd_cgas'], alpha=0.7, label='GROMACS', color='tab:blue')
ax5.plot(omm_data['time'], omm_data['rmsd_cgas'], alpha=0.7, label='OpenMM', color='tab:orange')
ax5.set_xlabel('Time (ns)'); ax5.set_ylabel('RMSD (Å)')
ax5.set_title('cGAS CA RMSD')
ax5.legend(); ax5.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, 'gmx_vs_openmm_200ns_v2.png'), dpi=150)
print(f"\nPlot saved to {OUTDIR}/gmx_vs_openmm_200ns_v2.png")

print("\n=== 200ns Statistics ===")
for label, data in [("GROMACS", gmx_data), ("OpenMM", omm_data)]:
    print(f"\n{label}:")
    print(f"  COM distance:    {data['com'].mean():.2f} ± {data['com'].std():.2f} Å  (min={data['com'].min():.2f}, max={data['com'].max():.2f})")
    print(f"  Rg TRIM41:       {data['rg_trim'].mean():.2f} ± {data['rg_trim'].std():.2f} Å")
    print(f"  Rg cGAS:         {data['rg_cgas'].mean():.2f} ± {data['rg_cgas'].std():.2f} Å")
    print(f"  RMSD all:        {data['rmsd_all'].mean():.2f} ± {data['rmsd_all'].std():.2f} Å")
    print(f"  RMSD cGAS:       {data['rmsd_cgas'].mean():.2f} ± {data['rmsd_cgas'].std():.2f} Å")
