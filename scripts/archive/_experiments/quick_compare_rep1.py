#!/usr/bin/env python3
"""Quick comparison of Hsap_WT rep1 vs Hsap_4mut rep1 (200ns each)."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("data/analysis/initial_comparison")
outdir.mkdir(parents=True, exist_ok=True)

systems = {
    "Hsap_WT": ("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd"),
    "Hsap_4mut": ("data/md_runs/Hsap_4mut/Hsap_4mut.prmtop", "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd"),
}

results = {}

for name, (top, traj) in systems.items():
    print(f"\n=== {name} ===")
    u = mda.Universe(top, traj)
    ref = mda.Universe(top, traj)
    
    prot = u.select_atoms("protein")
    ref_prot = ref.select_atoms("protein")
    
    print(f"  Protein atoms: {len(prot)}, residues: {len(prot.residues)}")
    
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    rmsd = rms.RMSD(u, ref, select="protein and name CA").run()
    time_ns = rmsd.results.rmsd[:, 1] / 1000.0
    rmsd_vals = rmsd.results.rmsd[:, 2]
    
    print(f"  Frames: {len(u.trajectory)}")
    print(f"  CA RMSD: mean={rmsd_vals.mean():.2f}Å, std={rmsd_vals.std():.2f}Å, max={rmsd_vals.max():.2f}Å")
    
    rmsf = rms.RMSF(prot.select_atoms("name CA")).run()
    rmsf_vals = rmsf.results.rmsf
    ca_resids = prot.select_atoms("name CA").resids
    
    print(f"  CA RMSF: mean={rmsf_vals.mean():.2f}Å, max={rmsf_vals.max():.2f}Å (res {ca_resids[np.argmax(rmsf_vals)]})")
    
    ca_atoms = prot.select_atoms("name CA")
    ca_indices = ca_atoms.indices
    cgas_ca = prot.select_atoms("name CA and resid 1-218")
    trim_ca = prot.select_atoms("name CA and resid 219-541")
    
    # Map global atom indices to local CA array indices
    ca_idx_map = {idx: i for i, idx in enumerate(ca_indices)}
    cgas_local = [ca_idx_map[idx] for idx in cgas_ca.indices]
    trim_local = [ca_idx_map[idx] for idx in trim_ca.indices]
    
    cgas_rmsf = rmsf_vals[cgas_local]
    trim_rmsf = rmsf_vals[trim_local]
    
    print(f"  cGAS CA RMSF: mean={cgas_rmsf.mean():.2f}Å, max={cgas_rmsf.max():.2f}Å")
    print(f"  TRIM41 CA RMSF: mean={trim_rmsf.mean():.2f}Å, max={trim_rmsf.max():.2f}Å")
    
    cgas = prot.select_atoms("resid 1-218")
    trim = prot.select_atoms("resid 219-541")
    
    com_dists = []
    for ts in u.trajectory:
        com_cgas = cgas.center_of_mass()
        com_trim = trim.center_of_mass()
        dist = np.linalg.norm(com_cgas - com_trim)
        com_dists.append(dist)
    
    com_dists = np.array(com_dists)
    print(f"  COM distance: mean={com_dists.mean():.2f}Å, std={com_dists.std():.2f}Å, range=[{com_dists.min():.2f}, {com_dists.max():.2f}]")
    
    rg_vals = []
    for ts in u.trajectory:
        rg_vals.append(prot.radius_of_gyration())
    rg_vals = np.array(rg_vals)
    print(f"  Rg: mean={rg_vals.mean():.2f}Å, std={rg_vals.std():.2f}Å")
    
    results[name] = {
        "rmsd": rmsd_vals,
        "rmsf": rmsf_vals,
        "ca_resids": ca_resids,
        "com_dists": com_dists,
        "rg": rg_vals,
        "time_ns": time_ns,
    }

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

ax = axes[0, 0]
for name, data in results.items():
    ax.plot(data["time_ns"], data["rmsd"], label=name, lw=0.5, alpha=0.8)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("CA RMSD (Å)")
ax.set_title("Backbone RMSD vs t=0")
ax.legend()

ax = axes[0, 1]
for name, data in results.items():
    ax.plot(data["time_ns"], data["com_dists"], label=name, lw=0.5, alpha=0.8)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("COM distance (Å)")
ax.set_title("cGAS-TRIM41 COM distance")
ax.legend()

ax = axes[1, 0]
for name, data in results.items():
    ax.plot(data["time_ns"], data["rg"], label=name, lw=0.5, alpha=0.8)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Rg (Å)")
ax.set_title("Radius of Gyration")
ax.legend()

ax = axes[1, 1]
for name, data in results.items():
    ax.plot(data["ca_resids"], data["rmsf"], label=name, lw=0.8, alpha=0.8)
ax.axvline(218.5, color="gray", ls="--", alpha=0.5, label="cGAS/TRIM41 boundary")
ax.set_xlabel("Residue ID (prmtop numbering)")
ax.set_ylabel("CA RMSF (Å)")
ax.set_title("Per-residue RMSF")
ax.legend()

fig.tight_layout()
fig.savefig(outdir / "comparison_rep1.png", dpi=300)
plt.close(fig)

print("\n" + "="*60)
print("SUMMARY: Hsap_WT rep1 vs Hsap_4mut rep1 (200ns)")
print("="*60)
print(f"{'Metric':<30} {'Hsap_WT':>12} {'Hsap_4mut':>12}")
print("-"*60)
for metric in ["rmsd", "com_dists", "rg"]:
    wt = results["Hsap_WT"][metric]
    mut = results["Hsap_4mut"][metric]
    print(f"{metric + ' mean (Å)':<30} {wt.mean():>12.2f} {mut.mean():>12.2f}")
    print(f"{metric + ' std (Å)':<30} {wt.std():>12.2f} {mut.std():>12.2f}")

print(f"\nPlots saved to {outdir / 'comparison_rep1.png'}")
