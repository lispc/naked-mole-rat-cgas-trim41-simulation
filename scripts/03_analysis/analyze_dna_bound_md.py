#!/usr/bin/env python3
"""Quick analysis of DNA-bound cGAS+SPRY MD trajectory."""
import json
import warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.distances import distance_array

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

BASE = Path(__file__).resolve().parent.parent.parent
OUTDIR = BASE / "data/analysis/dna_bound_md"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Load trajectory ────────────────────────────────────────────────────
pdb_path = BASE / "data/md_runs/cgas_dna_spry_protein/minimized.pdb"
dcd_path = BASE / "data/md_runs/cgas_dna_spry_protein/WT_rep1/cgas_dna_spry_WT_prod.dcd"

u = mda.Universe(str(pdb_path), str(dcd_path))
print(f"Atoms: {len(u.atoms)}, Residues: {len(u.residues)}")
print(f"Frames: {len(u.trajectory)}")

# Identify chains
for seg in sorted(set(u.atoms.segids)):
    ca = u.select_atoms(f"segid {seg} and name CA")
    if len(ca) > 0:
        print(f"  segid {seg}: {len(ca)} CA, resid {ca.resids[0]}-{ca.resids[-1]}")

# The PDBFixer output has renumbered residues. Identify cGAS vs SPRY by size.
# cGAS = larger chain (~468 CA), SPRY = smaller (~218 CA)
ca_by_seg = {}
for seg in sorted(set(u.atoms.segids)):
    ca = u.select_atoms(f"segid {seg} and name CA")
    if len(ca) > 100:
        ca_by_seg[seg] = len(ca)

# Sort by size: largest = cGAS, second = SPRY
sorted_segs = sorted(ca_by_seg.items(), key=lambda x: -x[1])
cgas_seg = sorted_segs[0][0]
spry_seg = sorted_segs[1][0]
print(f"\ncGAS: segid {cgas_seg} ({ca_by_seg[cgas_seg]} CA)")
print(f"SPRY: segid {spry_seg} ({ca_by_seg[spry_seg]} CA)")

# Align on cGAS CA
cgas_ca = u.select_atoms(f"segid {cgas_seg} and name CA")
u.trajectory[0]
ref_pos = cgas_ca.positions.copy()
align.AlignTraj(u, u, select=f"segid {cgas_seg} and name CA", in_memory=True).run()

# ── Analyze ────────────────────────────────────────────────────────────
n_frames = len(u.trajectory)
stride = max(1, n_frames // 200)  # ~200 frames for analysis
frame_indices = list(range(0, n_frames, stride))
n_sample = len(frame_indices)

times_ns = np.zeros(n_sample)
rmsd_cgas = np.zeros(n_sample)
rmsd_spry = np.zeros(n_sample)
com_dist = np.zeros(n_sample)
n_contacts = np.zeros(n_sample)

spry_ca = u.select_atoms(f"segid {spry_seg} and name CA")

for s, fi in enumerate(frame_indices):
    u.trajectory[fi]
    times_ns[s] = u.trajectory.time / 1000.0

    # RMSD from frame 0
    rmsd_cgas[s] = rms.rmsd(cgas_ca.positions, ref_pos, superposition=False)
    rmsd_spry[s] = rms.rmsd(spry_ca.positions, spry_ca.positions, superposition=False)

    # COM distance
    com_dist[s] = np.linalg.norm(cgas_ca.center_of_mass() - spry_ca.center_of_mass())

    # Interface contacts: CA-CA < 8 Å between cGAS and SPRY
    dist_mat = distance_array(cgas_ca.positions, spry_ca.positions)
    n_contacts[s] = np.sum(dist_mat < 8.0)

print(f"\nTrajectory summary (first half vs second half):")
mid = n_sample // 2
for label, lo, hi in [("First 25ns", 0, mid), ("Last 25ns", mid, n_sample)]:
    print(f"  {label}:")
    print(f"    RMSD cGAS: {np.mean(rmsd_cgas[lo:hi]):.1f} ± {np.std(rmsd_cgas[lo:hi]):.1f} Å")
    print(f"    RMSD SPRY: {np.mean(rmsd_spry[lo:hi]):.1f} ± {np.std(rmsd_spry[lo:hi]):.1f} Å")
    print(f"    COM dist:  {np.mean(com_dist[lo:hi]):.1f} ± {np.std(com_dist[lo:hi]):.1f} Å")
    print(f"    Contacts:  {np.mean(n_contacts[lo:hi]):.0f} ± {np.std(n_contacts[lo:hi]):.0f}")

# ── K315 position ──────────────────────────────────────────────────────
# Find K315 by sequence context in the cGAS chain
# K315 → in the Boltz-2 PDBFixer output, it should be a LYS somewhere in the chain
# Find it by checking residues near the middle of cGAS
cgas_lys = u.select_atoms(f"segid {cgas_seg} and resname LYS and name NZ")
print(f"\nK315 search: {len(cgas_lys)} LYS in cGAS")

# Compute distance from each LYS to SPRY COM at frame 0 and last frame
u.trajectory[0]
spry_com_0 = spry_ca.center_of_mass()
lys_dists_0 = []
for atom in cgas_lys:
    d = np.linalg.norm(atom.position - spry_com_0)
    lys_dists_0.append((atom.resid, d))

u.trajectory[-1]
spry_com_end = spry_ca.center_of_mass()
lys_dists_end = []
for atom in cgas_lys:
    d = np.linalg.norm(atom.position - spry_com_end)
    lys_dists_end.append((atom.resid, d))

# Sort by distance
lys_dists_0.sort(key=lambda x: x[1])
lys_dists_end.sort(key=lambda x: x[1])

print(f"\nClosest LYS to SPRY COM:")
print(f"  Frame 0:  ", end="")
for res, d in lys_dists_0[:5]:
    print(f"{res}:{d:.1f}Å ", end="")
print(f"\n  Frame 50ns:", end="")
for res, d in lys_dists_end[:5]:
    print(f"{res}:{d:.1f}Å ", end="")
print()

# ── Plots ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0, 0].plot(times_ns, rmsd_cgas, lw=0.5, alpha=0.8, color="steelblue", label="cGAS")
axes[0, 0].plot(times_ns, rmsd_spry, lw=0.5, alpha=0.8, color="coral", label="SPRY")
axes[0, 0].set_xlabel("Time (ns)")
axes[0, 0].set_ylabel("RMSD from frame 0 (Å)")
axes[0, 0].set_title("Backbone RMSD")
axes[0, 0].legend()

axes[0, 1].plot(times_ns, com_dist, lw=0.5, color="steelblue")
axes[0, 1].axhline(np.mean(com_dist), color="red", ls="--", alpha=0.5)
axes[0, 1].set_xlabel("Time (ns)")
axes[0, 1].set_ylabel("COM distance (Å)")
axes[0, 1].set_title(f"cGAS-SPRY COM distance (mean={np.mean(com_dist):.1f} Å)")

axes[1, 0].plot(times_ns, n_contacts, lw=0.5, color="steelblue")
axes[1, 0].axhline(np.mean(n_contacts), color="red", ls="--", alpha=0.5)
axes[1, 0].set_xlabel("Time (ns)")
axes[1, 0].set_ylabel("CA-CA contacts (<8 Å)")
axes[1, 0].set_title(f"Interface contacts (mean={np.mean(n_contacts):.0f})")

# Contact map at first and last frame
u.trajectory[0]
dist_0 = distance_array(cgas_ca.positions, spry_ca.positions)
u.trajectory[-1]
dist_end = distance_array(cgas_ca.positions, spry_ca.positions)
contact_0 = (dist_0 < 8.0).astype(float)
contact_end = (dist_end < 8.0).astype(float)

ax = axes[1, 1]
# Show contact difference: red = lost, blue = gained
diff = contact_end - contact_0
ax.imshow(diff.T, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1,
          extent=[cgas_ca.resids[0], cgas_ca.resids[-1],
                  spry_ca.resids[0], spry_ca.resids[-1]])
ax.set_xlabel("cGAS residue")
ax.set_ylabel("SPRY residue")
ax.set_title("Contact change (red=lost, blue=gained)")

fig.tight_layout()
fig.savefig(OUTDIR / "dna_bound_summary.png", dpi=300)
plt.close(fig)
print(f"\nSaved: {OUTDIR / 'dna_bound_summary.png'}")

# ── Export ─────────────────────────────────────────────────────────────
results = {
    "n_frames": n_frames,
    "rmsd_cgas_mean": float(np.mean(rmsd_cgas)),
    "rmsd_spry_mean": float(np.mean(rmsd_spry)),
    "com_dist_mean": float(np.mean(com_dist)),
    "com_dist_first_half": float(np.mean(com_dist[:mid])),
    "com_dist_second_half": float(np.mean(com_dist[mid:])),
    "contacts_mean": float(np.mean(n_contacts)),
    "contacts_first_half": float(np.mean(n_contacts[:mid])),
    "contacts_second_half": float(np.mean(n_contacts[mid:])),
    "closest_lys_frame0": [(int(r), float(d)) for r, d in lys_dists_0[:5]],
    "closest_lys_frame50": [(int(r), float(d)) for r, d in lys_dists_end[:5]],
}
with open(OUTDIR / "dna_bound_results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved: {OUTDIR / 'dna_bound_results.json'}")
print("Done.")
