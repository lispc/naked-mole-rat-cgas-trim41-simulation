#!/usr/bin/env python3
"""Compare Boltz-2 truncated prediction with AF3 reference structure."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from pathlib import Path

# Create merged AF3 PDB if not exists
af3_merged = Path("data/boltz_test_truncated/af3_merged_clean.pdb")

# Load structures
u_b = mda.Universe("data/boltz_test_truncated/boltz2_model.pdb")
u_a = mda.Universe(str(af3_merged))

print(f"Boltz-2: {len(u_b.atoms)} atoms")
print(f"AF3: {len(u_a.atoms)} atoms")

# Select CA atoms
b_cgas = u_b.select_atoms("chainID A and name CA")
b_trim = u_b.select_atoms("chainID B and name CA")
a_cgas = u_a.select_atoms("chainID A and name CA")
a_trim = u_a.select_atoms("chainID B and name CA")

print(f"cGAS CA: Boltz={len(b_cgas)}, AF3={len(a_cgas)}")
print(f"TRIM41 CA: Boltz={len(b_trim)}, AF3={len(a_trim)}")

# Kabsch alignment on cGAS CA
def kabsch(m, r):
    mc = m - m.mean(axis=0)
    rc = r - r.mean(axis=0)
    H = mc.T @ rc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1,:] *= -1
        R = Vt.T @ U.T
    t = r.mean(axis=0) - R @ m.mean(axis=0)
    return R, t

R, t = kabsch(b_cgas.positions, a_cgas.positions)
u_b.atoms.positions = (u_b.atoms.positions @ R.T) + t

# Re-select after transform
b_cgas = u_b.select_atoms("chainID A and name CA")
b_trim = u_b.select_atoms("chainID B and name CA")

cgas_rmsd = np.sqrt(np.mean((b_cgas.positions - a_cgas.positions)**2))
trim_rmsd = np.sqrt(np.mean((b_trim.positions - a_trim.positions)**2))
print(f"\ncGAS CA RMSD: {cgas_rmsd:.2f} Å")
print(f"TRIM41 CA RMSD: {trim_rmsd:.2f} Å")

# Per-residue distances
cgas_d = np.linalg.norm(b_cgas.positions - a_cgas.positions, axis=1)
trim_d = np.linalg.norm(b_trim.positions - a_trim.positions, axis=1)

# Interface contacts
def get_contacts(u, c1, c2, cutoff=5.0):
    A = u.select_atoms(c1)
    B = u.select_atoms(c2)
    contacts = set()
    for ra in A.residues:
        for rb in B.residues:
            d = np.min(np.linalg.norm(ra.atoms.positions[:,None] - rb.atoms.positions[None,:], axis=2))
            if d < cutoff:
                contacts.add((ra.resid, rb.resid))
    return contacts

b_contacts = get_contacts(u_b, "chainID A", "chainID B")
a_contacts = get_contacts(u_a, "chainID A", "chainID B")
shared = b_contacts & a_contacts
print(f"\nBoltz contacts: {len(b_contacts)}")
print(f"AF3 contacts: {len(a_contacts)}")
print(f"Shared: {len(shared)}")
print(f"Jaccard: {len(shared)/len(b_contacts | a_contacts):.3f}")

# Active sites
# Active sites: Boltz uses 1-based numbering, AF3 uses original numbering
# D431 = resid 232 in truncated (431-199=232)
sites_boltz = {"D431": 232, "K479": 280, "L495": 296, "K498": 299}
sites_af3 = {"D431": 431, "K479": 479, "L495": 495, "K498": 498}
print("\nActive site distances:")
for name in sites_boltz:
    b_site = u_b.select_atoms(f"chainID A and resid {sites_boltz[name]}")
    a_site = u_a.select_atoms(f"chainID A and resid {sites_af3[name]}")
    b_t = u_b.select_atoms("chainID B")
    a_t = u_a.select_atoms("chainID B")
    b_d = np.min(np.linalg.norm(b_site.positions[:,None] - b_t.positions[None,:], axis=2))
    a_d = np.min(np.linalg.norm(a_site.positions[:,None] - a_t.positions[None,:], axis=2))
    print(f"  {name}: Boltz={b_d:.1f}Å, AF3={a_d:.1f}Å, Δ={b_d-a_d:+.1f}Å")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax = axes[0,0]
ax.plot(b_cgas.resids, cgas_d, 'b-', lw=0.8, label='cGAS CT')
ax.axhline(2, color='g', ls='--', alpha=0.5)
ax.axhline(5, color='r', ls='--', alpha=0.5)
ax.set_xlabel('Residue')
ax.set_ylabel('CA distance (Å)')
ax.set_title('cGAS CT: per-residue CA distance')
ax.legend()

ax = axes[0,1]
ax.plot(b_trim.resids, trim_d, 'r-', lw=0.8, label='TRIM41 SPRY')
ax.axhline(2, color='g', ls='--', alpha=0.5)
ax.axhline(5, color='r', ls='--', alpha=0.5)
ax.set_xlabel('Residue')
ax.set_ylabel('CA distance (Å)')
ax.set_title('TRIM41 SPRY: per-residue CA distance')
ax.legend()

ax = axes[1,0]
ax.bar(['Shared', 'Boltz-only', 'AF3-only'], 
       [len(shared), len(b_contacts-a_contacts), len(a_contacts-b_contacts)], 
       color=['green', 'blue', 'orange'], alpha=0.7)
ax.set_ylabel('Contact pairs')
jac = len(shared)/len(b_contacts | a_contacts)
ax.set_title(f'Interface contacts (Jaccard={jac:.2f})')

ax = axes[1,1]
names = list(sites_boltz.keys())
b_vals = [np.min(np.linalg.norm(u_b.select_atoms(f"chainID A and resid {sites_boltz[n]}").positions[:,None] - u_b.select_atoms("chainID B").positions[None,:], axis=2)) for n in names]
a_vals = [np.min(np.linalg.norm(u_a.select_atoms(f"chainID A and resid {sites_af3[n]}").positions[:,None] - u_a.select_atoms("chainID B").positions[None,:], axis=2)) for n in names]
x = np.arange(len(names))
ax.bar(x-0.2, b_vals, 0.4, label='Boltz-2', color='blue', alpha=0.7)
ax.bar(x+0.2, a_vals, 0.4, label='AF3', color='orange', alpha=0.7)
ax.set_xticks(x)
ax.set_xticklabels(names)
ax.set_ylabel('Distance to TRIM41 (Å)')
ax.set_title('Active site distances')
ax.legend()

plt.tight_layout()
outdir = Path("data/boltz_test_truncated/comparison")
outdir.mkdir(exist_ok=True)
plt.savefig(outdir / "boltz2_vs_af3_comparison.png", dpi=300)
print(f"\nPlot saved: {outdir / 'boltz2_vs_af3_comparison.png'}")

# Assessment
print("\n=== Overall Assessment ===")
if trim_rmsd < 5:
    print(f"TRIM41 RMSD {trim_rmsd:.2f}Å: VERY similar placement")
elif trim_rmsd < 10:
    print(f"TRIM41 RMSD {trim_rmsd:.2f}Å: MODERATELY similar")
else:
    print(f"TRIM41 RMSD {trim_rmsd:.2f}Å: DIVERGENT placement")

if jac > 0.7:
    print(f"Jaccard {jac:.2f}: HIGHLY conserved interface")
elif jac > 0.4:
    print(f"Jaccard {jac:.2f}: PARTIALLY conserved interface")
else:
    print(f"Jaccard {jac:.2f}: DIFFERENT interface")
