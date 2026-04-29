import time
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

t0 = time.time()
print("Loading universe...")
u = mda.Universe("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd")
ref = mda.Universe("data/md_runs/Hsap_WT/Hsap_WT.prmtop", "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd")
print(f"  Loaded in {time.time()-t0:.1f}s, frames: {len(u.trajectory)}")

t1 = time.time()
print("Aligning...")
align.AlignTraj(u, ref, select="protein and name CA", in_memory=False).run()
print(f"  Aligned in {time.time()-t1:.1f}s")

t2 = time.time()
print("RMSD...")
rms.RMSD(u, ref, select="protein and name CA").run()
print(f"  RMSD in {time.time()-t2:.1f}s")

t3 = time.time()
print("COM distance (all frames)...")
prot = u.select_atoms("protein")
cgas = prot.select_atoms("resid 1-218")
trim = prot.select_atoms("resid 219-541")
com_dists = []
for ts in u.trajectory:
    com_dists.append(np.linalg.norm(cgas.center_of_mass() - trim.center_of_mass()))
print(f"  COM in {time.time()-t3:.1f}s, mean={np.mean(com_dists):.2f}")

print(f"Total: {time.time()-t0:.1f}s")
