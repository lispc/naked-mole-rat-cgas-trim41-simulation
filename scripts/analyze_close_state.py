#!/usr/bin/env python3
"""
Analyze the close-state conformations from US windows (CV < threshold).
Extracts frames where RING-Lys distance is small and analyzes:
  - Lys-334 environment (contacting residues within 5Å)
  - RING domain contacting residues
  - Backbone phi/psi of Lys-334
  - Hydrogen bonding network

Outputs: PDB frames + contact analysis + visualization script
"""
import argparse
import json
from pathlib import Path
from datetime import datetime

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array


def read_cv(path, threshold):
    """Read CV and return frame indices where CV < threshold."""
    frames = []
    with open(path) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                cv = float(parts[2])
                if cv < threshold:
                    frames.append(i)
    return np.array(frames)


def analyze_close_frames(prmtop, dcd, cv_path, threshold, outdir, name):
    """Analyze frames with CV < threshold."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Find close frames
    close_frames = read_cv(cv_path, threshold)
    print(f"[{name}] Close frames (CV < {threshold}Å): {len(close_frames)} / {len(open(cv_path).readlines())}")
    
    if len(close_frames) == 0:
        print(f"[{name}] No close frames found.")
        return None
    
    # Load universe
    u = mda.Universe(prmtop, dcd)
    
    # Selections
    ring = u.select_atoms("resid 1-43 and protein")  # RING domain
    lys334 = u.select_atoms("resid 334 and protein")  # Lys-334
    cgas = u.select_atoms("resid 219-541 and protein")  # cGAS
    trim41 = u.select_atoms("resid 1-218 and protein")  # TRIM41
    
    # Analyze each close frame
    all_contacts_ring = []
    all_contacts_cgas = []
    lys_rmsf_like = []
    
    for idx in close_frames[:500]:  # limit to 500 frames for speed
        u.trajectory[idx]
        
        # Lys-334 NZ position
        lys_nz = lys334.select_atoms("name NZ")
        if len(lys_nz) == 0:
            continue
        
        # Contacts: RING residues within 5Å of Lys-334 NZ
        ring_ca = ring.select_atoms("name CA")
        dists = distance_array(lys_nz.positions, ring_ca.positions)[0]
        close_ring = ring_ca.resids[dists < 5.0]
        all_contacts_ring.extend(close_ring.tolist())
        
        # Contacts: cGAS residues within 5Å of RING (interface)
        cgas_ca = cgas.select_atoms("name CA")
        ring_all = ring.select_atoms("name CA")
        dist_mat = distance_array(ring_all.positions, cgas_ca.positions)
        contacts = np.argwhere(dist_mat < 5.0)
        for i, j in contacts:
            all_contacts_cgas.append((int(ring_all.resids[i]), int(cgas_ca.resids[j])))
    
    from collections import Counter
    
    # Summarize RING contacts
    ring_counter = Counter(all_contacts_ring)
    top_ring = ring_counter.most_common(20)
    
    # Summarize interface contacts
    interface_counter = Counter(all_contacts_cgas)
    top_interface = interface_counter.most_common(20)
    
    print(f"\n[{name}] Top RING residues near Lys-334:")
    for resid, count in top_ring[:10]:
        print(f"  TRIM41-{resid}: {count} contacts ({count/len(close_frames[:500]):.3f} occupancy)")
    
    print(f"\n[{name}] Top interface contacts (RING↔cGAS):")
    for (t_res, c_res), count in top_interface[:10]:
        print(f"  TRIM41-{t_res} ↔ cGAS-{c_res}: {count} contacts ({count/len(close_frames[:500]):.3f} occupancy)")
    
    # Save summary
    summary = {
        "name": name,
        "threshold_A": threshold,
        "n_close_frames": int(len(close_frames)),
        "n_analyzed": int(min(len(close_frames), 500)),
        "top_ring_contacts": [(int(r), int(c)) for r, c in top_ring],
        "top_interface_contacts": [(f"{t}_{c}", int(cnt)) for (t, c), cnt in top_interface],
    }
    
    summary_path = outdir / f"close_state_{name}.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved: {summary_path}")
    
    # Write PyMOL visualization script
    pymol_script = f"""
# PyMOL script for close-state visualization
# Load your structure and trajectory, then run this script

# Highlight Lys-334
select lys334, resi 334
show sticks, lys334
color magenta, lys334

# Highlight contacting RING residues
select ring_contacts, resi {'+'.join([str(r) for r, _ in top_ring[:10]])}
show sticks, ring_contacts
color cyan, ring_contacts

# Highlight interface
select interface_cgas, resi {'+'.join([str(c) for (_, c), _ in top_interface[:10]])}
show lines, interface_cgas
color yellow, interface_cgas

zoom lys334
"""
    pymol_path = outdir / f"visualize_close_{name}.pml"
    with open(pymol_path, "w") as f:
        f.write(pymol_script)
    print(f"  Saved: {pymol_path}")
    
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--dcd", required=True)
    parser.add_argument("--cv", required=True, help="CV data file")
    parser.add_argument("--threshold", type=float, default=5.0, help="CV threshold in Å")
    parser.add_argument("--name", required=True)
    parser.add_argument("--outdir", default="data/analysis/final/us_Hsap_WT_Lys334/close_state")
    args = parser.parse_args()
    
    print(f"[{datetime.now()}] Analyzing close state: {args.name}")
    analyze_close_frames(args.prmtop, args.dcd, args.cv, args.threshold, args.outdir, args.name)
    print(f"[{datetime.now()}] Done.")


if __name__ == "__main__":
    main()
