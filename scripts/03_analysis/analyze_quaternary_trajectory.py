#!/usr/bin/env python3
"""
Trajectory analysis for quaternary MVP MD runs.

Uses prmtop-based topology with known resid mapping:
  - K315  = resid 718
  - Ub G76 = resid 707
  - E2 K85 = resid 765
"""

import numpy as np
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis import distances, align
import json

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")

def analyze_trajectory(prmtop_path, dcd_path, label, out_csv):
    """Analyze a single trajectory."""
    print(f"\n{'='*60}")
    print(f"Analyzing {label}")
    print(f"  prmtop: {prmtop_path}")
    print(f"  dcd: {dcd_path}")
    print(f"{'='*60}")

    u = mda.Universe(str(prmtop_path), str(dcd_path))
    n_frames = len(u.trajectory)
    print(f"Total frames: {n_frames}")

    # Reference structure for E2~Ub RMSD (5FER)
    u_fer = mda.Universe(str(BASE / "data/structures/quaternary_mvp/5FER.pdb"))
    fer_e2_ca = u_fer.select_atoms("chainID B and name CA")
    fer_ub_ca = u_fer.select_atoms("chainID C and name CA")

    # Selections in trajectory
    k315_nz = u.select_atoms("resid 718 and name NZ")
    ub_g76_c = u.select_atoms("resid 707 and name C")
    e2_k85_nz = u.select_atoms("resid 765 and name NZ")
    protein = u.select_atoms("protein")
    e2_ca = u.select_atoms("resid 681-765 and name CA")  # approximate E2 range
    ub_ca = u.select_atoms("resid 632-707 and name CA")  # approximate Ub range

    print(f"Selected atoms: K315 NZ={len(k315_nz)}, Ub G76 C={len(ub_g76_c)}, E2 K85 NZ={len(e2_k85_nz)}")
    print(f"E2 CA approx: {len(e2_ca)}, Ub CA approx: {len(ub_ca)}")

    results = {
        'time_ps': [],
        'd_k315_ub': [],
        'd_e2_ub': [],
        'd_k315_e2': [],
        'rg_total': [],
        'rmsd_e2': [],
        'rmsd_ub': [],
    }

    for ts in u.trajectory:
        t = ts.time  # ps
        
        # Key distances
        d1 = distances.distance_array(k315_nz.positions, ub_g76_c.positions)[0, 0]
        d2 = distances.distance_array(e2_k85_nz.positions, ub_g76_c.positions)[0, 0]
        d3 = distances.distance_array(k315_nz.positions, e2_k85_nz.positions)[0, 0]
        
        # Rg
        rg = protein.radius_of_gyration()
        
        # E2~Ub RMSD vs 5FER (no alignment, just raw distance)
        # Align E2 CA to 5FER E2 CA
        def calc_rmsd(pos1, pos2):
            n = min(len(pos1), len(pos2))
            diff = pos1[:n] - pos2[:n]
            return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        if len(e2_ca) > 0 and len(fer_e2_ca) > 0:
            rmsd_e2 = calc_rmsd(e2_ca.positions, fer_e2_ca.positions)
        else:
            rmsd_e2 = np.nan
            
        if len(ub_ca) > 0 and len(fer_ub_ca) > 0:
            rmsd_ub = calc_rmsd(ub_ca.positions, fer_ub_ca.positions)
        else:
            rmsd_ub = np.nan

        results['time_ps'].append(t)
        results['d_k315_ub'].append(d1)
        results['d_e2_ub'].append(d2)
        results['d_k315_e2'].append(d3)
        results['rg_total'].append(rg)
        results['rmsd_e2'].append(rmsd_e2)
        results['rmsd_ub'].append(rmsd_ub)

        if ts.frame % 50 == 0:
            print(f"  Frame {ts.frame}/{n_frames} ({t/1000:.1f} ns): K315->Ub={d1:.2f}Å, E2->Ub={d2:.2f}Å, Rg={rg:.1f}Å")

    # Save CSV
    import csv
    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time_ps', 'd_k315_ub', 'd_e2_ub', 'd_k315_e2', 'rg_total', 'rmsd_e2', 'rmsd_ub'])
        for i in range(len(results['time_ps'])):
            writer.writerow([results['time_ps'][i], results['d_k315_ub'][i],
                           results['d_e2_ub'][i], results['d_k315_e2'][i],
                           results['rg_total'][i], results['rmsd_e2'][i], results['rmsd_ub'][i]])

    print(f"\nSaved results: {out_csv}")
    
    # Summary statistics
    print(f"\nSummary for {label}:")
    print(f"  K315->Ub:  mean={np.mean(results['d_k315_ub']):.2f} Å, std={np.std(results['d_k315_ub']):.2f} Å, min={np.min(results['d_k315_ub']):.2f} Å, max={np.max(results['d_k315_ub']):.2f} Å")
    print(f"  E2->Ub:    mean={np.mean(results['d_e2_ub']):.2f} Å, std={np.std(results['d_e2_ub']):.2f} Å")
    print(f"  K315->E2:  mean={np.mean(results['d_k315_e2']):.2f} Å, std={np.std(results['d_k315_e2']):.2f} Å")
    print(f"  Rg:        mean={np.mean(results['rg_total']):.2f} Å, std={np.std(results['rg_total']):.2f} Å")
    print(f"  E2 RMSD:   mean={np.nanmean(results['rmsd_e2']):.2f} Å")
    print(f"  Ub RMSD:   mean={np.nanmean(results['rmsd_ub']):.2f} Å")


def main():
    # WT
    analyze_trajectory(
        BASE / "data/structures/quaternary_mvp/quaternary_mvp.prmtop",
        BASE / "data/md_runs/quaternary_mvp/WT_rep1/quaternary_mvp_WT_rep1.dcd",
        "WT_rep1",
        BASE / "data/md_runs/quaternary_mvp/WT_rep1/analysis.csv"
    )

    # 4mut
    analyze_trajectory(
        BASE / "data/structures/quaternary_mvp/quaternary_mvp_4mut.prmtop",
        BASE / "data/md_runs/quaternary_mvp/4mut_rep1/quaternary_mvp_4mut_rep1.dcd",
        "4mut_rep1",
        BASE / "data/md_runs/quaternary_mvp/4mut_rep1/analysis.csv"
    )


if __name__ == "__main__":
    main()
