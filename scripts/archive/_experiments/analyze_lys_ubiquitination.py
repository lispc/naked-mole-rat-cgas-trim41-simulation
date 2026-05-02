#!/usr/bin/env python3
"""Analyze cGAS lysine residues for ubiquitination potential."""

import sys
import json
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF


def analyze_lys_residues(prmtop, dcd, trim_range, cgas_range, ring_range,
                         system_name, replica_name, outdir="data/analysis"):
    print(f"\n{'='*60}")
    print(f"Analyzing {system_name} {replica_name}")
    print(f"{'='*60}")
    
    u = mda.Universe(prmtop, dcd)
    n_frames = len(u.trajectory)
    print(f"  Frames: {n_frames}")
    
    t_start, t_end = trim_range
    c_start, c_end = cgas_range
    r_start, r_end = ring_range
    
    ring_atoms = u.select_atoms(f"resid {r_start}-{r_end}")
    cgas_lys = u.select_atoms(f"resid {c_start}-{c_end} and resname LYS")
    
    print(f"  TRIM41 RING atoms: {len(ring_atoms)}")
    print(f"  cGAS Lys residues: {len(cgas_lys.residues)}")
    
    lys_resids = list(dict.fromkeys(cgas_lys.residues.resids))
    n_lys = len(lys_resids)
    
    # Align
    print("  Aligning trajectory...")
    ref = mda.Universe(prmtop, dcd)
    aligner = align.AlignTraj(u, ref, select="protein and name CA", in_memory=True)
    aligner.run()
    
    # RMSF
    print("  Computing RMSF...")
    protein_ca = u.select_atoms("protein and name CA")
    rmsf_calc = RMSF(protein_ca).run()
    ca_resids = protein_ca.residues.resids
    resid_to_rmsf_idx = {r: i for i, r in enumerate(ca_resids)}
    
    sasa_data = {r: [] for r in lys_resids}
    dist_data = {r: [] for r in lys_resids}
    
    print(f"  Processing {n_frames} frames...")
    for ts in u.trajectory:
        waters = u.select_atoms("resname WAT and name O")
        for resid in lys_resids:
            lys_atom = u.select_atoms(f"resid {resid} and name NZ")
            if len(lys_atom) == 0:
                continue
            
            dist_mat = mda.lib.distances.distance_array(lys_atom.positions, ring_atoms.positions)
            min_dist = dist_mat.min()
            dist_data[resid].append(min_dist)
            
            if len(waters) > 0:
                water_dists = mda.lib.distances.distance_array(lys_atom.positions, waters.positions)
                n_nearby = (water_dists < 4.0).sum()
                sasa_data[resid].append(n_nearby)
            else:
                sasa_data[resid].append(0)
        
        if ts.frame % 200 == 0 and ts.frame > 0:
            print(f"    Frame {ts.frame}/{n_frames}")
    
    results = []
    for resid in lys_resids:
        if len(dist_data[resid]) == 0:
            continue
        dist_arr = np.array(dist_data[resid])
        sasa_arr = np.array(sasa_data[resid])
        rmsf_val = rmsf_calc.results.rmsf[resid_to_rmsf_idx[resid]] if resid in resid_to_rmsf_idx else 0.0
        
        results.append({
            "resid": int(resid),
            "dist_to_ring_mean": float(np.mean(dist_arr)),
            "dist_to_ring_std": float(np.std(dist_arr)),
            "dist_to_ring_min": float(np.min(dist_arr)),
            "exposure_mean": float(np.mean(sasa_arr)),
            "exposure_std": float(np.std(sasa_arr)),
            "rmsf": float(rmsf_val),
        })
    
    results.sort(key=lambda x: x["dist_to_ring_mean"])
    
    print(f"\n  Top 10 closest Lys to RING:")
    print(f"  {'Resid':>6} {'Dist(Å)':>8} {'Min(Å)':>8} {'RMSF(Å)':>8} {'Exposure':>10}")
    for r in results[:10]:
        print(f"  {r['resid']:>6} {r['dist_to_ring_mean']:>8.1f} {r['dist_to_ring_min']:>8.1f} {r['rmsf']:>8.2f} {r['exposure_mean']:>10.1f}")
    
    output = {
        "system": system_name,
        "replica": replica_name,
        "n_frames": n_frames,
        "n_lys": len(results),
        "ring_range": list(ring_range),
        "lys_data": results,
    }
    
    import os
    os.makedirs(outdir, exist_ok=True)
    outpath = f"{outdir}/{system_name}_{replica_name}_lys_analysis.json"
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n  Saved: {outpath}")
    return output


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--dcd", required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--replica", default="rep1")
    parser.add_argument("--trim-range", nargs=2, type=int, default=[1, 218])
    parser.add_argument("--cgas-range", nargs=2, type=int, required=True)
    parser.add_argument("--ring-range", nargs=2, type=int, default=[1, 43])
    parser.add_argument("--outdir", default="data/analysis/final")
    args = parser.parse_args()
    
    analyze_lys_residues(
        args.prmtop, args.dcd,
        tuple(args.trim_range), tuple(args.cgas_range), tuple(args.ring_range),
        args.system, args.replica, args.outdir
    )
