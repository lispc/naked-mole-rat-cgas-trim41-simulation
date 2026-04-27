#!/usr/bin/env python3
"""Analyze Rosetta docking results and compare with previous best poses.

Usage:
    python scripts/analyze_docking_results.py \
        --scorefile structures/docking/rosetta/hgal_WT_global.sc \
        --output-dir structures/docking/rosetta/output_hgal_WT_global \
        --prev-best structures/docking/rosetta/input.pdb \
        --system Hgal_WT
"""
import argparse
import os
import numpy as np
from pathlib import Path


def read_pdb_atoms(filepath):
    """Read CA atoms from PDB."""
    atoms = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain = line[21].strip()
                resi = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[(chain, resi)] = np.array([x, y, z])
    return atoms


def calc_rmsd(atoms1, atoms2):
    common = set(atoms1.keys()) & set(atoms2.keys())
    if len(common) == 0:
        return None
    coords1 = np.array([atoms1[k] for k in sorted(common)])
    coords2 = np.array([atoms2[k] for k in sorted(common)])
    c1 = coords1.mean(axis=0)
    c2 = coords2.mean(axis=0)
    coords1 -= c1
    coords2 -= c2
    H = coords1.T @ coords2
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    coords2_rot = coords2 @ R
    rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2_rot)**2, axis=1)))
    return rmsd


def parse_scorefile(scorefile):
    """Parse Rosetta scorefile and return list of (name, score, I_sc)."""
    results = []
    with open(scorefile) as f:
        for line in f:
            if line.startswith("SCORE:") and "description" not in line:
                parts = line.split()
                if len(parts) < 3:
                    continue
                score = float(parts[1])
                desc = parts[-1]
                # Find I_sc if present
                isc = None
                for i, p in enumerate(parts):
                    if p == "I_sc":
                        if i + 1 < len(parts):
                            isc = float(parts[i + 1])
                        break
                results.append((desc, score, isc))
    return results


def analyze_docking(scorefile, output_dir, prev_best_pdb, system_name):
    results = parse_scorefile(scorefile)
    if not results:
        print(f"No results found in {scorefile}")
        return
    
    # Sort by total_score (lower is better)
    results.sort(key=lambda x: x[1])
    
    print(f"\n{'='*60}")
    print(f"  Rosetta Docking Analysis: {system_name}")
    print(f"{'='*60}")
    print(f"  Total decoys: {len(results)}")
    
    best_name, best_score, best_isc = results[0]
    worst_score = results[-1][1]
    
    print(f"\n  Score range: {best_score:.1f} to {worst_score:.1f} (span={worst_score-best_score:.1f})")
    isc_str = f"{best_isc:.1f}" if best_isc is not None else "N/A"
    print(f"  Best decoy: {best_name} (score={best_score:.1f}, I_sc={isc_str})")
    
    # Top 5
    print(f"\n  Top 5 decoys:")
    for i, (name, score, isc) in enumerate(results[:5]):
        isc_str = f"{isc:.1f}" if isc is not None else "N/A"
        print(f"    {i+1}. {name}: score={score:.1f}, I_sc={isc_str}")
    
    # Compare with previous best
    if prev_best_pdb and os.path.exists(prev_best_pdb):
        best_pdb = os.path.join(output_dir, f"{best_name}.pdb")
        if os.path.exists(best_pdb):
            prev_atoms = read_pdb_atoms(prev_best_pdb)
            new_atoms = read_pdb_atoms(best_pdb)
            rmsd = calc_rmsd(prev_atoms, new_atoms)
            
            print(f"\n  Comparison with previous best pose:")
            print(f"    Previous: {prev_best_pdb}")
            print(f"    New best: {best_pdb}")
            print(f"    CA-RMSD:  {rmsd:.2f} Å")
            
            if rmsd < 2.0:
                verdict = "PASS — Interface modes are essentially identical. Existing MD/US data likely valid."
            elif rmsd < 5.0:
                verdict = "CAUTION — Moderate interface differences. MD may have self-corrected. Assess interface contacts."
            else:
                verdict = "FAIL — Different interface mode. MD/US must be rerun with new pose."
            
            print(f"    Verdict:  {verdict}")
            
            return {
                "system": system_name,
                "n_decoys": len(results),
                "best_score": best_score,
                "best_isc": best_isc,
                "score_span": worst_score - best_score,
                "rmsd_to_prev": rmsd,
                "verdict": verdict,
                "best_pdb": best_pdb,
            }
    
    return {
        "system": system_name,
        "n_decoys": len(results),
        "best_score": best_score,
        "best_isc": best_isc,
        "score_span": worst_score - best_score,
        "best_pdb": os.path.join(output_dir, f"{best_name}.pdb"),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scorefile", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--prev-best", default=None, help="Previous best pose for RMSD comparison")
    parser.add_argument("--system", required=True)
    args = parser.parse_args()
    
    result = analyze_docking(args.scorefile, args.output_dir, args.prev_best, args.system)
    
    if result:
        import json
        out_json = os.path.join(args.output_dir, f"{args.system}_analysis.json")
        with open(out_json, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"\n  Saved: {out_json}")


if __name__ == "__main__":
    main()
