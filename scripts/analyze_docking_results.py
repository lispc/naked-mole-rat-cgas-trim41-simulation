#!/usr/bin/env python3
"""
Analyze Rosetta docking results:
1. Parse scorefiles to find best decoy (lowest total_score)
2. Calculate CA-RMSD between new best pose and previous best pose
3. Classify: PASS (<2Å), GRAY (2-5Å), FAIL (>5Å)

Usage:
  python scripts/analyze_docking_results.py
"""
import numpy as np
from pathlib import Path


# ---------------------------------------------------------------------------
# Kabsch RMSD (no Biopython/MDAnalysis dependency)
# ---------------------------------------------------------------------------

def parse_pdb_coords(pdb_path, atom_name="CA", chain=None):
    """Extract (x,y,z) coordinates for specified atom type from PDB."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            if len(line) < 54:
                continue
            rec_atom = line[12:16].strip()
            rec_chain = line[21:22].strip() if len(line) > 21 else ""
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if rec_atom == atom_name:
                if chain is None or rec_chain == chain:
                    coords.append([x, y, z])
    return np.array(coords, dtype=float)


def kabsch_rmsd(P, Q):
    """Calculate RMSD between two sets of points after optimal rotation."""
    # Center
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    # Covariance matrix
    H = Pc.T @ Qc
    # SVD
    U, S, Vt = np.linalg.svd(H)
    # Rotation matrix
    d = np.linalg.det(Vt.T @ U.T)
    if d < 0:
        Vt[-1, :] *= -1
    R = Vt.T @ U.T
    # Rotate P
    P_rot = Pc @ R
    # RMSD
    rmsd = np.sqrt(((P_rot - Qc) ** 2).sum() / len(P))
    return float(rmsd)


def ca_rmsd(path_a, path_b):
    """CA-RMSD between two PDBs (all chains)."""
    coords_a = parse_pdb_coords(path_a, "CA")
    coords_b = parse_pdb_coords(path_b, "CA")
    if len(coords_a) != len(coords_b):
        n = min(len(coords_a), len(coords_b))
        print(f"  Warning: mismatched CA counts ({len(coords_a)} vs {len(coords_b)}), using first {n}")
        coords_a = coords_a[:n]
        coords_b = coords_b[:n]
    return kabsch_rmsd(coords_a, coords_b)


# ---------------------------------------------------------------------------
# Scorefile parsing
# ---------------------------------------------------------------------------

def parse_scorefile(sc_path):
    """Parse Rosetta scorefile, return list of (total_score, desc) sorted by score."""
    entries = []
    with open(sc_path) as f:
        for line in f:
            if line.startswith("SCORE:") and "description" not in line:
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    total_score = float(parts[1])
                    desc = parts[-1]
                    entries.append((total_score, desc))
                except ValueError:
                    continue
    entries.sort(key=lambda x: x[0])
    return entries


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_system(name, new_sc, new_pdb_dir, old_sc_or_pdb, old_is_pdb=False):
    print(f"\n{'='*60}")
    print(f"System: {name}")
    print(f"{'='*60}")

    # New docking
    new_entries = parse_scorefile(new_sc)
    if not new_entries:
        print("  ERROR: No entries in new scorefile")
        return None
    new_best_score, new_best_desc = new_entries[0]
    print(f"  New best: {new_best_desc}  total_score={new_best_score:.3f}  (n={len(new_entries)})")

    # Top 5 new
    print("  New top 5:")
    for score, desc in new_entries[:5]:
        marker = " <<< BEST" if desc == new_best_desc else ""
        print(f"    {desc}: {score:.3f}{marker}")

    new_best_pdb = Path(new_pdb_dir) / f"{new_best_desc}.pdb"
    if not new_best_pdb.exists():
        # Try with _0001 suffix
        alt = Path(new_pdb_dir) / f"{new_best_desc}_0001.pdb"
        if alt.exists():
            new_best_pdb = alt
        else:
            print(f"  ERROR: PDB not found: {new_best_pdb}")
            return None

    # Old docking
    if old_is_pdb:
        old_best_pdb = Path(old_sc_or_pdb)
        old_best_score = None
        old_best_desc = Path(old_sc_or_pdb).name
        print(f"  Old best: {old_best_desc}  (LightDock best_pose)")
    else:
        old_entries = parse_scorefile(old_sc_or_pdb)
        if not old_entries:
            print("  ERROR: No entries in old scorefile")
            return None
        old_best_score, old_best_desc = old_entries[0]
        print(f"  Old best: {old_best_desc}  total_score={old_best_score:.3f}  (n={len(old_entries)})")
        old_pdb_dir = Path(old_sc_or_pdb).parent
        old_best_pdb = old_pdb_dir / f"{old_best_desc}.pdb"
        if not old_best_pdb.exists():
            alt = old_pdb_dir / f"{old_best_desc}_0001.pdb"
            if alt.exists():
                old_best_pdb = alt
            else:
                print(f"  ERROR: Old PDB not found: {old_best_pdb}")
                return None

    # CA-RMSD
    print(f"\n  Calculating CA-RMSD...")
    print(f"    New: {new_best_pdb}")
    print(f"    Old: {old_best_pdb}")
    rmsd = ca_rmsd(str(new_best_pdb), str(old_best_pdb))

    # Verdict
    if rmsd < 2.0:
        verdict = "PASS"
        action = "Keep existing MD data"
    elif rmsd < 5.0:
        verdict = "GRAY"
        action = "Analyze interface contacts; may need MD rerun"
    else:
        verdict = "FAIL"
        action = "Rerun MD with new pose"

    print(f"\n  >>> CA-RMSD = {rmsd:.2f} Å  →  {verdict}")
    print(f"  >>> Action: {action}")

    return {
        "system": name,
        "new_best": new_best_desc,
        "new_score": float(new_best_score),
        "old_best": old_best_desc,
        "old_score": float(old_best_score) if old_best_score else None,
        "rmsd": rmsd,
        "verdict": verdict,
        "action": action,
    }


def main():
    base = Path("structures/docking/rosetta")

    systems = [
        {
            "name": "Hgal_WT",
            "new_sc": base / "hgal_WT_global.sc",
            "new_pdb_dir": base / "output_hgal_WT_global",
            "old_sc_or_pdb": "structures/docking/lightdock/Hgal_domain/best_pose.pdb",
            "old_is_pdb": True,
        },
        {
            "name": "Hgal_4mut_rev",
            "new_sc": base / "hgal_4mut_rev_global.sc",
            "new_pdb_dir": base / "output_hgal_4mut_rev_global",
            "old_sc_or_pdb": base / "output_global/global.sc",
            "old_is_pdb": False,
        },
        {
            "name": "Hsap_WT",
            "new_sc": base / "hsap_WT_global.sc",
            "new_pdb_dir": base / "output_hsap_WT_global",
            "old_sc_or_pdb": base / "output_hsap_global/hsap_global.sc",
            "old_is_pdb": False,
        },
        {
            "name": "Hsap_4mut",
            "new_sc": base / "hsap_4mut_global.sc",
            "new_pdb_dir": base / "output_hsap_4mut_global",
            "old_sc_or_pdb": base / "hsap_4mut_global.sc",  # Same file? Need to check
            "old_is_pdb": False,
        },
    ]

    # Special case: Hsap_4mut old may have its own scorefile
    old_hsap_4mut_sc = base / "output_hsap_4mut_global/hsap_4mut_global.sc"
    if old_hsap_4mut_sc.exists():
        systems[3]["old_sc_or_pdb"] = old_hsap_4mut_sc
    else:
        # Old Hsap_4mut was already nstruct=100, but we need to check if it's the same run
        # For now, compare new best within the same scorefile (self-consistency check)
        print("WARNING: Hsap_4mut old scorefile not found. Will compare new best vs 2nd best as sanity check.")
        systems[3]["old_sc_or_pdb"] = base / "hsap_4mut_global.sc"

    results = []
    for sys in systems:
        if sys["name"] == "Hsap_4mut" and not old_hsap_4mut_sc.exists():
            # Self-comparison: compare best vs 2nd best
            print(f"\n{'='*60}")
            print(f"System: Hsap_4mut (self-comparison: best vs 2nd best)")
            print(f"{'='*60}")
            entries = parse_scorefile(sys["new_sc"])
            if len(entries) >= 2:
                best_score, best_desc = entries[0]
                second_score, second_desc = entries[1]
                best_pdb = sys["new_pdb_dir"] / f"{best_desc}.pdb"
                second_pdb = sys["new_pdb_dir"] / f"{second_desc}.pdb"
                if not best_pdb.exists():
                    alt = sys["new_pdb_dir"] / f"{best_desc}_0001.pdb"
                    if alt.exists():
                        best_pdb = alt
                if not second_pdb.exists():
                    alt = sys["new_pdb_dir"] / f"{second_desc}_0001.pdb"
                    if alt.exists():
                        second_pdb = alt
                if best_pdb.exists() and second_pdb.exists():
                    rmsd = ca_rmsd(str(best_pdb), str(second_pdb))
                    print(f"  Best: {best_desc} ({best_score:.3f})")
                    print(f"  2nd:  {second_desc} ({second_score:.3f})")
                    print(f"  Self CA-RMSD = {rmsd:.2f} Å")
                    results.append({
                        "system": "Hsap_4mut",
                        "new_best": best_desc,
                        "new_score": best_score,
                        "old_best": second_desc,
                        "old_score": second_score,
                        "rmsd": rmsd,
                        "verdict": "SELF_CHECK",
                        "action": "Old docking was already nstruct=100; self-consistency check only",
                    })
            continue

        res = analyze_system(
            sys["name"],
            sys["new_sc"],
            sys["new_pdb_dir"],
            sys["old_sc_or_pdb"],
            sys["old_is_pdb"],
        )
        if res:
            results.append(res)

    # Summary table
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'System':<16} {'RMSD(Å)':<10} {'Verdict':<8} {'Action'}")
    print("-" * 60)
    for r in results:
        print(f"{r['system']:<16} {r['rmsd']:<10.2f} {r['verdict']:<8} {r['action']}")

    # Save JSON
    import json
    out = Path("data/analysis/docking_comparison_0428.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
