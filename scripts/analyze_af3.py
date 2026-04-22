#!/usr/bin/env python3
"""Analyze AlphaFold3 prediction results: extract ipTM, pTM, pLDDT, interface residues."""
import json
import sys
from pathlib import Path
import numpy as np

JOBS = {
    "job1_Hsap_WT": {"cgas": "Hsap_cGAS_WT", "trim41": "TRIM41_WT"},
    "job2_Hsap_4mut": {"cgas": "Hsap_cGAS_4mut", "trim41": "TRIM41_WT"},
    "job3_Hgal_WT": {"cgas": "Hgal_cGAS_WT", "trim41": "TRIM41_WT"},
    "job4_Hgal_4mut_rev": {"cgas": "Hgal_cGAS_4mut_rev", "trim41": "TRIM41_WT"},
}

# Mutation positions (1-based) in each sequence
MUT_POS = {
    "Hsap_cGAS_WT": [463, 479, 495, 498],
    "Hsap_cGAS_4mut": [463, 479, 495, 498],
    "Hgal_cGAS_WT": [463, 511, 527, 530],
    "Hgal_cGAS_4mut_rev": [463, 511, 527, 530],
}

# cGAS chain lengths (for chain assignment)
CGAS_LEN = {
    "Hsap_cGAS_WT": 522,
    "Hsap_cGAS_4mut": 522,
    "Hgal_cGAS_WT": 554,
    "Hgal_cGAS_4mut_rev": 554,
}


def parse_pdb(pdb_path):
    """Parse AF3 PDB, extract CA atoms with B-factors (pLDDT) and chain IDs."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain = line[21]
                resi = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                bfactor = float(line[60:66])
                atoms.append({"chain": chain, "resi": resi, "x": x, "y": y, "z": z, "plddt": bfactor})
    return atoms


def compute_interface_residues(atoms, cgas_len, cutoff=5.0):
    """Find interface residues within cutoff Angstrom between chains."""
    chain_atoms = {}
    for a in atoms:
        chain_atoms.setdefault(a["chain"], []).append(a)
    
    chains = sorted(chain_atoms.keys())
    if len(chains) < 2:
        return set(), set()
    
    # Assume first chain = cGAS, second = TRIM41 (AF3 outputs A, B, ...)
    cgas_chain = chains[0]
    trim_chain = chains[1]
    
    cgas_interface = set()
    trim_interface = set()
    
    for a in chain_atoms[cgas_chain]:
        for b in chain_atoms[trim_chain]:
            dist = np.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)
            if dist <= cutoff:
                cgas_interface.add(a["resi"])
                trim_interface.add(b["resi"])
    
    return cgas_interface, trim_interface


def analyze_job(job_dir: Path):
    conf_path = job_dir / "confidence.json"
    pdb_path = job_dir / "ranked_0.pdb"
    
    if not conf_path.exists():
        print(f"  [SKIP] {job_dir.name}: confidence.json not found")
        return None
    if not pdb_path.exists():
        print(f"  [SKIP] {job_dir.name}: ranked_0.pdb not found")
        return None
    
    with open(conf_path) as f:
        conf = json.load(f)
    
    # AF3 server format
    iptm = conf.get("iptm", None)
    ptm = conf.get("ptm", None)
    
    atoms = parse_pdb(pdb_path)
    job_name = job_dir.name
    cgas_name = JOBS[job_name]["cgas"]
    cgas_len = CGAS_LEN[cgas_name]
    mut_pos = MUT_POS[cgas_name]
    
    cgas_iface, trim_iface = compute_interface_residues(atoms, cgas_len)
    
    # pLDDT of interface residues
    cgas_iface_plddt = [a["plddt"] for a in atoms if a["chain"] == "A" and a["resi"] in cgas_iface]
    trim_iface_plddt = [a["plddt"] for a in atoms if a["chain"] == "B" and a["resi"] in trim_iface]
    
    # pLDDT of mutation residues
    mut_plddt = [a["plddt"] for a in atoms if a["chain"] == "A" and a["resi"] in mut_pos]
    
    # Are mutations at interface?
    mut_at_iface = [p in cgas_iface for p in mut_pos]
    
    result = {
        "job": job_name,
        "iptm": iptm,
        "ptm": ptm,
        "cgas_interface_residues": sorted(cgas_iface),
        "trim41_interface_residues": sorted(trim_iface),
        "cgas_interface_size": len(cgas_iface),
        "trim41_interface_size": len(trim_iface),
        "cgas_iface_mean_plddt": np.mean(cgas_iface_plddt) if cgas_iface_plddt else 0,
        "trim41_iface_mean_plddt": np.mean(trim_iface_plddt) if trim_iface_plddt else 0,
        "mutation_positions": mut_pos,
        "mutation_plddt": mut_plddt,
        "mutations_at_interface": mut_at_iface,
    }
    return result


def main():
    af3_dir = Path("structures/af3_raw")
    print("=" * 60)
    print("AlphaFold3 Prediction Quality Report")
    print("=" * 60)
    
    all_results = []
    for job_name in JOBS:
        job_dir = af3_dir / job_name
        if not job_dir.exists():
            print(f"\n{job_name}: DIRECTORY NOT FOUND")
            continue
        
        print(f"\n--- {job_name} ---")
        res = analyze_job(job_dir)
        if res is None:
            continue
        all_results.append(res)
        
        print(f"  ipTM: {res['iptm']:.3f}" if res['iptm'] else "  ipTM: N/A")
        print(f"  pTM:  {res['ptm']:.3f}" if res['ptm'] else "  pTM: N/A")
        print(f"  Interface size: cGAS={res['cgas_interface_size']}, TRIM41={res['trim41_interface_size']}")
        print(f"  Interface mean pLDDT: cGAS={res['cgas_iface_mean_plddt']:.1f}, TRIM41={res['trim41_iface_mean_plddt']:.1f}")
        print(f"  Mutation pLDDT: {res['mutation_plddt']}")
        print(f"  Mutations at interface (5A): {res['mutations_at_interface']}")
        
        # Quality verdict
        if res['iptm'] and res['iptm'] >= 0.70:
            print(f"  ✅ ACCEPT: ipTM >= 0.70")
        elif res['iptm'] and res['iptm'] >= 0.60:
            print(f"  ⚠️  MARGINAL: ipTM 0.60-0.70")
        else:
            print(f"  ❌ REJECT: ipTM < 0.60")
    
    # Summary table
    print("\n" + "=" * 60)
    print("Summary Table")
    print("=" * 60)
    print(f"{'Job':<25} {'ipTM':>8} {'pTM':>8} {'Interface':>10} {'Status':>10}")
    print("-" * 65)
    for r in all_results:
        status = "ACCEPT" if r['iptm'] and r['iptm'] >= 0.70 else ("MARGINAL" if r['iptm'] and r['iptm'] >= 0.60 else "REJECT")
        print(f"{r['job']:<25} {r['iptm']:>8.3f} {r['ptm']:>8.3f} {r['cgas_interface_size']:>5}+{r['trim41_interface_size']:>4} {status:>10}")
    
    # Save JSON
    import json as json_mod
    out_path = Path("data/analysis/af3_analysis.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json_mod.dump(all_results, f, indent=2)
    print(f"\nDetailed results saved to: {out_path}")


if __name__ == "__main__":
    main()
