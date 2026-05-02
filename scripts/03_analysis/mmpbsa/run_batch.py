#!/usr/bin/env python3
"""
Fix MM-GBSA by re-running with -rp and -lp to get Delta binding energies.
"""

import os
import sys
import subprocess
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
AMBERHOME = Path("/home/scroll/miniforge3/envs/cgas-md")
MMPBSA_PY = AMBERHOME / "bin" / "MMPBSA.py"
OUTDIR = BASE / "data" / "analysis" / "mmpbsa"

JOBS = [
    ("Hsap_WT_rep2",
     BASE / "data/md_runs/Hsap_WT/Hsap_WT_protein.prmtop",
     BASE / "data/md_runs/Hsap_WT/Hsap_WT_receptor.prmtop",
     BASE / "data/md_runs/Hsap_WT/Hsap_WT_ligand.prmtop",
     BASE / "data/analysis/mmpbsa/Hsap_WT_rep2_prot.nc"),
    ("Hsap_WT_rep3",
     BASE / "data/md_runs/Hsap_WT/Hsap_WT_protein.prmtop",
     BASE / "data/md_runs/Hsap_WT/Hsap_WT_receptor.prmtop",
     BASE / "data/md_runs/Hsap_WT/Hsap_WT_ligand.prmtop",
     BASE / "data/analysis/mmpbsa/Hsap_WT_rep3_prot.nc"),
    ("Hsap_WT_S305phos_rep1",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_protein.prmtop",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_receptor.prmtop",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_ligand.prmtop",
     BASE / "data/analysis/mmpbsa/Hsap_WT_S305phos_rep1_prot.nc"),
    ("Hsap_WT_S305phos_rep2",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_protein.prmtop",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_receptor.prmtop",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_ligand.prmtop",
     BASE / "data/analysis/mmpbsa/Hsap_WT_S305phos_rep2_prot.nc"),
    ("Hsap_WT_S305phos_rep3",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_protein.prmtop",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_receptor.prmtop",
     BASE / "data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_ligand.prmtop",
     BASE / "data/analysis/mmpbsa/Hsap_WT_S305phos_rep3_prot.nc"),
]


def run_one(name, cp, rp, lp, nc):
    workdir = OUTDIR / f"{name}_delta_work"
    workdir.mkdir(parents=True, exist_ok=True)

    mmpbsa_in = workdir / "mmpbsa.in"
    mmpbsa_in.write_text("""&general
startframe=1, endframe=99999, interval=50,
keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=1,
print_res="all",
/
""")

    results = OUTDIR / f"{name}_delta_results.dat"
    decomp = OUTDIR / f"{name}_delta_decomp.dat"
    log = workdir / "mmpbsa.log"

    cmd = [
        sys.executable, str(MMPBSA_PY),
        "-i", str(mmpbsa_in),
        "-o", str(results),
        "-do", str(decomp),
        "-cp", str(cp), "-rp", str(rp), "-lp", str(lp),
        "-y", str(nc),
    ]

    env = {**dict(os.environ), "AMBERHOME": str(AMBERHOME)}
    env["OMP_NUM_THREADS"] = "4"

    print(f"[{name}] Starting...")
    print(f"[{name}] {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(workdir), env=env)

    log.write_text(result.stdout)
    if result.stderr:
        log.write_text(log.read_text() + "\n--- STDERR ---\n" + result.stderr)

    if result.returncode == 0:
        print(f"[{name}] ✅ Done")
        return True
    else:
        print(f"[{name}] ❌ Failed (exit {result.returncode})")
        err_tail = result.stderr[-500:] if len(result.stderr) > 500 else result.stderr
        print(f"[{name}] STDERR tail: {err_tail}")
        return False


def main():
    os.environ.setdefault("AMBERHOME", str(AMBERHOME))

    results = []
    for name, cp, rp, lp, nc in JOBS:
        if not all(f.exists() for f in [cp, rp, lp, nc]):
            print(f"[{name}] Skipping — missing input file")
            results.append((name, False))
            continue
        ok = run_one(name, cp, rp, lp, nc)
        results.append((name, ok))

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for name, ok in results:
        status = "✅ PASS" if ok else "❌ FAIL"
        print(f"  {status}  {name}")


if __name__ == "__main__":
    main()
