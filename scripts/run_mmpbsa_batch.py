#!/usr/bin/env python3
"""Batch MM-GBSA for WT and S305-phos (3 reps each)."""

import subprocess
import sys
from pathlib import Path

AMBERHOME = Path.home() / "miniforge3/envs/cgas-md"
MMPBSA = f"{AMBERHOME}/bin/python {AMBERHOME}/bin/MMPBSA.py"
OUTDIR = Path("data/analysis/mmpbsa")

def run_one(name, cp, rp, lp, nc):
    outdir = OUTDIR / name
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Write input
    inp = outdir / "mmpbsa.in"
    with open(inp, 'w') as f:
        f.write("""&general
  startframe=1, endframe=99999, interval=50,
  keep_files=0, use_sander=0,
/
&gb
  igb=5, saltcon=0.150,
/
""")
    
    cmd = [
        sys.executable, str(AMBERHOME / "bin/MMPBSA.py"),
        "-i", str(inp),
        "-o", str(outdir / "results.dat"),
        "-cp", cp,
        "-rp", rp,
        "-lp", lp,
        "-y", nc,
    ]
    
    log = outdir / "mmpbsa.log"
    print(f"\n[START] {name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    with open(log, 'w') as f:
        f.write(result.stdout)
        f.write(result.stderr)
    
    if result.returncode == 0:
        print(f"[DONE] {name}")
    else:
        print(f"[FAIL] {name} (see {log})")
    return result.returncode


def main():
    jobs = []
    for sys_name, sys_path in [
        ("Hsap_WT", "data/md_runs/Hsap_WT"),
        ("Hsap_WT_S305phos", "data/md_runs/Hsap_WT_S305phos"),
    ]:
        cp = f"{sys_path}/{sys_name}_protein.prmtop"
        rp = f"{sys_path}/{sys_name}_trim41.prmtop"
        lp = f"{sys_path}/{sys_name}_cgas.prmtop"
        
        for rep in [1, 2, 3]:
            name = f"{sys_name}_rep{rep}"
            nc = f"data/analysis/mmpbsa/{name}_prot.nc"
            jobs.append((name, cp, rp, lp, nc))
    
    # Run in parallel using multiprocessing
    from multiprocessing import Pool
    
    def run_job(args):
        return run_one(*args)
    
    with Pool(6) as pool:
        results = pool.map(run_job, jobs)
    
    print(f"\n{'='*60}")
    print(f"Batch complete: {sum(1 for r in results if r == 0)}/{len(results)} succeeded")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
