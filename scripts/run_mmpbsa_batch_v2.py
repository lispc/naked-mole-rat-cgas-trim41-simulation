#!/usr/bin/env python3
"""Batch MM-GBSA using ThreadPoolExecutor (no pickle issues)."""

import os
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

AMBERHOME = Path.home() / "miniforge3/envs/cgas-md"
OUTDIR = Path("data/analysis/mmpbsa")


def run_one(name, cp, rp, lp, nc):
    # Each job runs in its own working directory to avoid temp file conflicts
    workdir = (OUTDIR / name / "_work").resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    outdir = (OUTDIR / name).resolve()
    
    inp = workdir / "mmpbsa.in"
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
        str(AMBERHOME / "bin/python"),
        str(AMBERHOME / "bin/MMPBSA.py"),
        "-i", str(inp),
        "-o", str(outdir / "results.dat"),
        "-cp", str(Path(cp).resolve()),
        "-rp", str(Path(rp).resolve()),
        "-lp", str(Path(lp).resolve()),
        "-y", str(Path(nc).resolve()),
    ]
    
    log = outdir / "mmpbsa.log"
    print(f"[START] {name}")
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd=str(workdir),
                            env={**os.environ, 'AMBERHOME': str(AMBERHOME)})
    with open(log, 'w') as f:
        f.write(result.stdout)
        f.write(result.stderr)
    
    # Check for DELTA in results
    results_file = outdir / "results.dat"
    if results_file.exists() and "DELTA TOTAL" in results_file.read_text():
        print(f"[DONE] {name}")
        return (name, 0)
    else:
        print(f"[FAIL] {name}")
        return (name, 1)


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
    
    # Run in parallel (6 threads)
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {executor.submit(run_one, *job): job[0] for job in jobs}
        for future in as_completed(futures):
            name, code = future.result()
    
    print(f"\n{'='*60}")
    success = sum(1 for job in jobs if (OUTDIR / job[0] / "results.dat").exists() and "DELTA TOTAL" in (OUTDIR / job[0] / "results.dat").read_text())
    print(f"Batch complete: {success}/{len(jobs)} succeeded")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
