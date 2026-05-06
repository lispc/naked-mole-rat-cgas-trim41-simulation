#!/usr/bin/env python3
"""
Monitor Hgal MD completion and auto-launch next batch + final analysis.
"""
import subprocess
import time
import sys
from pathlib import Path
from datetime import datetime

REPS = ["rep1", "rep2", "rep3"]
LOG_FILES = [f"data/md_runs/Hgal_domain/{r}/Hgal_domain_{r}_prod.log" for r in REPS]
TARGET_NS = 200.0
CHECK_INTERVAL = 300  # 5 minutes

def check_progress():
    """Returns (all_done, progress_dict)"""
    progress = {}
    all_done = True
    for rep, log in zip(REPS, LOG_FILES):
        try:
            with open(log) as f:
                last = f.readlines()[-1].strip()
            # Parse: " 100000000   200.000     -1478000.00"
            parts = last.split()
            time_ns = float(parts[1])
            progress[rep] = time_ns
            if time_ns < TARGET_NS:
                all_done = False
        except Exception as e:
            progress[rep] = f"error: {e}"
            all_done = False
    return all_done, progress

def run_cmd(cmd, name, logfile):
    """Run a command in background, logging to file."""
    print(f"[{datetime.now()}] Starting: {name}")
    with open(logfile, "w") as f:
        proc = subprocess.Popen(
            cmd, shell=True, stdout=f, stderr=subprocess.STDOUT,
            executable="/bin/bash"
        )
    print(f"  PID: {proc.pid}, log: {logfile}")
    return proc

def launch_final_analysis():
    """Run full Hgal WT analysis with 200ns data."""
    cmd = """source ~/miniforge3/etc/profile.d/conda.sh && conda activate cgas-md && \
python scripts/analyze_system.py \
  --system Hgal_WT \
  --prmtop data/md_runs/Hgal_domain/Hgal_domain.prmtop \
  --trajectories \
    data/md_runs/Hgal_domain/rep1/Hgal_domain_rep1_prod.dcd \
    data/md_runs/Hgal_domain/rep2/Hgal_domain_rep2_prod.dcd \
    data/md_runs/Hgal_domain/rep3/Hgal_domain_rep3_prod.dcd \
  --replica-names rep1 rep2 rep3 \
  --cgas-range 219 573 \
  --active-sites '{"S463": 482, "E511": 530, "Y527": 546, "T530": 549}' \
  --dt-ns 0.1 \
  --outdir data/analysis/final_200ns"""
    return run_cmd(cmd, "Hgal_WT_final_analysis", "data/analysis/final_200ns/analysis.log")

def launch_new_md():
    """Launch next batch of MD simulations."""
    cmds = []
    
    # Hsap WT on GPU 0
    cmds.append(("""source ~/miniforge3/etc/profile.d/conda.sh && conda activate cgas-md && \
CUDA_VISIBLE_DEVICES=0 python scripts/run_md.py \
  --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
  --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
  --name Hsap_WT_rep1 --outdir data/md_runs/Hsap_WT/rep1 \
  --prod-ns 200 --platform CUDA""", "Hsap_WT_MD", "data/md_runs/Hsap_WT/rep1/md_launch.log"))
    
    # Hsap 4mut on GPU 1
    cmds.append(("""source ~/miniforge3/etc/profile.d/conda.sh && conda activate cgas-md && \
CUDA_VISIBLE_DEVICES=1 python scripts/run_md.py \
  --prmtop data/md_runs/Hsap_4mut/Hsap_4mut.prmtop \
  --pdb data/md_runs/Hsap_4mut/Hsap_4mut_minimized.pdb \
  --name Hsap_4mut_rep1 --outdir data/md_runs/Hsap_4mut/rep1 \
  --prod-ns 200 --platform CUDA""", "Hsap_4mut_MD", "data/md_runs/Hsap_4mut/rep1/md_launch.log"))
    
    # Hgal 4mut_rev on GPU 2
    cmds.append(("""source ~/miniforge3/etc/profile.d/conda.sh && conda activate cgas-md && \
CUDA_VISIBLE_DEVICES=2 python scripts/run_md.py \
  --prmtop data/md_runs/Hgal_4mut_rev/Hgal_4mut_rev.prmtop \
  --pdb data/md_runs/Hgal_4mut_rev/Hgal_4mut_rev_minimized.pdb \
  --name Hgal_4mut_rev_rep1 --outdir data/md_runs/Hgal_4mut_rev/rep1 \
  --prod-ns 200 --platform CUDA""", "Hgal_4mut_rev_MD", "data/md_runs/Hgal_4mut_rev/rep1/md_launch.log"))
    
    procs = []
    for cmd, name, log in cmds:
        Path(log).parent.mkdir(parents=True, exist_ok=True)
        procs.append(run_cmd(cmd, name, log))
    return procs

def main():
    print(f"[{datetime.now()}] Monitor started. Checking every {CHECK_INTERVAL}s.")
    print(f"Target: {TARGET_NS} ns for all {len(REPS)} replicas.")
    
    while True:
        all_done, progress = check_progress()
        print(f"[{datetime.now()}] Progress: {progress}")
        
        if all_done:
            print(f"[{datetime.now()}] ALL REPLICAS COMPLETE!")
            
            # Launch final analysis
            launch_final_analysis()
            
            # Launch new MD batch
            launch_new_md()
            
            print(f"[{datetime.now()}] All launched. Monitor exiting.")
            break
        
        time.sleep(CHECK_INTERVAL)

if __name__ == "__main__":
    main()
