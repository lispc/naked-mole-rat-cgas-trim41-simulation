#!/usr/bin/env python3
"""Launch multiple MD replicas with different random seeds.

Usage:
    python scripts/run_replicas.py \
        --prmtop data/md_runs/Hsap_WT/Hsap_WT.prmtop \
        --pdb data/md_runs/Hsap_WT/Hsap_WT_minimized.pdb \
        --system Hsap_WT \
        --n-replicas 2 \
        --prod-ns 200 \
        --gpu-start 0
"""
import argparse
import os
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description="Launch MD replicas")
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--n-replicas", type=int, default=2)
    parser.add_argument("--prod-ns", type=int, default=200)
    parser.add_argument("--gpu-start", type=int, default=0,
                        help="First GPU ID to use (0-3)")
    parser.add_argument("--outdir-base", default="data/md_runs")
    args = parser.parse_args()
    
    outdir = os.path.join(args.outdir_base, args.system)
    os.makedirs(outdir, exist_ok=True)
    
    n_gpus = 4
    
    for rep in range(1, args.n_replicas + 1):
        rep_name = f"{args.system}_rep{rep + 1}"  # rep2, rep3, ...
        gpu_id = (args.gpu_start + rep - 1) % n_gpus
        seed = int(time.time()) + rep * 1000
        
        cmd = [
            "python", "scripts/run_md.py",
            "--prmtop", args.prmtop,
            "--pdb", args.pdb,
            "--name", rep_name,
            "--outdir", outdir,
            "--prod-ns", str(args.prod_ns),
            "--gpu", str(gpu_id),
            "--seed", str(seed),
        ]
        
        log_file = os.path.join(outdir, f"{rep_name}.log")
        
        print(f"Launching {rep_name} on GPU {gpu_id} (seed={seed})")
        print(f"  Log: {log_file}")
        
        with open(log_file, 'w') as log:
            subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        
        time.sleep(2)  # stagger starts
    
    print(f"\nLaunched {args.n_replicas} replicas for {args.system}")
    print(f"Monitor with: tail -f {outdir}/{args.system}_rep*.log")


if __name__ == "__main__":
    main()
