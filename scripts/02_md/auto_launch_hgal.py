#!/usr/bin/env python3
"""
Auto-launcher for Hgal NEW system MD replicas.
Monitors existing MD jobs and launches new ones on idle GPUs.

Background task: launches remaining 5 reps when GPUs become available.
"""

import os
import subprocess
import time
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
PYTHON = "/home/scroll/miniforge3/envs/cgas-md/bin/python"
RUN_SCRIPT = BASE / "scripts/02_md/run_production.py"

QUEUE = [
    ("Hgal_WT", 2),
    ("Hgal_WT", 3),
    ("Hgal_4mut_rev", 1),
    ("Hgal_4mut_rev", 2),
    ("Hgal_4mut_rev", 3),
]


def get_running_md_pids():
    """Get all running MD python processes."""
    result = subprocess.run(["pgrep", "-f", "run_production.py"],
                          capture_output=True, text=True)
    pids = []
    for pid_str in result.stdout.strip().split("\n"):
        if pid_str.strip():
            try:
                pids.append(int(pid_str.strip()))
            except ValueError:
                pass
    return pids


def get_pid_gpu(pid):
    """Get GPU ID for a given PID."""
    result = subprocess.run(["nvidia-smi", "pmon", "-s", "um", "-c", "1"],
                          capture_output=True, text=True)
    for line in result.stdout.strip().split("\n")[2:]:
        parts = line.split()
        if len(parts) >= 2 and parts[1] != '-':
            gpu_id = int(parts[0])
            if int(parts[1]) == pid:
                return gpu_id
    return None


def get_progress(log_path):
    """Get last time_ns from a production log."""
    try:
        with open(log_path) as f:
            lines = f.readlines()
            for line in reversed(lines):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    return float(parts[1])
    except Exception:
        pass
    return 0.0


def find_idle_gpus():
    """Find GPUs with no active MD processes."""
    md_pids = get_running_md_pids()
    gpu_usage = {0: False, 1: False, 2: False, 3: False}

    for pid in md_pids:
        gpu = get_pid_gpu(pid)
        if gpu is not None:
            gpu_usage[gpu] = True

    return [gpu for gpu, used in gpu_usage.items() if not used]


def launch_rep(system, rep, gpu_id):
    outdir = BASE / "data/md_runs" / system / f"rep{rep}"
    prmtop = outdir / f"{system}.prmtop"
    pdb = outdir / f"{system}_minimized.pdb"
    name = f"{system}_rep{rep}"

    # Deterministic but different seeds per rep/system
    seed = 20250502 + rep * 1000 + (10000 if "4mut_rev" in system else 0)

    cmd = [
        PYTHON, str(RUN_SCRIPT),
        "--prmtop", str(prmtop),
        "--pdb", str(pdb),
        "--name", name,
        "--outdir", str(outdir),
        "--prod-ns", "200",
        "--platform", "CUDA",
        "--seed", str(seed),
    ]

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    print(f"[{time.strftime('%H:%M:%S')}] Launching {name} on GPU {gpu_id} (seed={seed})")

    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, env=env)
    return proc.pid


def main():
    print(f"[{time.strftime('%H:%M:%S')}] Hgal NEW auto-launcher started.")
    print(f"  Queue: {len(QUEUE)} replicas")
    print(f"  Checking every 60 seconds...")

    while QUEUE:
        idle_gpus = find_idle_gpus()

        for gpu_id in sorted(idle_gpus):
            if not QUEUE:
                break
            system, rep = QUEUE.pop(0)
            launch_rep(system, rep, gpu_id)
            time.sleep(5)

        # Progress report
        s305e_prog = {}
        for rep in [1, 2, 3]:
            log = BASE / f"data/md_runs/Hsap_WT_S305E/rep{rep}/Hsap_WT_S305E_rep{rep}_prod.log"
            s305e_prog[rep] = get_progress(log)

        hgal_wt1 = get_progress(BASE / "data/md_runs/Hgal_WT/rep1/Hgal_WT_rep1_prod.log")

        print(f"[{time.strftime('%H:%M:%S')}] S305E={s305e_prog} | "
              f"Hgal_WT_rep1={hgal_wt1:.1f}ns | Queue={len(QUEUE)} | Idle GPUs={idle_gpus}")

        time.sleep(60)

    print(f"[{time.strftime('%H:%M:%S')}] All replicas launched!")


if __name__ == "__main__":
    main()
