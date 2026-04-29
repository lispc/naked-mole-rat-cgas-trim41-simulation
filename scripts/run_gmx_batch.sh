#!/bin/bash
# Launch 4 GROMACS replicas in parallel (one per GPU)
set -e

source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx

echo "Launching GROMACS batch..."
echo "WT rep1  -> GPU 0"
echo "WT rep2  -> GPU 1"
echo "4mut rep1 -> GPU 2"
echo "4mut rep2 -> GPU 3"

nohup bash scripts/run_gmx_replica.sh Hsap_WT 1 0 > data/md_runs_gmx/Hsap_WT/rep1/run.log 2>&1 &
nohup bash scripts/run_gmx_replica.sh Hsap_WT 2 1 > data/md_runs_gmx/Hsap_WT/rep2/run.log 2>&1 &
nohup bash scripts/run_gmx_replica.sh Hsap_4mut 1 2 > data/md_runs_gmx/Hsap_4mut/rep1/run.log 2>&1 &
nohup bash scripts/run_gmx_replica.sh Hsap_4mut 2 3 > data/md_runs_gmx/Hsap_4mut/rep2/run.log 2>&1 &

echo "All replicas launched in background."
echo "Monitor with: tail -f data/md_runs_gmx/*/rep*/run.log"
