#!/bin/bash
# Launch 3 replicas of S305-phos OpenMM MD (from minimized structure)
set -e

PROJECT=/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation
PRMTOP=$PROJECT/data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop
PDB=$PROJECT/data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_minimized.pdb
SCRIPT=$PROJECT/scripts/run_md.py
export PATH="/home/scroll/miniforge3/bin:$PATH"

for rep in 1 2 3; do
    OUTDIR=$PROJECT/data/md_runs/Hsap_WT_S305phos/rep${rep}
    mkdir -p $OUTDIR
    
    GPU=$((rep - 1))
    SEED=$((1000 + rep))
    NAME="Hsap_WT_S305phos_rep${rep}"
    
    echo "Launching $NAME on CUDA:$GPU with seed $SEED"
    
    # Use CUDA_VISIBLE_DEVICES to isolate each replica to one GPU
    export CUDA_VISIBLE_DEVICES=$GPU
    nohup conda run -n cgas-md python $SCRIPT \
        --prmtop $PRMTOP \
        --pdb $PDB \
        --name $NAME \
        --outdir $OUTDIR \
        --prod-ns 200 \
        --platform CUDA \
        --seed $SEED \
        > $OUTDIR/${NAME}_run.log 2>&1 &
    
    echo "  PID: $! -> $OUTDIR/${NAME}_run.log"
done

echo "All 3 replicas launched."
