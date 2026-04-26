#!/bin/bash
# Launch Hsap_4mut US windows across 4 GPUs, 4 at a time

source /home/scroll/miniforge3/etc/profile.d/conda.sh && conda activate cgas-md

PRMTOP="data/md_runs/Hsap_4mut/Hsap_4mut.prmtop"
RST7="data/analysis/final/us_start_Hsap_4mut.rst7"
OUTDIR="data/analysis/final/us_Hsap_4mut_Lys334"
mkdir -p "$OUTDIR"

# Windows: name:center:k
declare -a WINDOWS=(
    "w04A:4.0:1000"
    "w05A:5.0:1000"
    "w06A:6.0:1000"
    "w07A:7.0:1000"
    "w08A:8.0:1000"
    "w09A:9.0:1000"
    "w10A:10.0:1000"
    "w11A:11.0:1000"
    "w12A:12.0:1000"
    "w13A:13.0:500"
    "w14A:14.0:500"
    "w15A:15.0:500"
    "w16A:16.0:500"
    "w17A:17.0:500"
    "w18A:18.0:500"
    "w19A:19.0:500"
    "w20A:20.0:500"
)

gpu=0
batch=0
for spec in "${WINDOWS[@]}"; do
    IFS=':' read -r name center k <<< "$spec"
    
    LOG="$OUTDIR/us_${name}.log"
    
    echo "[$(date)] Launching us_${name} (center=${center}A, k=${k}) on GPU ${gpu}"
    
    python scripts/run_us_simple.py \
        --prmtop "$PRMTOP" \
        --rst7 "$RST7" \
        --center "$center" \
        --k "$k" \
        --name "us_${name}" \
        --outdir "$OUTDIR" \
        --gpu "$gpu" \
        --em-steps 1000 \
        --prod-ns 20.0 \
        > "$LOG" 2>&1 &
    
    gpu=$((gpu + 1))
    
    # Wait for current batch of 4 to finish before next batch
    if [ "$gpu" -ge 4 ]; then
        gpu=0
        batch=$((batch + 1))
        echo "[$(date)] Waiting for batch ${batch} to complete..."
        wait
        echo "[$(date)] Batch ${batch} done."
    fi
done

# Wait for final batch
echo "[$(date)] Waiting for final batch..."
wait
echo "[$(date)] All Hsap_4mut US windows completed!"
