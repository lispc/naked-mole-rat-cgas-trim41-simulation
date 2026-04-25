#!/bin/bash
# Umbrella Sampling scheduler: auto-launch next window when GPU is free

PRMTOP="data/md_runs/Hsap_WT/Hsap_WT.prmtop"
RST7="data/analysis/final/us_start_Hsap_WT.rst7"
OUTDIR="data/analysis/final/us_Hsap_WT_Lys334"
SCRIPT="scripts/run_us_simple.py"

declare -a WINDOWS=(
  "4.0:1000:0"
  "5.0:1000:0"
  "6.0:1000:0"
  "7.0:1000:0"
  "9.0:1000:0"
  "10.0:1000:0"
  "11.0:1000:0"
  "13.0:500:1"
  "14.0:500:1"
  "15.0:500:1"
  "16.0:500:1"
  "17.0:500:1"
  "18.0:500:1"
  "19.0:500:1"
  "20.0:500:1"
)

# Also schedule Hgal 4mut_rev US after MD completes
# For now, just track Hsap windows

LOGFILE="$OUTDIR/scheduler.log"
mkdir -p "$OUTDIR"

source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md

echo "[$(date)] Scheduler started. Total remaining windows: ${#WINDOWS[@]}" >> "$LOGFILE"

idx=0
while [ $idx -lt ${#WINDOWS[@]} ]; do
  IFS=':' read -r center k gpu <<< "${WINDOWS[$idx]}"
  name=$(printf "us_w%02dA" ${center%.*})
  
  # Check if already running or completed
  if pgrep -f "$name" > /dev/null; then
    echo "[$(date)] $name already running, skipping" >> "$LOGFILE"
    idx=$((idx+1))
    continue
  fi
  if [ -f "$OUTDIR/${name}_cv.dat" ] && [ $(wc -l < "$OUTDIR/${name}_cv.dat" 2>/dev/null) -gt 10 ]; then
    echo "[$(date)] $name already completed, skipping" >> "$LOGFILE"
    idx=$((idx+1))
    continue
  fi
  
  # Check if target GPU is free (no US process on it)
  gpu_busy=$(nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv | grep -c ", python,")
  # Simpler: check if any run_us_simple on this GPU
  # We assign windows round-robin to GPU 0/1 for now
  
  # Check specific GPU utilization
  util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader -i $gpu | tr -d ' %')
  if [ "$util" -gt 10 ]; then
    echo "[$(date)] GPU $gpu busy (${util}%), waiting for $name" >> "$LOGFILE"
    sleep 60
    continue
  fi
  
  echo "[$(date)] Launching $name (center=${center}A, k=$k) on GPU $gpu" >> "$LOGFILE"
  nohup python "$SCRIPT" \
    --prmtop "$PRMTOP" \
    --rst7 "$RST7" \
    --center "$center" \
    --k "$k" \
    --name "$name" \
    --outdir "$OUTDIR" \
    --gpu "$gpu" \
    --em-steps 1000 \
    --prod-ns 20.0 >> "$OUTDIR/${name}.log" 2>&1 &
  
  sleep 30
  idx=$((idx+1))
done

echo "[$(date)] All Hsap windows scheduled. Waiting for completion..." >> "$LOGFILE"

# Wait for all US processes to finish
while pgrep -f "run_us_simple.py" > /dev/null; do
  sleep 60
done

echo "[$(date)] All US windows completed!" >> "$LOGFILE"
