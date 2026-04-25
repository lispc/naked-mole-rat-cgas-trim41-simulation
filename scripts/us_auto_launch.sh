#!/bin/bash
QUEUE="data/analysis/final/us_window_queue.txt"
OUTDIR="data/analysis/final/us_Hsap_WT_Lys334"
SCRIPT="scripts/run_us_simple.py"
PRMTOP="data/md_runs/Hsap_WT/Hsap_WT.prmtop"
RST7="data/analysis/final/us_start_Hsap_WT.rst7"
LOG="$OUTDIR/scheduler.log"

mkdir -p "$OUTDIR"
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md

echo "[$(date)] Auto-launcher started" >> "$LOG"

while [ -s "$QUEUE" ]; do
  # Check each GPU (0,1,2,3)
  for gpu in 0 1 2 3; do
    # Skip GPU 2,3 if Hgal MD is still running
    if [ "$gpu" -eq 2 ] || [ "$gpu" -eq 3 ]; then
      # Check if Hgal MD is running on this GPU
      md_pid=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader -i $gpu 2>/dev/null | grep "python" | awk -F', ' '{print $1}')
      if [ -n "$md_pid" ]; then
        # Check if it's run_md.py (MD) or run_us_simple.py (US)
        cmdline=$(ps -p "$md_pid" -o args= 2>/dev/null | head -c 50)
        if echo "$cmdline" | grep -q "run_md"; then
          echo "[$(date)] GPU $gpu busy with MD (pid $md_pid), skipping" >> "$LOG"
          continue
        fi
      fi
    fi
    
    # Check if GPU is free (<10% utilization and no python US process)
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader -i $gpu 2>/dev/null | tr -d ' %')
    if [ -z "$util" ]; then util=100; fi
    
    us_on_gpu=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader -i $gpu 2>/dev/null | grep -c "python")
    
    if [ "$util" -lt 15 ] && [ "$us_on_gpu" -eq 0 ]; then
      # GPU is free, launch next window
      next=$(head -1 "$QUEUE")
      if [ -z "$next" ]; then break; fi
      
      center=$(echo "$next" | awk '{print $1}')
      k=$(echo "$next" | awk '{print $2}')
      name=$(printf "us_w%02dA" ${center%.*})
      
      # Check if already done
      if [ -f "$OUTDIR/${name}_cv.dat" ] && [ $(wc -l < "$OUTDIR/${name}_cv.dat" 2>/dev/null) -gt 20 ]; then
        echo "[$(date)] $name already done, removing from queue" >> "$LOG"
        sed -i '1d' "$QUEUE"
        continue
      fi
      
      echo "[$(date)] Launching $name (center=${center}A, k=$k) on GPU $gpu" >> "$LOG"
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
      
      sed -i '1d' "$QUEUE"
      sleep 10
    fi
  done
  
  # Also check if Hgal MD has completed
  for log in data/md_runs/Hgal_4mut_rev/rep1/Hgal_4mut_rev_rep1_prod.log data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_prod.log; do
    if [ -f "$log" ]; then
      last_ns=$(tail -1 "$log" 2>/dev/null | awk '{print $2}')
      if [ -n "$last_ns" ] && [ "${last_ns%.*}" -ge 200 ] 2>/dev/null; then
        echo "[$(date)] MD detected at ${last_ns}ns in $log" >> "$LOG"
      fi
    fi
  done
  
  sleep 30
done

echo "[$(date)] Queue empty. Waiting for all US to finish..." >> "$LOG"
while pgrep -f "run_us_simple.py" > /dev/null; do
  sleep 30
done
echo "[$(date)] All US windows completed!" >> "$LOG"

# Final analysis trigger
echo "[$(date)] Ready for final systematic analysis." >> "$LOG"
