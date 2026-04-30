#!/bin/bash
# Monitor GROMACS replicas and auto-run analysis when all complete

PROJECT_ROOT="/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation"
LOGFILE="${PROJECT_ROOT}/data/md_runs_gmx/auto_analysis.log"

echo "[$(date)] Starting GROMACS monitor..." >> "$LOGFILE"

while true; do
    all_done=true
    for s in Hsap_WT Hsap_4mut; do
        for r in 1 2; do
            dir="${PROJECT_ROOT}/data/md_runs_gmx/${s}/rep${r}"
            log="${dir}/prod.log"
            if [ -f "$log" ]; then
                last_step=$(awk '/^\s+[0-9]+\s+[0-9]+\./{last=$1} END{print last}' "$log" 2>/dev/null)
                if [ -z "$last_step" ] || [ "$last_step" -lt 100000000 ]; then
                    all_done=false
                    echo "[$(date)] ${s} rep${r}: ${last_step:-unknown}/100000000" >> "$LOGFILE"
                else
                    echo "[$(date)] ${s} rep${r}: DONE (${last_step})" >> "$LOGFILE"
                fi
            else
                all_done=false
                echo "[$(date)] ${s} rep${r}: no log yet" >> "$LOGFILE"
            fi
        done
    done
    
    if [ "$all_done" = true ]; then
        echo "[$(date)] All replicas complete! Starting analysis..." >> "$LOGFILE"
        
        source /home/scroll/miniforge3/etc/profile.d/conda.sh
        conda activate cgas-md
        cd "$PROJECT_ROOT"
        python scripts/batch_analyze_hsap_gmx.py >> "$LOGFILE" 2>&1
        echo "[$(date)] Analysis complete. Output: data/analysis/hsap_batch_gmx/" >> "$LOGFILE"
        
        # Also create a completion marker
        touch "${PROJECT_ROOT}/data/md_runs_gmx/.analysis_complete"
        exit 0
    fi
    
    # Check every 10 minutes
    sleep 600
done
