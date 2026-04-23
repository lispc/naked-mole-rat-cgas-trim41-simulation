#!/bin/bash
# MD monitoring script with physics/business logic checks

LOGFILE="data/md_runs/Hgal_domain/monitor.log"
OUTDIR="data/md_runs/Hgal_domain"
TARGET_TEMP=300
TARGET_NS=200

echo "========================================" >> "$LOGFILE"
echo "MD Status Check: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOGFILE"
echo "========================================" >> "$LOGFILE"

# --- Hardware Check ---
echo "--- Hardware ---" >> "$LOGFILE"
nvidia-smi --query-gpu=index,temperature.gpu,power.draw,utilization.gpu --format=csv,noheader >> "$LOGFILE"

# --- Process Check ---
echo "--- Processes ---" >> "$LOGFILE"
for i in 1 2 3; do
    pid=$(pgrep -f "Hgal_domain_rep${i}" | head -1)
    if [ -n "$pid" ]; then
        cpu=$(ps -p "$pid" -o %cpu --noheaders 2>/dev/null | tr -d ' ')
        echo "rep${i}: PID=$pid CPU=${cpu}%" >> "$LOGFILE"
    else
        echo "rep${i}: MISSING PROCESS!" >> "$LOGFILE"
    fi
done

# --- Physics/Business Logic Checks for each replica ---
for rep in rep1 rep2 rep3; do
    echo "" >> "$LOGFILE"
    echo "--- $rep Physics Check ---" >> "$LOGFILE"
    
    # 1. Heating check: temperature ramp 0->300K
    heat_log="$OUTDIR/$rep/Hgal_domain_${rep}_heating.log"
    if [ -s "$heat_log" ]; then
        first_temp=$(tail -n +2 "$heat_log" | head -1 | awk '{print $4}')
        last_temp=$(tail -1 "$heat_log" | awk '{print $4}')
        echo "Heating: start=${first_temp}K end=${last_temp}K" >> "$LOGFILE"
        if (( $(echo "$last_temp < 180 || $last_temp > 320" | bc -l) )); then
            echo "  WARNING: Final heating temp abnormal (expected 200-300K for short heating)!" >> "$LOGFILE"
        else
            echo "  OK: Heating temp reached ${last_temp}K (short 100ps ramp for 116k atoms)" >> "$LOGFILE"
        fi
    fi
    
    # 2. NPT check: volume convergence
    npt_log="$OUTDIR/$rep/Hgal_domain_${rep}_npt.log"
    if [ -s "$npt_log" ]; then
        first_vol=$(tail -n +2 "$npt_log" | head -1 | awk '{print $4}')
        last_vol=$(tail -1 "$npt_log" | awk '{print $4}')
        vol_change=$(echo "scale=2; ($last_vol - $first_vol) / $first_vol * 100" | bc)
        echo "NPT: start_vol=${first_vol} end_vol=${last_vol} change=${vol_change}%" >> "$LOGFILE"
    fi
    
    # 3. Production check: energy stability and drift
    prod_log="$OUTDIR/$rep/Hgal_domain_${rep}_prod.log"
    if [ -s "$prod_log" ]; then
        n_lines=$(wc -l < "$prod_log")
        if [ "$n_lines" -gt 1 ]; then
            # Get first and last energy values
            first_energy=$(tail -n +2 "$prod_log" | head -1 | awk '{print $3}')
            last_energy=$(tail -1 "$prod_log" | awk '{print $3}')
            last_ns=$(tail -1 "$prod_log" | awk '{print $2}')
            
            # Check for NaN or extremely high energy
            if echo "$last_energy" | grep -qi "nan\|inf"; then
                echo "  CRITICAL: Energy is NaN/Inf!" >> "$LOGFILE"
                continue
            fi
            
            # Check energy drift (kJ/mol per ns)
            if [ -n "$first_energy" ] && [ -n "$last_energy" ] && [ -n "$last_ns" ]; then
                energy_drift=$(echo "scale=2; ($last_energy - $first_energy) / $last_ns" | bc)
                abs_drift=$(echo "$energy_drift" | sed 's/-//')
                echo "Production: ${last_ns}ns | Energy: ${last_energy} kJ/mol | Drift: ${energy_drift} kJ/mol/ns" >> "$LOGFILE"
                
                # Warnings based on drift magnitude
                if (( $(echo "$abs_drift > 5000" | bc -l) )); then
                    echo "  WARNING: Severe energy drift (>5000 kJ/mol/ns)!" >> "$LOGFILE"
                elif (( $(echo "$abs_drift > 1000" | bc -l) )); then
                    echo "  CAUTION: Moderate energy drift (>1000 kJ/mol/ns)" >> "$LOGFILE"
                else
                    echo "  OK: Energy stable" >> "$LOGFILE"
                fi
                
                # Check if energy is physically reasonable for this system
                # Expected range for solvated protein: -1.5M to -1.0M kJ/mol
                if (( $(echo "$last_energy < -2000000 || $last_energy > -500000" | bc -l) )); then
                    echo "  WARNING: Energy outside expected range (-2M to -0.5M)!" >> "$LOGFILE"
                fi
            fi
            
            # 4. DCD trajectory file check
            dcd_file="$OUTDIR/$rep/Hgal_domain_${rep}_prod.dcd"
            if [ -f "$dcd_file" ]; then
                dcd_size=$(stat -c %s "$dcd_file" 2>/dev/null)
                expected_size_per_frame=1400652  # ~1.4MB per frame for 116k atoms
                n_frames=$((dcd_size / expected_size_per_frame))
                echo "  DCD: ${n_frames} frames, $(echo "scale=1; $dcd_size/1048576" | bc) MB" >> "$LOGFILE"
            else
                echo "  WARNING: DCD file missing!" >> "$LOGFILE"
            fi
        else
            echo "Production: log started, no data yet" >> "$LOGFILE"
        fi
    else
        echo "Production: log empty" >> "$LOGFILE"
    fi
done

# --- Speed & ETA Calculation ---
echo "" >> "$LOGFILE"
echo "--- ETA Estimate ---" >> "$LOGFILE"
rep1_log="$OUTDIR/rep1/Hgal_domain_rep1_prod.log"
if [ -s "$rep1_log" ]; then
    state_file="$OUTDIR/.monitor_state"
    current_time=$(date +%s)
    last_ns=$(tail -1 "$rep1_log" | awk '{print $2}')
    
    if [ -f "$state_file" ]; then
        read prev_time prev_ns < "$state_file"
        time_diff=$((current_time - prev_time))
        ns_diff=$(echo "$last_ns - $prev_ns" | bc -l)
        if (( time_diff > 0 )) && (( $(echo "$ns_diff > 0" | bc -l) )); then
            ns_per_day=$(echo "scale=1; $ns_diff / $time_diff * 86400" | bc -l)
            remaining=$(echo "200 - $last_ns" | bc -l)
            hours_left=$(echo "scale=1; $remaining / $ns_per_day * 24" | bc -l)
            eta=$(date -d "+${hours_left} hours" '+%m-%d %H:%M' 2>/dev/null)
            echo "Speed: ${ns_per_day} ns/day" >> "$LOGFILE"
            echo "Rep1 remaining: ${remaining}ns | ETA: $eta (~${hours_left}h)" >> "$LOGFILE"
            pct=$(echo "scale=1; $last_ns / 2" | bc -l)
            echo "Overall progress: ${pct}%" >> "$LOGFILE"
        fi
    fi
    echo "$current_time $last_ns" > "$state_file"
fi

echo "" >> "$LOGFILE"
