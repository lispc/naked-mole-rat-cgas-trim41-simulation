#!/bin/bash
# Parallel MM-GBSA for S305E (3 reps) - v3 using protein-only .nc trajectories

set -e
BASE="/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation"
AMBERHOME="/home/scroll/miniforge3/envs/cgas-md"
MMPBSA_PY="$AMBERHOME/bin/MMPBSA.py"
OUTDIR="$BASE/data/analysis/mmpbsa"
mkdir -p "$OUTDIR"

CP="$BASE/data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_complex.prmtop"
RP="$BASE/data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_receptor_new.prmtop"
LP="$BASE/data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_ligand_new.prmtop"

export AMBERHOME
export OMP_NUM_THREADS=4

run_rep() {
    local rep=$1
    local name="S305E_rep${rep}"
    local nc="$BASE/data/md_runs/Hsap_WT_S305E/rep${rep}/protein.nc"
    local workdir="$OUTDIR/${name}_work_v3"
    mkdir -p "$workdir"

    cat > "$workdir/mmpbsa.in" <<'EOF'
&general
startframe=1, endframe=99999, interval=50,
keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
EOF

    local results="$OUTDIR/${name}_results_v3.dat"
    local log="$workdir/mmpbsa.log"

    echo "[$name] Starting MM-GBSA..."
    cd "$workdir"
    if "$MMPBSA_PY" -i "$workdir/mmpbsa.in" -o "$results" \
        -cp "$CP" -rp "$RP" -lp "$LP" -y "$nc" > "$log" 2>&1; then
        echo "[$name] Done"
    else
        echo "[$name] Failed (exit $?)"
    fi
}

for rep in 1 2 3; do
    ( run_rep "$rep" ) &
done
wait

echo "All S305E MM-GBSA reps complete."
