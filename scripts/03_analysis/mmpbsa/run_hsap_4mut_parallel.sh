#!/bin/bash
# Parallel MM-GBSA for Hsap_4mut rep2 + rep3
# Rep1 already completed; rep2/rep3 run in parallel.

set -e
BASE="/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation"
AMBERHOME="/home/scroll/miniforge3/envs/cgas-md"
MMPBSA_PY="$AMBERHOME/bin/MMPBSA.py"
OUTDIR="$BASE/data/analysis/mmpbsa"

CP="$BASE/data/md_runs/Hsap_4mut/Hsap_4mut_protein.prmtop"
RP="$BASE/data/md_runs/Hsap_4mut/Hsap_4mut_receptor.prmtop"
LP="$BASE/data/md_runs/Hsap_4mut/Hsap_4mut_ligand.prmtop"

export AMBERHOME
export OMP_NUM_THREADS=4

run_rep() {
    local name=$1
    local nc=$2
    local workdir="$OUTDIR/${name}_delta_work"
    mkdir -p "$workdir"

    cat > "$workdir/mmpbsa.in" <<'EOF'
&general
startframe=1, endframe=99999, interval=50,
keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=1,
print_res="all",
/
EOF

    local results="$OUTDIR/${name}_delta_results.dat"
    local decomp="$OUTDIR/${name}_delta_decomp.dat"
    local log="$workdir/mmpbsa.log"

    echo "[$name] Starting..."
    if "$MMPBSA_PY" -i "$workdir/mmpbsa.in" -o "$results" -do "$decomp" \
        -cp "$CP" -rp "$RP" -lp "$LP" -y "$nc" > "$log" 2>&1; then
        echo "[$name] ✅ Done"
    else
        echo "[$name] ❌ Failed (exit $?)"
    fi
}

# Launch rep2 and rep3 in parallel
run_rep "Hsap_4mut_rep2" "$OUTDIR/Hsap_4mut_rep2_prot.nc" &
run_rep "Hsap_4mut_rep3" "$OUTDIR/Hsap_4mut_rep3_prot.nc" &

wait
echo "All parallel jobs finished."
