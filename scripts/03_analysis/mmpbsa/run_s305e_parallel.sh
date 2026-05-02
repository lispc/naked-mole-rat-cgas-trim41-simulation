#!/bin/bash
# Parallel MM-GBSA for S305E (3 reps)
#
# PREREQUISITE: Convert DCD trajectories to protein-only NetCDF first:
#   python3 -c "
#     import MDAnalysis as mda
#     for rep in [1,2,3]:
#       u = mda.Universe('data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E.prmtop',
#                        f'data/md_runs/Hsap_WT_S305E/rep{rep}/Hsap_WT_S305E_rep{rep}_prod.dcd')
#       protein = u.select_atoms('protein')
#       with mda.Writer(f'data/md_runs/Hsap_WT_S305E/rep{rep}/protein.nc', n_atoms=len(protein)) as w:
#         for ts in u.trajectory: w.write(protein)
#   "
#
# PREREQUISITE: Generate split prmtops with ante-MMPBSA.py:
#   export AMBERHOME=/home/scroll/miniforge3/envs/cgas-md
#   $AMBERHOME/bin/ante-MMPBSA.py \
#     -p data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E.prmtop \
#     -c data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_complex.prmtop \
#     -r data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_receptor_new.prmtop \
#     -l data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_ligand_new.prmtop \
#     -s ":WAT,Na+,Cl-" -n ":1-218"

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
    local workdir="$OUTDIR/${name}_work"
    mkdir -p "$workdir"

    cat > "$workdir/mmpbsa.in" <<'EOFIN'
&general
startframe=1, endframe=99999, interval=50,
keep_files=0,
/
&gb
igb=5, saltcon=0.150,
/
EOFIN

    local results="$OUTDIR/${name}_results.dat"
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
