#!/bin/bash
# Batch MM-GBSA for WT and S305-phos (3 reps each)
# Must be run from project root

set -e
AMBERHOME=/home/scroll/miniforge3/envs/cgas-md
MMPBSA="$AMBERHOME/bin/python $AMBERHOME/bin/MMPBSA.py"
OUTDIR=data/analysis/mmpbsa

mkdir -p $OUTDIR

# Function to extract protein trajectory
extract_prot() {
    local prmtop=$1
    local dcd=$2
    local outnc=$3
    local name=$(basename $dcd .dcd)
    
    if [ -f "$outnc" ]; then
        echo "Using existing: $outnc"
        return
    fi
    
    cat > /tmp/${name}_extract.in <<EOF
parm $prmtop
trajin $dcd
strip :WAT,Cl-,Na+
trajout $outnc netcdf
run
EOF
    $AMBERHOME/bin/cpptraj -i /tmp/${name}_extract.in > /dev/null 2>&1
    echo "Extracted: $outnc"
}

# Function to run MMPBSA
run_gbsa() {
    local name=$1
    local prmtop=$2
    local nc=$3
    local interval=$4
    
    cat > $OUTDIR/${name}_mmpbsa.in <<EOF
&general
  startframe=1, endframe=99999, interval=$interval,
  keep_files=0, use_sander=0,
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=1,
  print_res="all",
/
EOF

    $MMPBSA \
        -i $OUTDIR/${name}_mmpbsa.in \
        -o $OUTDIR/${name}_results.dat \
        -do $OUTDIR/${name}_decomp.dat \
        -cp $prmtop \
        -y $nc \
        > $OUTDIR/${name}_mmpbsa.log 2>&1
    
    echo "✅ Completed: $name"
}

# Create protein-only prmtops if needed
for sys in Hsap_WT Hsap_WT_S305phos; do
    prot_prmtop=data/md_runs/$sys/${sys}_protein.prmtop
    if [ ! -f "$prot_prmtop" ]; then
        echo "Creating protein prmtop for $sys..."
        $AMBERHOME/bin/python -c "
import parmed as pmd
parm = pmd.load_file('data/md_runs/$sys/${sys}.prmtop')
parm.strip(':WAT,Cl-,Na+')
parm.save('$prot_prmtop')
" 2>/dev/null
    fi
done

echo "Starting MM-GBSA batch..."

# WT (3 reps, interval=50)
for rep in 1 2 3; do
    name=Hsap_WT_rep${rep}
    dcd=data/md_runs/Hsap_WT/rep${rep}/${name}_prod.dcd
    nc=$OUTDIR/${name}_prot.nc
    prmtop=data/md_runs/Hsap_WT/Hsap_WT_protein.prmtop
    
    extract_prot data/md_runs/Hsap_WT/Hsap_WT.prmtop $dcd $nc
    run_gbsa $name $prmtop $nc 50
done

# S305-phos (3 reps, interval=50)
for rep in 1 2 3; do
    name=Hsap_WT_S305phos_rep${rep}
    dcd=data/md_runs/Hsap_WT_S305phos/rep${rep}/${name}_prod.dcd
    nc=$OUTDIR/${name}_prot.nc
    prmtop=data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_protein.prmtop
    
    extract_prot data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop $dcd $nc
    run_gbsa $name $prmtop $nc 50
done

echo "All done!"
