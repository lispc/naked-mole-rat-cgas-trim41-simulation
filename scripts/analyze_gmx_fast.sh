#!/bin/bash
# Fast GROMACS analysis using native gmx tools
set -e

PROJECT_ROOT="/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation"
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx

SYSTEM=$1
REP=$2
DIR="${PROJECT_ROOT}/data/md_runs_gmx/${SYSTEM}/rep${REP}"
NAME="${SYSTEM}_rep${REP}"
OUT="${PROJECT_ROOT}/data/analysis/hsap_batch_gmx/${NAME}"
mkdir -p "$OUT"

cd "$DIR"

echo "=== Analyzing ${NAME} ==="

# 1. Create index file with cGAS, TRIM41, and CA subgroups
echo "ri 1-218" > /tmp/make_ndx_${NAME}.in
echo "name 19 cGAS" >> /tmp/make_ndx_${NAME}.in
echo "ri 219-541" >> /tmp/make_ndx_${NAME}.in
echo "name 20 TRIM41" >> /tmp/make_ndx_${NAME}.in
echo "19 & a CA" >> /tmp/make_ndx_${NAME}.in
echo "name 21 cGAS_CA" >> /tmp/make_ndx_${NAME}.in
echo "20 & a CA" >> /tmp/make_ndx_${NAME}.in
echo "name 22 TRIM41_CA" >> /tmp/make_ndx_${NAME}.in
echo "q" >> /tmp/make_ndx_${NAME}.in

gmx make_ndx -f prod.tpr -o ${NAME}.ndx < /tmp/make_ndx_${NAME}.in 2>/dev/null || true

# 2. RMSD (protein CA vs first frame)
gmx rms -s prod.tpr -f prod.xtc -n ${NAME}.ndx -o ${OUT}/rmsd.xvg -xvg none 2>/dev/null << 'SEL'
21
21
SEL

# 3. Rg (protein)
gmx gyrate -s prod.tpr -f prod.xtc -n ${NAME}.ndx -o ${OUT}/rg.xvg -xvg none 2>/dev/null << 'SEL'
1
SEL

# 4. COM distance (cGAS vs TRIM41)
gmx distance -s prod.tpr -f prod.xtc -n ${NAME}.ndx -select 'com of group "cGAS" plus com of group "TRIM41"' -o ${OUT}/com_dist.xvg -xvg none 2>/dev/null || echo "com_dist failed"

# 5. Min CA-CA distance
gmx mindist -s prod.tpr -f prod.xtc -n ${NAME}.ndx -od ${OUT}/mindist.xvg -xvg none -group 2>/dev/null << 'SEL'
21
22
SEL

echo "Done: ${OUT}"
