#!/bin/bash
# Run one replica through EM -> NVT -> NPT -> Production
set -e

PROJECT_ROOT="/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation"
SYSTEM=$1
REP=$2
GPU_ID=$3

BASE_DIR="${PROJECT_ROOT}/data/md_runs_gmx/${SYSTEM}/rep${REP}"
MDP_DIR="${PROJECT_ROOT}/data/md_runs_gmx/mdp"
NAME="${SYSTEM}_rep${REP}"

mkdir -p "${BASE_DIR}"
cd "${BASE_DIR}"

echo "========================================"
echo "Starting ${NAME} on GPU ${GPU_ID}"
echo "Working dir: $(pwd)"
echo "========================================"

# ---- Energy Minimization ----
echo "[${NAME}] Step 1: Energy Minimization"
if [ ! -f "em.gro" ]; then
    gmx grompp -f "${MDP_DIR}/em.mdp" -c "${NAME}.gro" -p "${NAME}.top" -o em.tpr -maxwarn 10 2>&1 | tee em_grompp.log
    gmx mdrun -deffnm em -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu 2>&1 | tee em.log
    echo "[${NAME}] EM completed"
else
    echo "[${NAME}] EM already done, skipping"
fi

# ---- NVT Equilibration ----
echo "[${NAME}] Step 2: NVT Equilibration"
if [ ! -f "nvt.gro" ]; then
    gmx grompp -f "${MDP_DIR}/nvt.mdp" -c em.gro -p "${NAME}.top" -o nvt.tpr -maxwarn 10 2>&1 | tee nvt_grompp.log
    gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu -pme gpu 2>&1 | tee nvt.log
    echo "[${NAME}] NVT completed"
else
    echo "[${NAME}] NVT already done, skipping"
fi

# ---- NPT Equilibration ----
echo "[${NAME}] Step 3: NPT Equilibration"
if [ ! -f "npt.gro" ]; then
    gmx grompp -f "${MDP_DIR}/npt.mdp" -c nvt.gro -p "${NAME}.top" -t nvt.cpt -o npt.tpr -maxwarn 10 2>&1 | tee npt_grompp.log
    gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu -pme gpu 2>&1 | tee npt.log
    echo "[${NAME}] NPT completed"
else
    echo "[${NAME}] NPT already done, skipping"
fi

# ---- Production Run ----
echo "[${NAME}] Step 4: Production MD 200ns"
if [ ! -f "prod.gro" ]; then
    gmx grompp -f "${MDP_DIR}/prod.mdp" -c npt.gro -p "${NAME}.top" -t npt.cpt -o prod.tpr -maxwarn 10 2>&1 | tee prod_grompp.log
    gmx mdrun -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu -pme gpu 2>&1 | tee prod.log
    echo "[${NAME}] Production completed"
else
    echo "[${NAME}] Production already done, skipping"
fi

echo "[${NAME}] All steps completed!"
