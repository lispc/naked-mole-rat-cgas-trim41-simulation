#!/bin/bash
# Prepare GROMACS 2026 native topology (amber19sb.ff + opc water)
# This avoids the parmed CMAP conversion bug entirely.
#
# Usage:
#   ./prepare_gromacs_2026_native.sh Hsap_WT 1
#
# Steps:
#   1. Extract protein coordinates from OpenMM NPT final frame
#   2. gmx pdb2gmx -ff amber19sb -water opc (correct CMAP, correct topology)
#   3. gmx solvate + gmx genion (rebuild water box and ions)
#   4. EM + NVT + NPT equilibration
#   5. Production run with fixed prod.mdp

set -e

PROJECT_ROOT="/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation"
SYSTEM=$1
REP=$2
GPU_ID=${3:-0}

# Source conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate gmx

BASE_DIR="${PROJECT_ROOT}/data/md_runs_gmx2026/${SYSTEM}/rep${REP}"
MDP_DIR="${PROJECT_ROOT}/data/md_runs_gmx/mdp"
NAME="${SYSTEM}_rep${REP}"

mkdir -p "${BASE_DIR}"
cd "${BASE_DIR}"

echo "========================================"
echo "Preparing ${NAME} with GROMACS 2026 native ff19SB"
echo "Working dir: $(pwd)"
echo "========================================"

# ---- Step 0: Extract protein from OpenMM NPT last frame ----
# Use MDAnalysis to extract last frame from DCD, save as PDB
if [ ! -f "${NAME}_protein.pdb" ]; then
    echo "[${NAME}] Extracting protein from OpenMM NPT trajectory..."
    PRMTOP="${PROJECT_ROOT}/data/md_runs/${SYSTEM}/${SYSTEM}.prmtop"
    DCD="${PROJECT_ROOT}/data/md_runs/${SYSTEM}/rep${REP}/${SYSTEM}_rep${REP}_npt.dcd"
    OUT_PDB="${NAME}_protein.pdb"
    conda run -n cgas-md python3 -c "
import MDAnalysis as mda
u = mda.Universe('${PRMTOP}', '${DCD}')
u.trajectory[-1]
cgas = u.select_atoms('protein and resid 1-218')
trim41 = u.select_atoms('protein and resid 219-541')
with open('${OUT_PDB}', 'w') as f:
    for atom in cgas:
        f.write(f'ATOM  {atom.id:5d} {atom.name:4s} {atom.resname:3s} A{atom.resid:4d}    '
                f'{atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}'
                f'{1.0:6.2f}{0.0:6.2f}          {atom.element:>2s}\\n')
    f.write('TER\\n')
    for atom in trim41:
        f.write(f'ATOM  {atom.id:5d} {atom.name:4s} {atom.resname:3s} B{atom.resid:4d}    '
                f'{atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}'
                f'{1.0:6.2f}{0.0:6.2f}          {atom.element:>2s}\\n')
    f.write('END\\n')
print(f'Saved {len(cgas) + len(trim41)} protein atoms to ${OUT_PDB}')
" 2>&1 | tee extract_protein.log
fi

# ---- Step 1: gmx pdb2gmx (native ff19SB + OPC) ----
echo "[${NAME}] Step 1: gmx pdb2gmx (amber19sb.ff + opc)"
if [ ! -f "${NAME}_processed.gro" ]; then
    gmx pdb2gmx -f "${NAME}_protein.pdb" \
                  -o "${NAME}_processed.gro" \
                  -p "${NAME}_native.top" \
                  -ff amber19sb \
                  -water opc \
                  -ignh \
                  <<EOF 2>&1 | tee pdb2gmx.log
1
EOF
    echo "[${NAME}] pdb2gmx completed"
fi

# ---- Step 2: Solvate ----
echo "[${NAME}] Step 2: Solvating with OPC water"
if [ ! -f "${NAME}_solvated.gro" ]; then
    # Use dodecahedron box with 1.2nm buffer (similar to OpenMM)
    gmx editconf -f "${NAME}_processed.gro" \
                 -o "${NAME}_boxed.gro" \
                 -d 1.2 -bt dodecahedron
    
    gmx solvate -cp "${NAME}_boxed.gro" \
                -cs tip4p.gro \
                -o "${NAME}_solvated.gro" \
                -p "${NAME}_native.top"
    echo "[${NAME}] Solvation completed"
fi

# ---- Step 3: Add ions (neutralize + 150mM NaCl) ----
echo "[${NAME}] Step 3: Adding ions"
if [ ! -f "${NAME}_ionized.gro" ]; then
    gmx grompp -f "${MDP_DIR}/ions.mdp" \
               -c "${NAME}_solvated.gro" \
               -p "${NAME}_native.top" \
               -o ions.tpr -maxwarn 10 2>&1 | tee ions_grompp.log
    
    # Replace solvent with ions: neutralize, then add Na+ and Cl-
    gmx genion -s ions.tpr \
               -o "${NAME}_ionized.gro" \
               -p "${NAME}_native.top" \
               -pname NA -nname CL \
               -neutral \
               -conc 0.15 \
               <<EOF 2>&1 | tee genion.log
13
EOF
    echo "[${NAME}] Ionization completed"
fi

# ---- Step 4: EM ----
echo "[${NAME}] Step 4: Energy Minimization"
if [ ! -f "em.gro" ]; then
    gmx grompp -f "${MDP_DIR}/em.mdp" \
               -c "${NAME}_ionized.gro" \
               -p "${NAME}_native.top" \
               -o em.tpr -maxwarn 10 2>&1 | tee em_grompp.log
    gmx mdrun -deffnm em -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu 2>&1 | tee em.log
    echo "[${NAME}] EM completed"
fi

# ---- Step 5: NVT Equilibration ----
echo "[${NAME}] Step 5: NVT Equilibration"
if [ ! -f "nvt.gro" ]; then
    gmx grompp -f "${MDP_DIR}/nvt.mdp" \
               -c em.gro \
               -p "${NAME}_native.top" \
               -o nvt.tpr -maxwarn 10 2>&1 | tee nvt_grompp.log
    gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu -pme gpu 2>&1 | tee nvt.log
    echo "[${NAME}] NVT completed"
fi

# ---- Step 6: NPT Equilibration ----
echo "[${NAME}] Step 6: NPT Equilibration"
if [ ! -f "npt.gro" ]; then
    gmx grompp -f "${MDP_DIR}/npt.mdp" \
               -c nvt.gro \
               -p "${NAME}_native.top" \
               -t nvt.cpt \
               -o npt.tpr -maxwarn 10 2>&1 | tee npt_grompp.log
    gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu -pme gpu 2>&1 | tee npt.log
    echo "[${NAME}] NPT completed"
fi

# ---- Step 7: Production Run ----
echo "[${NAME}] Step 7: Production MD 200ns"
if [ ! -f "prod.gro" ]; then
    gmx grompp -f "${MDP_DIR}/prod.mdp" \
               -c npt.gro \
               -p "${NAME}_native.top" \
               -t npt.cpt \
               -o prod.tpr -maxwarn 10 2>&1 | tee prod_grompp.log
    gmx mdrun -deffnm prod -ntmpi 1 -ntomp 8 -gpu_id ${GPU_ID} -nb gpu -pme gpu 2>&1 | tee prod.log
    echo "[${NAME}] Production completed"
fi

echo "[${NAME}] All steps completed with GROMACS 2026 native ff19SB!"
