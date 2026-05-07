#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate boltz
export LD_LIBRARY_PATH="/home/scroll/miniforge3/envs/boltz/lib/python3.11/site-packages/nvidia/cu13/lib:$LD_LIBRARY_PATH"
cd /home/scroll/personal/naked-mole-rat-cgas-trim41-simulation
CUDA_VISIBLE_DEVICES=1 boltz predict data/boltz_cgas_dna_trim41.yaml --out_dir data/boltz_cgas_dna_trim41 --diffusion_samples 5 --use_msa_server 2>&1 | tee data/boltz_cgas_dna_trim41/boltz.log
