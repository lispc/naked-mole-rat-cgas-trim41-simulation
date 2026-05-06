#!/usr/bin/env python3
"""Restart production MD from a checkpoint file.

Usage:
    python restart_production.py --prmtop Hgal_WT.prmtop \
                                 --pdb Hgal_WT_minimized.pdb \
                                 --checkpoint Hgal_WT_rep2_prod_111ns.chk \
                                 --name Hgal_WT_rep2 \
                                 --outdir data/md_runs/Hgal_WT/rep2 \
                                 --prod-ns 89 \
                                 --platform CUDA \
                                 --seed 20252502
"""
import argparse
import os
import sys
import time

import openmm as mm
from openmm import app


def get_platform(platform_name='auto'):
    """Select best OpenMM platform."""
    available = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    
    if platform_name == 'auto':
        for name in ['CUDA', 'OpenCL', 'CPU']:
            if name in available:
                platform_name = name
                break
        else:
            platform_name = 'Reference'
    
    platform = mm.Platform.getPlatformByName(platform_name)
    
    props = {}
    if platform_name == 'CUDA':
        props['CudaPrecision'] = 'mixed'
        props['CudaDeviceIndex'] = '0'
    elif platform_name == 'OpenCL':
        props['OpenCLPrecision'] = 'single'
    
    return platform, props


def run_production(prmtop, pdb, checkpoint, stage_name, outdir,
                   prod_ns=200, step_size=0.002, report_interval=50000, seed=None):
    """Restart production NVT MD from checkpoint."""
    print(f"\n[{'='*60}]")
    print(f"Production NVT MD (restart)")
    print(f"  Duration: {prod_ns} ns")
    print(f"  Step size: {step_size} ps")
    print(f"  Total steps: {int(prod_ns * 1000 / step_size):,}")
    print(f"  Report every: {report_interval} steps ({report_interval * step_size} ps)")
    print(f"  Checkpoint: {checkpoint}")
    print(f"[{'='*60}]\n")
    
    # Load topology and coordinates
    print("[Setup] Loading system...")
    prmtop_obj = app.AmberPrmtopFile(prmtop)
    pdb_obj = app.PDBFile(pdb)
    topology = prmtop_obj.topology
    
    system = prmtop_obj.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * mm.unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    
    # NVT integrator
    integrator = mm.LangevinMiddleIntegrator(
        300 * mm.unit.kelvin,
        1.0 / mm.unit.picosecond,
        step_size * mm.unit.picosecond,
    )
    if seed is not None:
        integrator.setRandomNumberSeed(seed)
    
    platform, props = get_platform('CUDA')
    context = mm.Context(system, integrator, platform, props)
    
    # Load checkpoint
    print(f"[Setup] Loading checkpoint from {checkpoint}...")
    with open(checkpoint, 'rb') as f:
        context.loadCheckpoint(f.read())
    
    state = context.getState(getPositions=True, getEnergy=True, getVelocities=True)
    start_step = state.getStepCount()
    start_time_ns = start_step * step_size / 1000
    print(f"  Resuming from step {start_step:,} ({start_time_ns:.3f} ns)")
    
    # Reporters
    dcd_file = os.path.join(outdir, f"{stage_name}_restart.dcd")
    log_file = os.path.join(outdir, f"{stage_name}_restart.log")
    
    print(f"[Setup] DCD output: {dcd_file}")
    dcd = app.DCDFile(open(dcd_file, 'wb'), topology, step_size * mm.unit.picosecond)
    
    print(f"[Setup] Log output: {log_file}")
    log = open(log_file, 'w')
    log.write(f"# Restart from {checkpoint}\n")
    log.write(f"# Starting step: {start_step:,}\n")
    log.write(f"# {'Step':>10s}  {'Time(ns)':>8s}  {'Potential(kJ/mol)':>18s}\n")
    
    # Simulation
    steps = int(prod_ns * 1000 / step_size)
    chk_interval = int(1.0 / step_size)
    
    start = time.time()
    last_report = start
    
    for step in range(0, steps, report_interval):
        integrator.step(report_interval)
        
        state = context.getState(getPositions=True, getEnergy=True)
        e = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        current_step = start_step + step + report_interval
        time_ns = current_step * step_size / 1000
        
        log.write(f"{current_step:10d}  {time_ns:8.3f}  {e:14.2f}\n")
        log.flush()
        
        dcd.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
        
        if current_step % (1000 * chk_interval) == 0:
            chk_file = os.path.join(outdir, f"{stage_name}_{int(time_ns)}ns.chk")
            with open(chk_file, 'wb') as f:
                f.write(context.createCheckpoint())
        
        if time.time() - last_report > 60:
            frac = (step + report_interval) / steps
            elapsed = time.time() - start
            eta = elapsed / frac - elapsed if frac > 0 else 0
            ns_per_hour = (step + report_interval) * step_size / 1000 / (elapsed / 3600)
            print(f"  [{frac*100:5.1f}%] {time_ns:.2f} ns | "
                  f"{ns_per_hour:.2f} ns/h | ETA: {eta/3600:.1f} h")
            last_report = time.time()
    
    # Final checkpoint
    final_chk = os.path.join(outdir, f"{stage_name}_final.chk")
    with open(final_chk, 'wb') as f:
        f.write(context.createCheckpoint())
    
    log.close()
    print(f"\n[Done] Production complete: {time_ns:.3f} ns")
    print(f"  Final checkpoint: {final_chk}")
    print(f"  DCD: {dcd_file}")
    print(f"  Log: {log_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Restart production MD from checkpoint')
    parser.add_argument('--prmtop', required=True, help='Amber prmtop file')
    parser.add_argument('--pdb', required=True, help='Starting PDB coordinates')
    parser.add_argument('--checkpoint', required=True, help='Checkpoint file to restart from')
    parser.add_argument('--name', default='prod', help='Run name')
    parser.add_argument('--outdir', default='.', help='Output directory')
    parser.add_argument('--prod-ns', type=float, default=200, help='Production length (ns)')
    parser.add_argument('--platform', default='auto', choices=['auto', 'CUDA', 'OpenCL', 'CPU', 'Reference'])
    parser.add_argument('--seed', type=int, default=None)
    args = parser.parse_args()
    
    run_production(
        args.prmtop, args.pdb, args.checkpoint, args.name, args.outdir,
        prod_ns=args.prod_ns, seed=args.seed,
    )
