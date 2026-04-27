#!/usr/bin/env python3
"""Production MD with OpenMM for cGAS-TRIM41 systems.

Stages:
  1. Heating:     NVT, 0K -> 300K, 100ps (50k steps @ 2fs)
  2. NPT equil:   NPT, 300K, 100ps (50k steps @ 2fs)
  3. NVT prod:    NVT, 300K, 200ns (100M steps @ 2fs)

Output: DCD trajectory, checkpoint files, energy log

Usage:
    python run_md.py --prmtop data/md_runs/Hgal_domain/Hgal_domain.prmtop \
                     --pdb data/md_runs/Hgal_domain/Hgal_domain_minimized.pdb \
                     --name Hgal_domain_rep1 \
                     --outdir data/md_runs/Hgal_domain \
                     --prod-ns 200
"""

import argparse
import os
import sys
import time

import openmm as mm
from openmm import app


def get_platform(platform_name='auto'):
    """Select best OpenMM platform with appropriate properties.
    
    Args:
        platform_name: 'auto', 'CUDA', 'OpenCL', 'CPU', 'Reference'
    
    Returns:
        (platform, properties_dict)
    """
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
        # RTX 3090 supports mixed precision
        props['CudaPrecision'] = 'mixed'
        props['CudaDeviceIndex'] = '0'
        # Optional: deterministic forces for reproducibility
        # props['DeterministicForces'] = 'true'
    elif platform_name == 'OpenCL':
        # Apple Silicon only supports single precision
        props['OpenCLPrecision'] = 'single'
    
    return platform, props


def setup_system(prmtop, pdb, platform_name='auto'):
    """Create OpenMM system and context."""
    print(f"[Setup] Loading topology from {prmtop}")
    prmtop_obj = app.AmberPrmtopFile(prmtop)
    
    print(f"[Setup] Loading coordinates from {pdb}")
    pdb_obj = app.PDBFile(pdb)
    
    print("[Setup] Creating system...")
    system = prmtop_obj.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * mm.unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    
    platform, props = get_platform(platform_name)
    print(f"[Setup] Platform: {platform.getName()}")
    if props:
        print(f"[Setup] Properties: {props}")
    
    return system, prmtop_obj.topology, pdb_obj.positions, platform, props


def run_heating(system, topology, positions, platform, props, outdir, name,
                target_temp=300, steps=50000, step_size=0.002, report_interval=1000, seed=None):
    """Stage 1: Heating from 0K to target_temp (NVT)."""
    
    stage_name = f"{name}_heating"
    print(f"\n{'='*60}")
    print(f"  Stage 1: HEATING (0K -> {target_temp}K)")
    print(f"{'='*60}")
    
    # Set temperature ramp
    integrator = mm.LangevinMiddleIntegrator(
        target_temp * mm.unit.kelvin,
        1.0 / mm.unit.picosecond,
        step_size * mm.unit.picosecond,
    )
    if seed is not None:
        integrator.setRandomNumberSeed(seed)
    
    context = mm.Context(system, integrator, platform, props)
    context.setPositions(positions)
    
    # Initial energy
    state = context.getState(getEnergy=True)
    e_init = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
    print(f"  Initial energy: {e_init:.1f} kJ/mol")
    
    # Setup reporters
    dcd_file = os.path.join(outdir, f"{stage_name}.dcd")
    log_file = os.path.join(outdir, f"{stage_name}.log")
    
    dcd = app.DCDFile(open(dcd_file, 'wb'), topology, step_size * mm.unit.picosecond)
    
    with open(log_file, 'w') as log:
        log.write("# Step  Time(ps)  Potential(kJ/mol)  Temperature(K)\n")
        
        start = time.time()
        temp_increment = target_temp / (steps / report_interval)
        current_temp = 0.0
        
        for step in range(0, steps, report_interval):
            # Update temperature
            current_temp = min(target_temp, current_temp + temp_increment)
            integrator.setTemperature(current_temp * mm.unit.kelvin)
            
            # Run
            integrator.step(report_interval)
            
            # Report
            state = context.getState(getEnergy=True, getPositions=True)
            e = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
            t = state.getKineticEnergy()
            # Approximate temperature from kinetic energy
            dof = system.getNumParticles() * 3 - system.getNumConstraints()
            temp = (2 * t.value_in_unit(mm.unit.kilojoule_per_mole) / 
                    (dof * mm.unit.MOLAR_GAS_CONSTANT_R.value_in_unit(mm.unit.kilojoule_per_mole / mm.unit.kelvin)))
            
            time_ps = (step + report_interval) * step_size
            log.write(f"{step + report_interval:8d}  {time_ps:8.2f}  {e:14.2f}  {temp:10.2f}\n")
            
            if step % (10 * report_interval) == 0:
                dcd.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
        
        elapsed = time.time() - start
        ns_day = (steps * step_size / 1000) / (elapsed / 86400)
        print(f"  Completed in {elapsed:.1f}s ({ns_day:.1f} ns/day)")
    
    # Save final state
    final_state = context.getState(getPositions=True, getVelocities=True, 
                                    getParameters=True, enforcePeriodicBox=True)
    checkpoint = os.path.join(outdir, f"{stage_name}.chk")
    with open(checkpoint, 'wb') as f:
        f.write(context.createCheckpoint())
    
    final_positions = final_state.getPositions()
    del context
    
    print(f"  Output: {dcd_file}")
    print(f"  Checkpoint: {checkpoint}")
    return final_positions


def run_npt(system, topology, positions, platform, props, outdir, name,
            steps=50000, step_size=0.002, report_interval=1000, seed=None):
    """Stage 2: NPT equilibration at 300K, 1 bar."""
    
    stage_name = f"{name}_npt"
    print(f"\n{'='*60}")
    print(f"  Stage 2: NPT EQUILIBRATION (300K, 1bar)")
    print(f"{'='*60}")
    
    # Add barostat
    barostat = mm.MonteCarloBarostat(1.0 * mm.unit.bar, 300 * mm.unit.kelvin, 25)
    system.addForce(barostat)
    
    integrator = mm.LangevinMiddleIntegrator(
        300 * mm.unit.kelvin,
        1.0 / mm.unit.picosecond,
        step_size * mm.unit.picosecond,
    )
    if seed is not None:
        integrator.setRandomNumberSeed(seed)
    
    context = mm.Context(system, integrator, platform, props)
    context.setPositions(positions)
    
    # Initial energy
    state = context.getState(getEnergy=True)
    e_init = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
    print(f"  Initial energy: {e_init:.1f} kJ/mol")
    
    dcd_file = os.path.join(outdir, f"{stage_name}.dcd")
    log_file = os.path.join(outdir, f"{stage_name}.log")
    
    dcd = app.DCDFile(open(dcd_file, 'wb'), topology, step_size * mm.unit.picosecond)
    
    with open(log_file, 'w') as log:
        log.write("# Step  Time(ps)  Potential(kJ/mol)  Volume(nm^3)\n")
        
        start = time.time()
        for step in range(0, steps, report_interval):
            integrator.step(report_interval)
            
            state = context.getState(getEnergy=True, getPositions=True)
            e = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
            vol = state.getPeriodicBoxVolume().value_in_unit(mm.unit.nanometer**3)
            
            time_ps = (step + report_interval) * step_size
            log.write(f"{step + report_interval:8d}  {time_ps:8.2f}  {e:14.2f}  {vol:12.3f}\n")
            
            if step % (10 * report_interval) == 0:
                dcd.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
        
        elapsed = time.time() - start
        ns_day = (steps * step_size / 1000) / (elapsed / 86400)
        print(f"  Completed in {elapsed:.1f}s ({ns_day:.1f} ns/day)")
    
    final_state = context.getState(getPositions=True, getVelocities=True,
                                    getParameters=True, enforcePeriodicBox=True)
    checkpoint = os.path.join(outdir, f"{stage_name}.chk")
    with open(checkpoint, 'wb') as f:
        f.write(context.createCheckpoint())
    
    final_positions = final_state.getPositions()
    box = final_state.getPeriodicBoxVectors()
    del context
    
    # Remove barostat for production
    system.removeForce(system.getNumForces() - 1)
    
    print(f"  Final box: {[v.value_in_unit(mm.unit.nanometer) for v in box]}")
    print(f"  Output: {dcd_file}")
    print(f"  Checkpoint: {checkpoint}")
    return final_positions, box


def run_production(system, topology, positions, box, platform, props, outdir, name,
                   prod_ns=200, step_size=0.002, report_interval=50000, seed=None):
    """Stage 3: NVT production run."""
    
    steps = int(prod_ns * 1000 / step_size)  # ns -> steps
    stage_name = f"{name}_prod"
    
    print(f"\n{'='*60}")
    print(f"  Stage 3: NVT PRODUCTION ({prod_ns}ns)")
    print(f"{'='*60}")
    print(f"  Total steps: {steps}")
    print(f"  Report every: {report_interval} steps ({report_interval * step_size}ps)")
    
    integrator = mm.LangevinMiddleIntegrator(
        300 * mm.unit.kelvin,
        1.0 / mm.unit.picosecond,
        step_size * mm.unit.picosecond,
    )
    if seed is not None:
        integrator.setRandomNumberSeed(seed)
    
    context = mm.Context(system, integrator, platform, props)
    context.setPositions(positions)
    context.setPeriodicBoxVectors(*box)
    
    dcd_file = os.path.join(outdir, f"{stage_name}.dcd")
    log_file = os.path.join(outdir, f"{stage_name}.log")
    
    # Write every 50ps (25k steps @ 2fs) to keep file size manageable
    dcd = app.DCDFile(open(dcd_file, 'wb'), topology, step_size * mm.unit.picosecond)
    
    # Checkpoint every 1ns
    chk_interval = int(1.0 / step_size)
    
    with open(log_file, 'w') as log:
        log.write("# Step  Time(ns)  Potential(kJ/mol)\n")
        
        start = time.time()
        last_report = start
        
        for step in range(0, steps, report_interval):
            integrator.step(report_interval)
            
            # Energy log (every report_interval)
            state = context.getState(getEnergy=True, getPositions=True)
            e = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
            time_ns = (step + report_interval) * step_size / 1000
            log.write(f"{step + report_interval:10d}  {time_ns:8.3f}  {e:14.2f}\n")
            log.flush()
            
            # DCD write (every report_interval)
            dcd.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
            
            # Checkpoint (every 1ns)
            if (step + report_interval) % chk_interval == 0:
                chk_file = os.path.join(outdir, f"{stage_name}_{int(time_ns)}ns.chk")
                with open(chk_file, 'wb') as f:
                    f.write(context.createCheckpoint())
            
            # Progress report
            if time.time() - last_report > 60:  # Every minute
                elapsed = time.time() - start
                frac = (step + report_interval) / steps
                total_est = elapsed / frac
                remaining = total_est - elapsed
                ns_done = time_ns
                ns_day = ns_done / (elapsed / 86400)
                print(f"  {ns_done:.1f}ns / {prod_ns}ns ({frac*100:.1f}%) | "
                      f"{ns_day:.1f} ns/day | ETA: {remaining/3600:.1f}h")
                last_report = time.time()
        
        elapsed = time.time() - start
        ns_day = prod_ns / (elapsed / 86400)
        print(f"\n  Completed in {elapsed/3600:.1f}h ({ns_day:.1f} ns/day)")
    
    # Final checkpoint
    final_chk = os.path.join(outdir, f"{stage_name}_final.chk")
    with open(final_chk, 'wb') as f:
        f.write(context.createCheckpoint())
    
    del context
    
    print(f"  Output: {dcd_file}")
    print(f"  Log: {log_file}")
    print(f"  Final checkpoint: {final_chk}")


def main():
    parser = argparse.ArgumentParser(description='Production MD for cGAS-TRIM41')
    parser.add_argument('--prmtop', required=True, help='Amber prmtop file')
    parser.add_argument('--pdb', required=True, help='Starting PDB coordinates (minimized)')
    parser.add_argument('--name', default='prod', help='Run name')
    parser.add_argument('--outdir', default='.', help='Output directory')
    parser.add_argument('--prod-ns', type=float, default=200, help='Production length (ns)')
    parser.add_argument('--skip-heat', action='store_true', help='Skip heating')
    parser.add_argument('--skip-npt', action='store_true', help='Skip NPT equilibration')
    parser.add_argument('--platform', default='auto', choices=['auto', 'CUDA', 'OpenCL', 'CPU', 'Reference'],
                        help='OpenMM platform (auto=best available)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for Langevin integrator (replica independence)')
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"{'='*60}")
    print(f"  cGAS-TRIM41 Production MD")
    print(f"{'='*60}")
    print(f"  System:   {args.prmtop}")
    print(f"  Platform: {args.platform}")
    print(f"  Output:   {args.outdir}")
    print()
    
    # Setup
    system, topology, positions, platform, props = setup_system(
        args.prmtop, args.pdb, args.platform
    )
    
    current_positions = positions
    current_box = None
    
    # Stage 1: Heating
    if not args.skip_heat:
        current_positions = run_heating(
            system, topology, current_positions, platform, props,
            args.outdir, args.name, seed=args.seed
        )
    
    # Stage 2: NPT equilibration
    if not args.skip_npt:
        current_positions, current_box = run_npt(
            system, topology, current_positions, platform, props,
            args.outdir, args.name, seed=args.seed
        )
    
    # Stage 3: Production
    if current_box is None:
        # Get box from system if skipped NPT
        state = mm.Context(system, mm.VerletIntegrator(0.001), platform, props).getState(getParameters=True)
        current_box = state.getPeriodicBoxVectors()
    
    run_production(
        system, topology, current_positions, current_box, platform, props,
        args.outdir, args.name, prod_ns=args.prod_ns, seed=args.seed
    )
    
    print(f"\n{'='*60}")
    print(f"  ✅ All stages complete!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
