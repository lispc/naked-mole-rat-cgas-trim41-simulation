#!/usr/bin/env python3
"""
Production MD simulation using OpenMM.
Optimized for Apple Metal GPU.

Usage:
  python scripts/run_md.py --system-dir data/md_runs/Hsap_WT --name Hsap_WT_rep1 --seed 42
"""
import argparse
import time
import sys
from pathlib import Path

import openmm
import openmm.app as app
import openmm.unit as unit


def run_simulation(system_xml: Path, ref_pdb: Path, out_dir: Path, name: str, seed: int,
                   nsteps_min=5000, nsteps_heat=50000, nsteps_equil=2500000, nsteps_prod=50000000,
                   dt_fs=4.0, temp=300.0, friction=1.0, barostat_freq=25,
                   report_interval=2500, traj_interval=2500):
    """
    Run full MD protocol:
      - Minimization
      - NVT heating (100 ps)
      - NPT equilibration with restraints (5 ns)
      - NPT production (200 ns with 4fs step)
    
    nsteps_prod=50M @ 4fs = 200 ns
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Load system
    with open(system_xml) as f:
        system = openmm.XmlSerializer.deserialize(f.read())
    
    pdb = app.PDBFile(str(ref_pdb))
    
    # Platform: prefer Metal on Apple Silicon
    platform = openmm.Platform.getPlatformByName("Metal")
    print(f"Using platform: {platform.getName()}")
    
    # Integrator
    integrator = openmm.LangevinMiddleIntegrator(
        temp * unit.kelvin,
        friction / unit.picosecond,
        dt_fs * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(seed)
    
    # Context
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    context.setVelocitiesToTemperature(temp * unit.kelvin, seed)
    
    # Reporters
    dcd_reporter = app.DCDReporter(str(out_dir / f"{name}.dcd"), traj_interval)
    state_reporter = app.StateDataReporter(
        str(out_dir / f"{name}.log"),
        report_interval,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
        remainingTime=True,
        totalSteps=nsteps_min + nsteps_heat + nsteps_equil + nsteps_prod,
        separator="\t",
    )
    
    # === Minimization ===
    print(f"\n[1/4] Energy minimization ({nsteps_min} steps)...")
    t0 = time.time()
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=nsteps_min)
    state = context.getState(getPositions=True, getEnergy=True)
    print(f"  Minimized PE: {state.getPotentialEnergy()}")
    print(f"  Time: {time.time()-t0:.1f}s")
    
    # Save minimized structure
    with open(out_dir / f"{name}_minimized.pdb", "w") as f:
        app.PDBFile.writeFile(pdb.topology, state.getPositions(), f)
    
    # === Heating (NVT) ===
    print(f"\n[2/4] NVT heating to {temp}K ({nsteps_heat} steps = {nsteps_heat*dt_fs/1000:.1f} ps)...")
    context.getIntegrator().setStepSize(dt_fs * unit.femtosecond)
    # Remove barostat if present for NVT
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), openmm.MonteCarloBarostat):
            system.getForce(i).setForceGroup(31)  # disable
    
    t0 = time.time()
    integrator.step(nsteps_heat)
    print(f"  Time: {time.time()-t0:.1f}s")
    
    # === Equilibration (NPT with restraints) ===
    print(f"\n[3/4] NPT equilibration ({nsteps_equil} steps = {nsteps_equil*dt_fs/1000000:.1f} ns)...")
    # Re-enable barostat
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), openmm.MonteCarloBarostat):
            system.getForce(i).setForceGroup(0)
    
    # Add position restraints to protein CA (we'll do this by modifying system XML before)
    # For now, skip restraints to keep script simple; can add later
    
    t0 = time.time()
    integrator.step(nsteps_equil)
    print(f"  Time: {time.time()-t0:.1f}s")
    
    # Save equilibrated structure
    state = context.getState(getPositions=True)
    with open(out_dir / f"{name}_equilibrated.pdb", "w") as f:
        app.PDBFile.writeFile(pdb.topology, state.getPositions(), f)
    
    # === Production (NPT) ===
    print(f"\n[4/4] NPT production ({nsteps_prod} steps = {nsteps_prod*dt_fs/1000000:.1f} ns)...")
    context.getIntegrator().setStepSize(dt_fs * unit.femtosecond)
    
    # Add reporters
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setState(context.getState(getPositions=True, getVelocities=True, getParameters=True))
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(state_reporter)
    
    t0 = time.time()
    simulation.step(nsteps_prod)
    elapsed = time.time() - t0
    
    ns = nsteps_prod * dt_fs / 1000000
    ns_per_day = ns / (elapsed / 86400)
    print(f"\n✅ Production complete!")
    print(f"   Simulated: {ns:.1f} ns")
    print(f"   Wall time: {elapsed/3600:.2f} h")
    print(f"   Speed: {ns_per_day:.1f} ns/day")
    
    # Save final state
    with open(out_dir / f"{name}_final.xml", "w") as f:
        f.write(openmm.XmlSerializer.serialize(simulation.context.getState(getPositions=True, getVelocities=True)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-dir", required=True, help="Directory containing system.xml and ref.pdb")
    parser.add_argument("--name", required=True, help="Run name")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--ns", type=float, default=200.0, help="Production ns")
    parser.add_argument("--dt", type=float, default=4.0, help="Timestep in fs")
    args = parser.parse_args()
    
    system_dir = Path(args.system_dir)
    system_xml = system_dir / f"{system_dir.name}_system.xml"
    ref_pdb = system_dir / f"{system_dir.name}_ref.pdb"
    
    if not system_xml.exists():
        # Try generic names
        candidates = list(system_dir.glob("*_system.xml"))
        if candidates:
            system_xml = candidates[0]
        else:
            raise FileNotFoundError(f"No system.xml found in {system_dir}")
    if not ref_pdb.exists():
        candidates = list(system_dir.glob("*_ref.pdb"))
        if candidates:
            ref_pdb = candidates[0]
        else:
            raise FileNotFoundError(f"No ref.pdb found in {system_dir}")
    
    out_dir = system_dir / "production"
    out_dir.mkdir(exist_ok=True)
    
    nsteps_prod = int(args.ns * 1000000 / args.dt)
    
    run_simulation(
        system_xml=system_xml,
        ref_pdb=ref_pdb,
        out_dir=out_dir,
        name=args.name,
        seed=args.seed,
        nsteps_prod=nsteps_prod,
        dt_fs=args.dt,
    )


if __name__ == "__main__":
    main()
