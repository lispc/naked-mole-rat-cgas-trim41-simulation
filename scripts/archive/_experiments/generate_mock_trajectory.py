#!/usr/bin/env python3
"""
Generate mock MD trajectories for pipeline testing.
Runs a short NVT simulation (default 10ps) from a minimized structure
to produce a multi-frame DCD file for analysis script validation.
"""
import argparse
import os
import sys
from datetime import datetime

try:
    from openmm import app
    import openmm as mm
    from openmm import unit
except ImportError:
    print("Error: OpenMM not available. Activate cgas-md environment.")
    sys.exit(1)


def generate_mock(prmtop_path, coord_path, out_dcd, out_log,
                  duration_ps=10.0, dt_fs=2.0, save_interval=10,
                  temperature=300.0, platform_name='CUDA', gpu_index='3'):
    """
    Run short NVT simulation and write DCD.
    
    Parameters
    ----------
    duration_ps : float
        Total simulation time in picoseconds.
    dt_fs : float
        Timestep in femtoseconds.
    save_interval : int
        Save frame every N steps.
    coord_path : str
        Path to coordinate file (.pdb or .rst7/.inpcrd)
    """
    print(f"[{datetime.now()}] Loading topology: {prmtop_path}")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    
    if coord_path.lower().endswith('.pdb'):
        inpcrd = app.PDBFile(coord_path)
    else:
        inpcrd = app.AmberInpcrdFile(coord_path)

    print(f"[{datetime.now()}] Creating system...")
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )

    integrator = mm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1.0 / unit.picosecond,
        dt_fs * unit.femtoseconds,
    )

    platform = mm.Platform.getPlatformByName(platform_name)
    properties = {}
    if platform_name == 'CUDA':
        properties['CudaDeviceIndex'] = gpu_index
        properties['CudaPrecision'] = 'mixed'

    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)
    
    # Handle periodic box vectors for both PDB and inpcrd formats
    if hasattr(inpcrd, 'boxVectors') and inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    elif hasattr(inpcrd, 'getTopology'):
        ucd = inpcrd.getTopology().getUnitCellDimensions()
        if ucd is not None:
            a, b, c = ucd.value_in_unit(unit.nanometer)
            box_vectors = [
                mm.Vec3(a, 0, 0),
                mm.Vec3(0, b, 0),
                mm.Vec3(0, 0, c),
            ]
            simulation.context.setPeriodicBoxVectors(*box_vectors)

    # Quick energy report
    state = simulation.context.getState(getEnergy=True)
    print(f"  Initial energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

    # Add reporters
    n_steps = int((duration_ps * 1000) / dt_fs)
    simulation.reporters.append(app.DCDReporter(out_dcd, save_interval))
    simulation.reporters.append(
        app.StateDataReporter(
            out_log, save_interval * 10,
            step=True, time=True, potentialEnergy=True, temperature=True,
            speed=True, remainingTime=True, totalSteps=n_steps,
        )
    )

    print(f"[{datetime.now()}] Running {duration_ps} ps = {n_steps} steps (save every {save_interval})")
    simulation.step(n_steps)

    state = simulation.context.getState(getEnergy=True)
    print(f"[{datetime.now()}] Final energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")
    print(f"[{datetime.now()}] Wrote {out_dcd} and {out_log}")


def main():
    parser = argparse.ArgumentParser(description="Generate mock MD trajectory")
    parser.add_argument('--prmtop', required=True)
    parser.add_argument('--coord', required=True, help='Coordinate file (.pdb or .rst7/.inpcrd)')
    parser.add_argument('--out-dcd', required=True)
    parser.add_argument('--out-log', required=True)
    parser.add_argument('--duration-ps', type=float, default=10.0)
    parser.add_argument('--dt-fs', type=float, default=2.0)
    parser.add_argument('--save-interval', type=int, default=10)
    parser.add_argument('--platform', default='CUDA', choices=['CUDA', 'OpenCL', 'CPU'])
    parser.add_argument('--gpu', default='3')
    args = parser.parse_args()

    generate_mock(
        args.prmtop, args.coord, args.out_dcd, args.out_log,
        duration_ps=args.duration_ps, dt_fs=args.dt_fs,
        save_interval=args.save_interval,
        platform_name=args.platform, gpu_index=args.gpu,
    )


if __name__ == '__main__':
    main()
