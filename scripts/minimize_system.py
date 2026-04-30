#!/usr/bin/env python3
"""Energy minimize a system using OpenMM."""
import argparse
import openmm as mm
from openmm import app

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prmtop', required=True)
    parser.add_argument('--pdb', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--max-iterations', type=int, default=500)
    args = parser.parse_args()
    
    prmtop = app.AmberPrmtopFile(args.prmtop)
    pdb = app.PDBFile(args.pdb)
    
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*mm.unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )
    
    platform = mm.Platform.getPlatformByName('CUDA')
    props = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': '0'}
    
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    context = mm.Context(system, integrator, platform, props)
    context.setPositions(pdb.positions)
    
    init_state = context.getState(getEnergy=True)
    print(f"Initial energy: {init_state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole):.1f} kJ/mol")
    
    print(f"Minimizing (max {args.max_iterations} iterations)...")
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=args.max_iterations)
    
    final_state = context.getState(getEnergy=True, getPositions=True)
    print(f"Final energy: {final_state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole):.1f} kJ/mol")
    
    # Save minimized PDB
    app.PDBFile.writeFile(prmtop.topology, final_state.getPositions(), open(args.out, 'w'))
    print(f"Saved minimized structure to {args.out}")

if __name__ == '__main__':
    main()
