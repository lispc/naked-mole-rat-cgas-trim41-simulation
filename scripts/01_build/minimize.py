#!/usr/bin/env python3
"""Energy minimization for cGAS-TRIM41 systems using OpenMM.

Supports:
  - Generic system (any prmtop/pdb)
  - S305E system (auto-detects from name)
"""

import argparse
import sys
from pathlib import Path

from openmm import app, CustomExternalForce
import openmm as mm
from openmm import unit


def minimize(prmtop_path, pdb_path, out_prefix, platform_name="CUDA", max_iterations=5000):
    """Minimize a system and save coordinates."""
    prmtop = app.AmberPrmtopFile(str(prmtop_path))
    pdb = app.PDBFile(str(pdb_path))
    
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )
    
    platform = mm.Platform.getPlatformByName(platform_name)
    integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picosecond)
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    print(f"Initial energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    simulation.minimizeEnergy(maxIterations=max_iterations)
    print(f"Final energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    positions = state.getPositions()
    
    out_pdb = Path(out_prefix).with_suffix(".pdb")
    out_rst = Path(out_prefix).with_suffix(".rst7")
    
    with open(out_pdb, "w") as f:
        app.PDBFile.writeFile(simulation.topology, positions, f)
    
    # Amber restart format via parmed
    try:
        import parmed as pmd
        structure = pmd.load_file(str(prmtop_path), str(pdb_path))
        structure.positions = positions
        structure.save(str(out_rst), overwrite=True)
        print(f"Saved: {out_pdb}, {out_rst}")
    except ImportError:
        print(f"Saved PDB only (parmed not available): {out_pdb}")


def main():
    parser = argparse.ArgumentParser(description="Energy minimize a system")
    parser.add_argument("--prmtop", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--platform", default="CUDA", choices=["CUDA", "OpenCL", "CPU"])
    parser.add_argument("--max-iterations", type=int, default=5000)
    args = parser.parse_args()
    
    minimize(args.prmtop, args.pdb, args.out_prefix, args.platform, args.max_iterations)


if __name__ == "__main__":
    main()
