#!/usr/bin/env python3
"""Quick energy minimization for quaternary MVP complex."""
import sys
from pathlib import Path
from openmm import app
import openmm as mm
from openmm import unit

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
prmtop_path = BASE / "data/structures/quaternary_mvp/quaternary_mvp.prmtop"
inpcrd_path = BASE / "data/structures/quaternary_mvp/quaternary_mvp.inpcrd"
out_pdb = BASE / "data/structures/quaternary_mvp/quaternary_mvp_minimized.pdb"

print("Loading system...")
prmtop = app.AmberPrmtopFile(str(prmtop_path))
inpcrd = app.AmberInpcrdFile(str(inpcrd_path))

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds,
)

integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond
)

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '2', 'CudaPrecision': 'mixed'}

simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print(f"Initial energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")

print("Minimizing (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)

print(f"Final energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")

# Save minimized structure
state = simulation.context.getState(getPositions=True)
with open(out_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

print(f"Saved minimized structure: {out_pdb}")
