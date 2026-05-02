#!/usr/bin/env python3
"""Energy minimization for S305E system."""

from openmm import app, unit, CustomExternalForce, LangevinIntegrator, LocalEnergyMinimizer, Platform

prmtop = app.AmberPrmtopFile('data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E.prmtop')
pdb = app.PDBFile('data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_solvated.pdb')

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds
)

# Add weak position restraints to protein backbone
force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
force.addGlobalParameter("k", 100.0)  # kJ/mol/nm^2
for atom in pdb.topology.atoms():
    if atom.residue.chain.id in ['A', 'B'] and atom.name in ['N', 'CA', 'C', 'O']:
        force.addParticle(atom.index, pdb.positions[atom.index])
system.addForce(force)

integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
platform = Platform.getPlatformByName('CUDA')
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

# Initial energy
state = simulation.context.getState(getEnergy=True)
init_e = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Initial energy: {init_e:.1f} kJ/mol")

# Minimize
print("Minimizing...")
LocalEnergyMinimizer.minimize(simulation.context, maxIterations=5000)

# Final energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
final_e = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
print(f"Final energy: {final_e:.1f} kJ/mol")
print(f"Delta: {final_e - init_e:.1f} kJ/mol")

# Save
minimized_pdb = 'data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E_minimized.pdb'
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
print(f"Saved: {minimized_pdb}")
