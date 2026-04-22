#!/usr/bin/env python3
"""Verify OpenMM installation and Metal GPU backend."""
import openmm
import openmm.app as app
import openmm.unit as unit

print(f"OpenMM version: {openmm.__version__}")
print(f"Number of platforms: {openmm.Platform.getNumPlatforms()}")
for i in range(openmm.Platform.getNumPlatforms()):
    p = openmm.Platform.getPlatform(i)
    print(f"  Platform {i}: {p.getName()}")

# Quick benchmark: alanine dipeptide in vacuum
print("\nRunning quick benchmark (alanine dipeptide, 10k steps)...")
from openmmtools.testsystems import AlanineDipeptideVacuum

try:
    testsystem = AlanineDipeptideVacuum()
    system = testsystem.system
    integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond)
    
    # Try Metal first
    platform = openmm.Platform.getPlatformByName("Metal")
    context = openmm.Context(system, integrator, platform)
    context.setPositions(testsystem.positions)
    
    import time
    t0 = time.time()
    integrator.step(10000)
    dt = time.time() - t0
    print(f"Metal backend: 10k steps in {dt:.2f}s ({10000/dt:.0f} steps/s)")
    state = context.getState(getEnergy=True)
    print(f"Potential energy: {state.getPotentialEnergy()}")
    print("\n✅ OpenMM Metal backend is working!")
except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
