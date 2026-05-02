#!/usr/bin/env python3
"""Verify OpenMM installation and GPU backend (CUDA/OpenCL/Metal)."""
import openmm
import openmm.app as app
import openmm.unit as unit

print(f"OpenMM version: {openmm.__version__}")
print(f"Number of platforms: {openmm.Platform.getNumPlatforms()}")

available = []
for i in range(openmm.Platform.getNumPlatforms()):
    p = openmm.Platform.getPlatform(i)
    print(f"  Platform {i}: {p.getName()}")
    available.append(p.getName())

# Quick benchmark: alanine dipeptide in vacuum
print("\nRunning quick benchmark (alanine dipeptide, 10k steps)...")
from openmmtools.testsystems import AlanineDipeptideVacuum

testsystem = AlanineDipeptideVacuum()
system = testsystem.system

# Try platforms in order of preference
for platform_name in ['CUDA', 'OpenCL', 'Metal', 'CPU']:
    if platform_name not in available:
        continue
    
    integrator = openmm.LangevinMiddleIntegrator(
        300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picosecond
    )
    platform = openmm.Platform.getPlatformByName(platform_name)
    
    props = {}
    if platform_name == 'CUDA':
        props['CudaPrecision'] = 'mixed'
    elif platform_name == 'OpenCL':
        props['OpenCLPrecision'] = 'single'
    
    try:
        context = openmm.Context(system, integrator, platform, props)
        context.setPositions(testsystem.positions)
        
        import time
        t0 = time.time()
        integrator.step(10000)
        dt = time.time() - t0
        
        state = context.getState(getEnergy=True)
        print(f"\n✅ {platform_name} backend: 10k steps in {dt:.2f}s ({10000/dt:.0f} steps/s)")
        print(f"   Potential energy: {state.getPotentialEnergy()}")
        del context
        break  # Stop at first successful GPU platform
    except Exception as e:
        print(f"\n❌ {platform_name} failed: {e}")
        del integrator

print("\nVerification complete.")
