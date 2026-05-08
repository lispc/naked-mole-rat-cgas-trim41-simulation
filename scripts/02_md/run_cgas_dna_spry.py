#!/usr/bin/env python3
"""
Build and run MD for DNA-bound cGAS+SPRY system (Route B).

Uses Boltz-2 model 0 prediction as starting structure.
Force field: ff19SB (protein) + OL21 (DNA) + OPC (water)

Protocol: minimize → heat (100ps) → NPT equil (100ps) → NVT prod (100ns)
"""
import argparse
import sys
import time
import numpy as np
from pathlib import Path

import openmm as mm
from openmm import app
from openmm import unit

BASE = Path(__file__).resolve().parent.parent.parent

# ── Force field ─────────────────────────────────────────────────────────
FORCEFIELDS = [
    "amber19/protein.ff19SB.xml",  # protein
    "amber19/DNA.OL21.xml",        # DNA
    "amber19/tip4pew.xml",         # water (4-point, compatible with Modeller)
]


def build_system(pdb_path, buffer_nm=1.0, ionic_strength=0.15):
    """Build OpenMM system from PDB with protein+DNA+water+ions."""
    print(f"Reading PDB: {pdb_path}")
    pdb = app.PDBFile(str(pdb_path))

    # Detect force field
    ff = app.ForceField(*FORCEFIELDS)

    # Modeller: add hydrogens and solvent
    print("Adding hydrogens and solvent...")
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)

    # Add water (truncated octahedron, 1.0 nm buffer)
    modeller.addSolvent(
        ff,
        model="tip4pew",
        boxShape="octahedron",
        padding=buffer_nm * unit.nanometers,
        ionicStrength=ionic_strength * unit.molar,
        positiveIon="Na+",
        negativeIon="Cl-",
    )

    print(f"  Atoms: {modeller.topology.getNumAtoms()}")
    print(f"  Residues: {modeller.topology.getNumResidues()}")

    # Create system
    print("Creating system...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        ewaldErrorTolerance=0.0005,
    )

    return system, modeller


def run_md(system, modeller, outdir, name, prod_ns=100, gpu_idx=0):
    """Run full MD protocol: minimize → heat → NPT → NVT prod."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    temperature = 300 * unit.kelvin

    # Platform
    platform = mm.Platform.getPlatformByName("CUDA")
    props = {"CudaDeviceIndex": str(gpu_idx), "Precision": "mixed"}

    # Integrator
    integrator = mm.LangevinMiddleIntegrator(
        temperature, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )

    # ── 1. Minimization ─────────────────────────────────────────────
    print(f"\n[1/4] Minimization ({name})")
    simulation = app.Simulation(modeller.topology, system, integrator, platform, props)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=5000)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    e_pot = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"  Post-min PE: {e_pot:.1f} kJ/mol")

    # Save minimized PDB
    pos = state.getPositions()
    with open(outdir / f"{name}_minimized.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, pos, f)

    # ── 2. Heating (NVT, 0→300K, 100ps) ────────────────────────────
    print(f"\n[2/4] Heating ({name})")
    n_heat = 50000  # 100 ps / 2 fs
    simulation.context.setVelocitiesToTemperature(0.1 * unit.kelvin)
    dcd_heat = outdir / f"{name}_heating.dcd"
    simulation.reporters.append(app.DCDReporter(str(dcd_heat), 5000))
    simulation.reporters.append(
        app.StateDataReporter(
            str(outdir / f"{name}_heating.log"), 5000,
            step=True, temperature=True, potentialEnergy=True, volume=True
        )
    )

    for i in range(n_heat):
        frac = (i + 1) / n_heat
        t = temperature * frac
        integrator.setTemperature(t)
        simulation.step(1)

    simulation.reporters.clear()

    # ── 3. NPT Equilibration (100ps) ───────────────────────────────
    print(f"\n[3/4] NPT Equilibration ({name})")
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature, 25))
    simulation.context.reinitialize(preserveState=True)

    n_npt = 50000
    dcd_npt = outdir / f"{name}_npt.dcd"
    simulation.reporters.append(app.DCDReporter(str(dcd_npt), 5000))
    simulation.reporters.append(
        app.StateDataReporter(
            str(outdir / f"{name}_npt.log"), 5000,
            step=True, temperature=True, potentialEnergy=True, density=True, volume=True
        )
    )
    for _ in range(n_npt):
        simulation.step(1)
    simulation.reporters.clear()

    # ── 4. NVT Production (prod_ns ns) ─────────────────────────────
    print(f"\n[4/4] NVT Production ({name}, {prod_ns}ns)")
    # Remove barostat for NVT
    forces = list(system.getForces())
    for i, f in enumerate(forces):
        if isinstance(f, mm.MonteCarloBarostat):
            system.removeForce(i)
            break
    simulation.context.reinitialize(preserveState=True)

    n_prod = int(prod_ns * 500000)  # 2fs timestep → 500k steps/ns
    checkpoint_interval = int(5e6)  # every 10ns
    report_interval = 50000          # every 100ps

    dcd_prod = outdir / f"{name}_prod.dcd"
    log_prod = outdir / f"{name}_prod.log"
    simulation.reporters.append(app.DCDReporter(str(dcd_prod), 5000))
    simulation.reporters.append(
        app.StateDataReporter(
            str(log_prod), report_interval,
            step=True, time=True, speed=True, temperature=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            volume=True,
        )
    )

    start_time = time.time()
    for block in range(n_prod // checkpoint_interval):
        simulation.step(checkpoint_interval)
        chk = outdir / f"{name}_prod_{(block+1)*10}ns.chk"
        simulation.saveCheckpoint(str(chk))
        elapsed = time.time() - start_time
        ns_done = (block + 1) * 10
        ns_day = ns_done / (elapsed / 86400)
        print(f"  {ns_done}ns complete, {ns_day:.1f} ns/day")

    # Final chunk
    remaining = n_prod % checkpoint_interval
    if remaining > 0:
        simulation.step(remaining)

    simulation.saveCheckpoint(str(outdir / f"{name}_prod_final.chk"))
    print(f"\nDone: {name}")


def main():
    parser = argparse.ArgumentParser(description="Build and run DNA-bound cGAS+SPRY MD")
    parser.add_argument("--pdb", required=True, help="Input PDB (Boltz-2 model)")
    parser.add_argument("--name", required=True, help="System name")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--prod-ns", type=float, default=100, help="Production length (ns)")
    parser.add_argument("--gpu", type=int, default=0, help="GPU device index")
    parser.add_argument("--buffer", type=float, default=1.0, help="Solvent buffer (nm)")
    args = parser.parse_args()

    print("=" * 60)
    print(f"Route B: DNA-bound cGAS+SPRY MD - {args.name}")
    print("=" * 60)

    system, modeller = build_system(args.pdb, args.buffer)
    run_md(system, modeller, args.outdir, args.name, args.prod_ns, args.gpu)


if __name__ == "__main__":
    main()
