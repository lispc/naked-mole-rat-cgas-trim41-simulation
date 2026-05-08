#!/usr/bin/env python3
"""Run production MD for DNA-bound cGAS+SPRY system (protein only, DNA stripped)."""
import sys
import time
import argparse
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import app
from openmm import unit

BASE = Path(__file__).resolve().parent.parent.parent


def run_md(pdb_path, outdir, name, prod_ns=50, gpu_idx=0):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading system from {pdb_path}...", flush=True)
    pdb = app.PDBFile(str(pdb_path))

    # Force field
    ff = app.ForceField("amber19/protein.ff19SB.xml", "amber19/tip4pew.xml")

    # Build system from minimized PDB
    modeller = app.Modeller(pdb.topology, pdb.positions)
    print(f"  Atoms: {modeller.topology.getNumAtoms()}", flush=True)

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        ewaldErrorTolerance=0.0005,
    )
    print(f"  System created", flush=True)

    temperature = 300 * unit.kelvin
    platform = mm.Platform.getPlatformByName("CUDA")
    props = {"CudaDeviceIndex": str(gpu_idx), "Precision": "mixed"}
    integrator = mm.LangevinMiddleIntegrator(
        temperature, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )

    simulation = app.Simulation(modeller.topology, system, integrator, platform, props)
    simulation.context.setPositions(modeller.positions)

    # ── Re-minimize ──
    print("Re-minimizing...", flush=True)
    simulation.minimizeEnergy(maxIterations=2000)

    # ── Heating (100ps) ──
    print("Heating...", flush=True)
    n_heat = 50000
    simulation.context.setVelocitiesToTemperature(0.1 * unit.kelvin)
    simulation.reporters.append(
        app.StateDataReporter(
            str(outdir / f"{name}_heating.log"), 5000,
            step=True, temperature=True, potentialEnergy=True
        )
    )
    for i in range(n_heat):
        frac = (i + 1) / n_heat
        integrator.setTemperature(temperature * frac)
        simulation.step(1)
    simulation.reporters.clear()

    # ── NPT Equilibration (200ps) ──
    print("NPT equil...", flush=True)
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature, 25))
    simulation.context.reinitialize(preserveState=True)
    n_npt = 100000
    simulation.reporters.append(
        app.StateDataReporter(
            str(outdir / f"{name}_npt.log"), 10000,
            step=True, temperature=True, density=True, potentialEnergy=True
        )
    )
    for _ in range(n_npt):
        simulation.step(1)
    simulation.reporters.clear()

    # ── NVT Production ──
    print(f"NVT production ({prod_ns}ns)...", flush=True)
    # Remove barostat
    forces = list(system.getForces())
    for i, f in enumerate(forces):
        if isinstance(f, mm.MonteCarloBarostat):
            system.removeForce(i)
            break
    simulation.context.reinitialize(preserveState=True)

    n_prod = int(prod_ns * 500000)
    report_interval = 50000  # every 100ps
    checkpoint_interval = int(5e6)  # every 10ns

    simulation.reporters.append(app.DCDReporter(str(outdir / f"{name}_prod.dcd"), 5000))
    simulation.reporters.append(
        app.StateDataReporter(
            str(outdir / f"{name}_prod.log"), report_interval,
            step=True, time=True, speed=True, temperature=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
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
        print(f"  {ns_done}ns, {ns_day:.1f} ns/day", flush=True)

    remaining = n_prod % checkpoint_interval
    if remaining > 0:
        simulation.step(remaining)

    simulation.saveCheckpoint(str(outdir / f"{name}_prod_final.chk"))
    print(f"Done: {outdir}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prod-ns", type=float, default=50)
    parser.add_argument("--gpu", type=int, default=0)
    args = parser.parse_args()
    run_md(args.pdb, args.outdir, args.name, args.prod_ns, args.gpu)


if __name__ == "__main__":
    main()
