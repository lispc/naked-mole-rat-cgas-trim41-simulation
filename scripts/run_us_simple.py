#!/usr/bin/env python3
import argparse
import os
import numpy as np
from openmm import app, Platform, CustomCentroidBondForce, LangevinMiddleIntegrator
from openmm.unit import kelvin, picosecond, femtosecond, nanometer

def run_us(prmtop_path, rst7_path, center_A, k, name, outdir, gpu_id, em_steps, prod_ns):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    print(f"US: {name}  center={center_A:.1f}A  k={k}")
    
    prmtop = app.AmberPrmtopFile(prmtop_path)
    inpcrd = app.AmberInpcrdFile(rst7_path)
    
    # Find RING CA and Lys-334 NZ
    topology = prmtop.topology
    ring_ca = []
    lys_nz = None
    for atom in topology.atoms():
        res = atom.residue
        if res.index + 1 <= 43 and atom.name == "CA":
            ring_ca.append(atom.index)
        if res.index + 1 == 334 and atom.name == "NZ":
            lys_nz = atom.index
    
    print(f"  RING CA: {len(ring_ca)}  |  Lys-334 NZ: {lys_nz}")
    
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=app.HBonds,
    )
    
    # Add restraint
    r0 = center_A / 10.0
    force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
    force.addGlobalParameter("k", k)
    force.addGlobalParameter("r0", r0)
    force.addGroup(ring_ca)
    force.addGroup([lys_nz])
    force.addBond([0, 1])
    system.addForce(force)
    
    platform = Platform.getPlatformByName("CUDA")
    integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtosecond)
    sim = app.Simulation(topology, system, integrator, platform)
    sim.context.setPositions(inpcrd.positions)
    
    # Set periodic box vectors from inpcrd if available
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        print(f"  Box: {inpcrd.boxVectors}")
    
    # EM
    print(f"  EM ({em_steps} steps)...")
    sim.minimizeEnergy(maxIterations=em_steps)
    
    # Short warm-up (100ps NVT) 
    print(f"  Warm-up 100ps...")
    sim.step(int(0.1 * 1000 / 0.002))
    
    # Production
    prod_steps = int(prod_ns * 1000 / 0.002)
    dcd_freq = int(1000.0 / 0.002)  # 1ns = 500k steps
    log_freq = int(0.5 / 0.002)
    
    os.makedirs(outdir, exist_ok=True)
    sim.reporters.append(app.DCDReporter(f"{outdir}/{name}.dcd", dcd_freq))
    sim.reporters.append(app.StateDataReporter(
        f"{outdir}/{name}.log", log_freq, step=True, time=True,
        potentialEnergy=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=prod_steps, separator="\t"
    ))
    
    cv = open(f"{outdir}/{name}_cv.dat", "w")
    cv.write("#step\ttime_ns\tcv_A\n")
    
    print(f"  Production {prod_ns}ns...")
    for i in range(0, prod_steps, log_freq):
        sim.step(log_freq)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        ring_com = np.mean(pos[ring_ca], axis=0)
        dist = np.linalg.norm(ring_com - pos[lys_nz]) * 10.0
        step = sim.context.getState().getStepCount()
        t = step * 0.002 / 1000.0
        cv.write(f"{step}\t{t:.3f}\t{dist:.2f}\n")
        cv.flush()
    
    cv.close()
    print(f"Done: {name}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--prmtop", required=True)
    p.add_argument("--rst7", required=True)
    p.add_argument("--center", type=float, required=True)
    p.add_argument("--k", type=float, default=200.0)
    p.add_argument("--name", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--gpu", type=int, default=0)
    p.add_argument("--em-steps", type=int, default=500)
    p.add_argument("--prod-ns", type=float, default=20.0)
    args = p.parse_args()
    run_us(args.prmtop, args.rst7, args.center, args.k, args.name, args.outdir,
           args.gpu, args.em_steps, args.prod_ns)
