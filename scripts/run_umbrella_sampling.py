#!/usr/bin/env python3
import argparse
import os
from datetime import datetime
import numpy as np
from openmm import app, Platform, CustomCentroidBondForce, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.unit import kelvin, picosecond, femtosecond, nanometer, kilojoule_per_mole


def find_atom_indices(prmtop, ring_range, lys_resid, lys_atom_name="NZ"):
    topology = prmtop.topology
    ring_indices = []
    lys_index = None
    for atom in topology.atoms():
        res = atom.residue
        if res.index + 1 >= ring_range[0] and res.index + 1 <= ring_range[1]:
            ring_indices.append(atom.index)
        if res.index + 1 == lys_resid and atom.name == lys_atom_name:
            lys_index = atom.index
    return ring_indices, lys_index


def make_system_with_restraint(prmtop, k, r0_nm, ring_indices, lys_index):
    system = prmtop.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=1.0*nanometer, constraints=app.HBonds,
    )
    force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
    force.addGlobalParameter("k", k)
    force.addGlobalParameter("r0", r0_nm)
    force.addGroup(ring_indices)
    force.addGroup([lys_index])
    force.addBond([0, 1])
    system.addForce(force)
    return system


def run_us_window(prmtop_path, coord_path, center_A, k, name, outdir, gpu_id=0, 
                  em_steps=500, npt_ns=1.0, prod_ns=20.0):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    print(f"[{datetime.now()}] US Window: {name}  center={center_A:.1f}A  k={k}")
    
    prmtop = app.AmberPrmtopFile(prmtop_path)
    if coord_path.endswith('.rst7') or coord_path.endswith('.ncrst'):
        positions = app.AmberInpcrdFile(coord_path).positions
    else:
        positions = app.PDBFile(coord_path).positions
    
    ring_indices, lys_index = find_atom_indices(prmtop, (1, 43), 334, "NZ")
    r0_nm = center_A / 10.0
    print(f"  RING: {len(ring_indices)} atoms  |  cGAS-Lys315 NZ (top resid 334): {lys_index}")
    
    platform = Platform.getPlatformByName("CUDA")
    
    # EM
    print(f"  EM ({em_steps} steps)...")
    system = make_system_with_restraint(prmtop, k, r0_nm, ring_indices, lys_index)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtosecond)
    sim = app.Simulation(prmtop.topology, system, integrator, platform)
    sim.context.setPositions(positions)
    sim.minimizeEnergy(maxIterations=em_steps)
    state = sim.context.getState(getPositions=True)
    pos_em = state.getPositions()
    del sim, integrator
    
    # NPT
    print(f"  NPT {npt_ns}ns...")
    system2 = make_system_with_restraint(prmtop, k, r0_nm, ring_indices, lys_index)
    system2.addForce(MonteCarloBarostat(1.0, 300))
    integrator2 = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtosecond)
    sim2 = app.Simulation(prmtop.topology, system2, integrator2, platform)
    sim2.context.setPositions(pos_em)
    sim2.step(int(npt_ns * 1000 / 0.002))
    state2 = sim2.context.getState(getPositions=True, getVelocities=True)
    pos_npt = state2.getPositions()
    vel_npt = state2.getVelocities()
    del sim2, integrator2
    
    # NVT Production
    print(f"  Production {prod_ns}ns...")
    system3 = make_system_with_restraint(prmtop, k, r0_nm, ring_indices, lys_index)
    integrator3 = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtosecond)
    sim3 = app.Simulation(prmtop.topology, system3, integrator3, platform)
    sim3.context.setPositions(pos_npt)
    sim3.context.setVelocities(vel_npt)
    
    os.makedirs(outdir, exist_ok=True)
    dcd_freq = int(0.1 / 0.002)
    log_freq = int(0.5 / 0.002)
    prod_steps = int(prod_ns * 1000 / 0.002)
    
    sim3.reporters.append(app.DCDReporter(f"{outdir}/{name}.dcd", dcd_freq))
    sim3.reporters.append(app.StateDataReporter(
        f"{outdir}/{name}.log", log_freq, step=True, time=True,
        potentialEnergy=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=prod_steps, separator="\t"
    ))
    
    cv_file = open(f"{outdir}/{name}_cv.dat", "w")
    cv_file.write("#step\ttime_ns\tcv_nm\tcv_A\n")
    
    for i in range(0, prod_steps, log_freq):
        sim3.step(log_freq)
        state = sim3.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometer)
        ring_com = np.mean(pos[ring_indices], axis=0)
        dist_nm = np.linalg.norm(ring_com - pos[lys_index])
        step = sim3.context.getState().getStepCount()
        time_ns = step * 0.002 / 1000.0
        cv_file.write(f"{step}\t{time_ns:.3f}\t{dist_nm:.4f}\t{dist_nm*10:.2f}\n")
        cv_file.flush()
    
    cv_file.close()
    print(f"[{datetime.now()}] Done: {name}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--prmtop", required=True)
    p.add_argument("--coord", required=True)
    p.add_argument("--center", type=float, required=True)
    p.add_argument("--k", type=float, default=1000.0)
    p.add_argument("--name", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--gpu", type=int, default=0)
    p.add_argument("--em-steps", type=int, default=500)
    p.add_argument("--npt-ns", type=float, default=1.0)
    p.add_argument("--prod-ns", type=float, default=20.0)
    args = p.parse_args()
    run_us_window(args.prmtop, args.coord, args.center, args.k, args.name, args.outdir,
                  args.gpu, args.em_steps, args.npt_ns, args.prod_ns)
