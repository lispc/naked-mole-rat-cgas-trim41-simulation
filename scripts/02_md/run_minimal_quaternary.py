#!/usr/bin/env python3
"""
Minimize + production MD for minimal quaternary (E2~Ub + cGAS).

The prmtop has all protein residues in a single chain. Residue ranges:
  E2 (UBE2D1): residues 1-145
  Ub: residues 146-221 (76 residues)
  cGAS: residues 222-544 (323 residues, starts with ASP200)

A COM distance flat-bottom restraint keeps E2~Ub and cGAS from drifting
too far apart while allowing K315 to sample distances to Ub-G76.
"""
import sys
import argparse
from pathlib import Path
from datetime import datetime
import openmm as mm
from openmm import app, unit
import numpy as np

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")


def log(msg, logfile=None):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    if logfile:
        with open(logfile, 'a') as f:
            f.write(line + '\n')


def make_integrator():
    return mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)


def build_simulation(topology, system, platform, properties):
    """Create a fresh Simulation with a new integrator."""
    return app.Simulation(topology, system, make_integrator(), platform, properties)


def get_protein_residues(topology):
    """Return protein residues (excluding water/ions)."""
    return [r for r in topology.residues() if r.name not in ('HOH', 'WAT', 'Na+', 'Cl-', 'NA', 'CL')]


def get_atom_range(topology, residues, start_res, end_res):
    """Get atom indices for residue range [start_res, end_res)."""
    indices = []
    for r in residues[start_res:end_res]:
        for atom in r.atoms():
            indices.append(atom.index)
    return indices


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', required=True)
    parser.add_argument('--prmtop', required=True)
    parser.add_argument('--inpcrd', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--prod-ns', type=float, default=50)
    parser.add_argument('--gpu', default='0')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logfile = outdir / f"{args.name}.log"

    log(f"=== Minimal Quaternary MD: {args.name} ===", logfile)

    # Load
    log("Loading prmtop/inpcrd...", logfile)
    prmtop = app.AmberPrmtopFile(args.prmtop)
    inpcrd = app.AmberInpcrdFile(args.inpcrd)
    topology = prmtop.topology

    # Find protein residues and component ranges
    protein_res = get_protein_residues(topology)
    log(f"Protein residues: {len(protein_res)}", logfile)

    # E2: ~145 residues, Ub: ~76, cGAS: ~323. Total = 544
    # Verify by checking Ub G76 and cGAS start
    ub_start = 145
    cgas_start = 221
    log(f"E2: [0-{ub_start}), Ub: [{ub_start}-{cgas_start}), cGAS: [{cgas_start}-{len(protein_res)})", logfile)

    e2ub_atoms = get_atom_range(topology, protein_res, 0, cgas_start)
    cgas_atoms = get_atom_range(topology, protein_res, cgas_start, len(protein_res))
    log(f"E2+Ub atoms: {len(e2ub_atoms)}, cGAS atoms: {len(cgas_atoms)}", logfile)

    # Find cGAS K315 NZ atom for analysis
    k315_nz_idx = None
    for r in protein_res[cgas_start:]:
        if r.name == 'LYS' and '315' in str(r.id):
            for atom in r.atoms():
                if atom.name == 'NZ':
                    k315_nz_idx = atom.index
                    break
            break
    if k315_nz_idx is not None:
        log(f"cGAS K315 NZ atom index: {k315_nz_idx}", logfile)

    # ---- Build system ----
    log("Building system...", logfile)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    # Flat-bottom COM distance restraint (30-60 Å, k=50 kJ/mol/nm²)
    restraint = mm.CustomCentroidBondForce(2,
        "0.5*k*max(0,distance(g1,g2)-upper)^2 + 0.5*k*max(0,lower-distance(g1,g2))^2")
    restraint.addGlobalParameter("k", 50 * unit.kilojoules_per_mole / unit.nanometer**2)
    restraint.addGlobalParameter("lower", 30 * unit.angstrom)
    restraint.addGlobalParameter("upper", 60 * unit.angstrom)

    masses = [system.getParticleMass(i).value_in_unit(unit.amu) for i in range(system.getNumParticles())]
    g1 = restraint.addGroup(e2ub_atoms, [masses[i] for i in e2ub_atoms])
    g2 = restraint.addGroup(cgas_atoms, [masses[i] for i in cgas_atoms])
    restraint.addBond([g1, g2], [])
    system.addForce(restraint)
    log("Added COM flat-bottom restraint: 30-60 Å", logfile)

    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': args.gpu, 'Precision': 'mixed'}

    # ---- Minimization with backbone restraints ----
    log("Phase 1: Minimization with backbone restraints...", logfile)

    # Add backbone position restraints
    bb_force = mm.CustomExternalForce("0.5*k_bb*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    bb_force.addGlobalParameter("k_bb", 100 * unit.kilojoules_per_mole / unit.nanometer**2)
    bb_force.addPerParticleParameter("x0")
    bb_force.addPerParticleParameter("y0")
    bb_force.addPerParticleParameter("z0")

    bb_indices = []
    for atom in topology.atoms():
        if atom.name in ('CA', 'C', 'N'):
            bb_indices.append(atom.index)

    sim = build_simulation(topology, system, platform, properties)
    sim.context.setPositions(inpcrd.positions)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    initial_pos = sim.context.getState(getPositions=True).getPositions()
    for idx in bb_indices:
        pos = initial_pos[idx]
        bb_force.addParticle(idx, [pos[0], pos[1], pos[2]])
    bb_force_idx = system.addForce(bb_force)

    sim = build_simulation(topology, system, platform, properties)
    sim.context.setPositions(inpcrd.positions)

    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    log(f"  Initial energy: {e0}", logfile)

    sim.minimizeEnergy(maxIterations=5000)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    log(f"  After minimization: {e1}", logfile)
    log(f"  Delta: {e1 - e0}", logfile)

    # Remove backbone restraints
    system.removeForce(bb_force_idx)
    log("  Removed backbone restraints", logfile)

    # Save minimized structure
    min_pos = sim.context.getState(getPositions=True).getPositions()
    pdb_path = outdir / f"{args.name}_minimized.pdb"
    with open(pdb_path, 'w') as f:
        app.PDBFile.writeFile(topology, min_pos, f)
    log(f"  Saved: {pdb_path}", logfile)

    # ---- Production ----
    log(f"Phase 2: Production ({args.prod_ns} ns)...", logfile)

    sim = build_simulation(topology, system, platform, properties)
    sim.context.setPositions(min_pos)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    # Reporters
    dcd_path = outdir / f"{args.name}.dcd"
    sim.reporters.append(app.DCDReporter(str(dcd_path), 1000))

    state_path = outdir / f"{args.name}_state.log"
    sim.reporters.append(app.StateDataReporter(
        str(state_path), 10000,
        step=True, time=True, potentialEnergy=True, temperature=True,
        speed=True,
    ))

    n_steps = int(args.prod_ns * 1e6 / 2.0)
    chunk_steps = int(10e6 / 2.0)  # ~10 ns chunks

    for chunk_start in range(0, n_steps, chunk_steps):
        steps = min(chunk_steps, n_steps - chunk_start)
        sim.step(steps)
        ns_done = (chunk_start + steps) * 2.0 / 1e6
        chk = outdir / f"{args.name}_{ns_done:.0f}ns.chk"
        sim.saveCheckpoint(str(chk))
        log(f"  {ns_done:.0f} ns checkpoint saved", logfile)

    log(f"Production complete: {args.prod_ns} ns", logfile)


if __name__ == '__main__':
    main()
