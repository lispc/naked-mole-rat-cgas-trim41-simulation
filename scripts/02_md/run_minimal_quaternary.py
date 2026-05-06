#!/usr/bin/env python3
"""
Minimize + production MD for quaternary v2: RING-E2~Ub + cGAS.

Topology (all protein in one chain after tleap):
  RING1:    residues 0-78    (~79 aa, 5FER chain A)
  RING2:    residues 79-158  (~80 aa, 5FER chain D)
  E2:       residues 159-303 (~145 aa, 5FER chain B)
  Ub:       residues 304-379 (~76 aa, 5FER chain C)
  cGAS:     residues 380-702 (~323 aa)

Key restraints:
  - Isopeptide bond: E2 K85 NZ ↔ Ub G76 C (harmonic, k=500, r0=1.5 Å)
  - COM flat-bottom: RING+E2+Ub ↔ cGAS (30-60 Å)
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
    return mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)


def build_sim(topology, system, platform, props):
    return app.Simulation(topology, system, make_integrator(), platform, props)


def find_key_atoms(topology):
    """Find key atoms by scanning protein residues.

    Returns dict with atom indices for Ub G76 C, E2 K85 NZ, cGAS K315 NZ,
    cGAS start residue index, and atom ranges for COM restraint.
    """
    protein = [r for r in topology.residues()
               if r.name not in ('HOH', 'WAT', 'Na+', 'Cl-', 'NA', 'CL')]

    # Find cGAS start: ASP after the last GLY of Ub
    # Pattern: ... GLY(Ub76) → ASP(cGAS200)
    cgas_start = None
    ub_g76_c = None
    for i in range(len(protein) - 1):
        if protein[i].name == 'GLY' and protein[i + 1].name == 'ASP':
            # Verify: ASP's id contains '200'
            if '200' in str(protein[i + 1].id) or '385' in str(protein[i + 1].id):
                cgas_start = i + 1
                # Ub G76 C
                for a in protein[i].atoms():
                    if a.name == 'C':
                        ub_g76_c = a.index
                break

    if cgas_start is None:
        # Fallback: use known indices
        cgas_start = 384
        ub_g76_c = 5964

    # E2 K85: LYS near position 247 in protein list
    e2_k85_nz = None
    for i in range(240, 260):
        if i >= len(protein):
            break
        r = protein[i]
        if r.name == 'LYS' and 'NZ' in [a.name for a in r.atoms()]:
            # First LYS after position 240 in E2 region
            for a in r.atoms():
                if a.name == 'NZ':
                    e2_k85_nz = a.index
                    break
            break

    # cGAS K315: offset 315-200=115 from cGAS start
    cgas_k315_nz = None
    k315_idx = cgas_start + 115  # ~499
    for offset in range(-5, 6):
        i = k315_idx + offset
        if 0 <= i < len(protein):
            r = protein[i]
            if r.name == 'LYS' and 'NZ' in [a.name for a in r.atoms()]:
                for a in r.atoms():
                    if a.name == 'NZ':
                        cgas_k315_nz = a.index
                        break
                break

    # Atom ranges for COM restraint
    ring_e2_ub_atoms = []
    cgas_atoms = []
    for atom in topology.atoms():
        # Find which protein residue this atom belongs to
        for i, r in enumerate(protein):
            if atom.residue is r:
                if i < cgas_start:
                    ring_e2_ub_atoms.append(atom.index)
                else:
                    cgas_atoms.append(atom.index)
                break

    return {
        'ub_g76_c': ub_g76_c,
        'e2_k85_nz': e2_k85_nz,
        'cgas_k315_nz': cgas_k315_nz,
        'cgas_start': cgas_start,
        'ring_e2_ub_atoms': ring_e2_ub_atoms,
        'cgas_atoms': cgas_atoms,
        'n_protein': len(protein),
    }


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

    log(f"=== Quaternary v2 MD: {args.name} ===", logfile)
    log(f"Production: {args.prod_ns} ns, GPU: {args.gpu}", logfile)

    # ---- Load ----
    log("Loading system...", logfile)
    prmtop = app.AmberPrmtopFile(args.prmtop)
    inpcrd = app.AmberInpcrdFile(args.inpcrd)
    topology = prmtop.topology

    # Find key atoms
    ka = find_key_atoms(topology)
    log(f"Protein residues: {ka['n_protein']}, cGAS start: {ka['cgas_start']}", logfile)
    log(f"Ub G76 C: {ka['ub_g76_c']}, E2 K85 NZ: {ka['e2_k85_nz']}, K315 NZ: {ka['cgas_k315_nz']}", logfile)
    log(f"RING+E2+Ub atoms: {len(ka['ring_e2_ub_atoms'])}, cGAS atoms: {len(ka['cgas_atoms'])}", logfile)

    assert ka['ub_g76_c'] is not None, "Ub G76 C not found!"
    assert ka['e2_k85_nz'] is not None, "E2 K85 NZ not found!"

    # ---- Build system ----
    log("Building system...", logfile)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    # Isopeptide bond restraint: E2 K85 NZ ↔ Ub G76 C
    # Harmonic, r0=1.5 Å (amide C-N bond), k=500 kJ/mol/nm²
    iso_force = mm.CustomBondForce("0.5*k_iso*(r-r0)^2")
    iso_force.addGlobalParameter("k_iso", 5000 * unit.kilojoules_per_mole / unit.nanometer**2)
    iso_force.addGlobalParameter("r0", 0.135 * unit.nanometer)  # 1.35 Å (amide C-N bond)
    iso_force.addBond(ka['e2_k85_nz'], ka['ub_g76_c'], [])
    system.addForce(iso_force)
    log(f"Isopeptide restraint: K85 NZ ↔ G76 C, r0=1.35 Å, k=5000", logfile)

    # COM flat-bottom restraint
    masses = [system.getParticleMass(i).value_in_unit(unit.amu)
              for i in range(system.getNumParticles())]
    com_force = mm.CustomCentroidBondForce(2,
        "0.5*k_com*max(0,distance(g1,g2)-upper)^2 + 0.5*k_com*max(0,lower-distance(g1,g2))^2")
    com_force.addGlobalParameter("k_com", 50 * unit.kilojoules_per_mole / unit.nanometer**2)
    com_force.addGlobalParameter("lower", 30 * unit.angstrom)
    com_force.addGlobalParameter("upper", 60 * unit.angstrom)
    g1 = com_force.addGroup(ka['ring_e2_ub_atoms'], [masses[i] for i in ka['ring_e2_ub_atoms']])
    g2 = com_force.addGroup(ka['cgas_atoms'], [masses[i] for i in ka['cgas_atoms']])
    com_force.addBond([g1, g2], [])
    system.addForce(com_force)
    log("COM flat-bottom restraint: 30-60 Å", logfile)

    platform = mm.Platform.getPlatformByName('CUDA')
    props = {'CudaDeviceIndex': args.gpu, 'Precision': 'mixed'}

    # ---- Minimization with backbone restraints ----
    log("Phase 1: Minimization with BB restraints...", logfile)

    bb_force = mm.CustomExternalForce("0.5*k_bb*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    bb_force.addGlobalParameter("k_bb", 100 * unit.kilojoules_per_mole / unit.nanometer**2)
    bb_force.addPerParticleParameter("x0")
    bb_force.addPerParticleParameter("y0")
    bb_force.addPerParticleParameter("z0")

    bb_indices = [a.index for a in topology.atoms() if a.name in ('CA', 'C', 'N')]

    sim = build_sim(topology, system, platform, props)
    sim.context.setPositions(inpcrd.positions)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    init_pos = sim.context.getState(getPositions=True).getPositions()
    for idx in bb_indices:
        pos = init_pos[idx]
        bb_force.addParticle(idx, [pos[0], pos[1], pos[2]])
    bb_idx = system.addForce(bb_force)

    sim = build_sim(topology, system, platform, props)
    sim.context.setPositions(inpcrd.positions)

    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    log(f"  Initial energy: {e0}", logfile)

    sim.minimizeEnergy(maxIterations=5000)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    log(f"  After min: {e1}  Δ={e1-e0}", logfile)

    system.removeForce(bb_idx)
    log("  Removed BB restraints", logfile)

    # Save minimized
    min_pos = sim.context.getState(getPositions=True).getPositions()
    pdb_path = outdir / f"{args.name}_minimized.pdb"
    with open(pdb_path, 'w') as f:
        app.PDBFile.writeFile(topology, min_pos, f)
    log(f"  Saved: {pdb_path}", logfile)

    # Check isopeptide distance after minimization
    iso_dist = np.linalg.norm(min_pos[ka['e2_k85_nz']].value_in_unit(unit.angstrom)
                              - min_pos[ka['ub_g76_c']].value_in_unit(unit.angstrom))
    log(f"  Isopeptide distance after min: {iso_dist:.2f} Å", logfile)

    # ---- Production ----
    log(f"Phase 2: Production ({args.prod_ns} ns)...", logfile)

    sim = build_sim(topology, system, platform, props)
    sim.context.setPositions(min_pos)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    dcd_path = outdir / f"{args.name}.dcd"
    sim.reporters.append(app.DCDReporter(str(dcd_path), 1000))

    state_path = outdir / f"{args.name}_state.log"
    sim.reporters.append(app.StateDataReporter(
        str(state_path), 10000,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True,
    ))

    n_steps = int(args.prod_ns * 1e6 / 2.0)
    chunk = int(10e6 / 2.0)

    for start in range(0, n_steps, chunk):
        steps = min(chunk, n_steps - start)
        sim.step(steps)
        ns_done = (start + steps) * 2.0 / 1e6
        chk = outdir / f"{args.name}_{ns_done:.0f}ns.chk"
        sim.saveCheckpoint(str(chk))
        log(f"  {ns_done:.0f} ns checkpoint", logfile)

    log(f"Production complete: {args.prod_ns} ns", logfile)


if __name__ == '__main__':
    main()
