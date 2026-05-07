#!/usr/bin/env python3
"""
Production MD for FULL quaternary: RING-E2~Ub + SPRY-cGAS.

Key atom indices (hardcoded from topology scan):
  E2 K85 NZ = 3395
  Ub G76 C  = 5964
  cGAS K315 NZ = 11306
  RING1 C-term CA = 1190
  SPRY N-term CA = 5971
  Ring+E2+Ub+SPRY atoms: 0-9398
  cGAS atoms: 9399-...

Restraints:
  1. Isopeptide: K85 NZ ↔ G76 C, harmonic k=5000, r0=1.35 Å
  2. CC linker: RING1 C-term CA ↔ SPRY N-term CA, flat-bottom 80-120 Å
  3. COM: RING+E2+Ub+SPRY ↔ cGAS, flat-bottom 30-60 Å
"""
import sys, argparse, numpy as np
from pathlib import Path
from datetime import datetime
import openmm as mm
from openmm import app, unit

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")


def log(msg, lf=None):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    if lf:
        with open(lf, 'a') as f:
            f.write(line + '\n')


def mk_int():
    return mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)


def mk_sim(top, sys, plat, props):
    return app.Simulation(top, sys, mk_int(), plat, props)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--prmtop', required=True)
    p.add_argument('--inpcrd', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--prod-ns', type=float, default=50)
    p.add_argument('--gpu', default='0')
    p.add_argument('--seed', type=int, default=42)
    a = p.parse_args()

    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    lf = od / f"{a.name}.log"
    log(f"=== FULL Quaternary MD: {a.name} === ({a.prod_ns}ns, GPU{a.gpu})", lf)

    # Load
    prmtop = app.AmberPrmtopFile(a.prmtop)
    inpcrd = app.AmberInpcrdFile(a.inpcrd)
    top = prmtop.topology

    # Build system
    sys = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)

    # Restraint 1: Isopeptide bond
    iso = mm.CustomBondForce("0.5*k_iso*(r-r0)^2")
    iso.addGlobalParameter("k_iso", 5000 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    iso.addGlobalParameter("r0", 0.135 * unit.nanometer)
    iso.addBond(3395, 5964, [])
    sys.addForce(iso)
    log("R1 Isopeptide: K85 NZ ↔ G76 C  k=5000 r0=1.35Å", lf)

    # Restraint 2: CC linker (RING C-term ↔ SPRY N-term)
    cc = mm.CustomBondForce("0.5*k_cc*max(0,r-cc_up)^2 + 0.5*k_cc*max(0,cc_lo-r)^2")
    cc.addGlobalParameter("k_cc", 100 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    cc.addGlobalParameter("cc_lo", 4.0 * unit.nanometer)   # 40 Å (was 80)
    cc.addGlobalParameter("cc_up", 10.0 * unit.nanometer)  # 100 Å (was 120)
    cc.addBond(1190, 5971, [])
    sys.addForce(cc)
    log("R2 CC linker: RING1 C-term ↔ SPRY N-term  flat-bottom 40-100Å", lf)

    # Restraint 3: COM
    masses = [sys.getParticleMass(i).value_in_unit(unit.amu) for i in range(sys.getNumParticles())]
    com = mm.CustomCentroidBondForce(2,
        "0.5*k_cm*max(0,distance(g1,g2)-cm_up)^2 + 0.5*k_cm*max(0,cm_lo-distance(g1,g2))^2")
    com.addGlobalParameter("k_cm", 50 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    com.addGlobalParameter("cm_lo", 3.0 * unit.nanometer)
    com.addGlobalParameter("cm_up", 6.0 * unit.nanometer)
    # Protein-only atoms: RING+E2+Ub+SPRY = 0-9398, cGAS = 9399-14781
    g1 = com.addGroup(list(range(0, 9399)), [masses[i] for i in range(0, 9399)])
    g2 = com.addGroup(list(range(9399, 14782)), [masses[i] for i in range(9399, 14782)])
    com.addBond([g1, g2], [])
    sys.addForce(com)
    log("R3 COM: RING+E2+Ub+SPRY ↔ cGAS  flat-bottom 30-60Å", lf)

    plat = mm.Platform.getPlatformByName('CUDA')
    props = {'CudaDeviceIndex': a.gpu, 'Precision': 'mixed'}

    # Minimization with BB restraints
    log("Phase 1: Minimization...", lf)
    bb = mm.CustomExternalForce("0.5*k_bb*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    bb.addGlobalParameter("k_bb", 100 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    bb.addPerParticleParameter("x0");
    bb.addPerParticleParameter("y0");
    bb.addPerParticleParameter("z0")
    bb_idx_list = [a.index for a in top.atoms() if a.name in ('CA', 'C', 'N')]

    sim = mk_sim(top, sys, plat, props)
    sim.context.setPositions(inpcrd.positions)
    ip = sim.context.getState(getPositions=True).getPositions()
    for idx in bb_idx_list:
        p = ip[idx];
        bb.addParticle(idx, [p[0], p[1], p[2]])
    bb_fidx = sys.addForce(bb)

    sim2 = mk_sim(top, sys, plat, props)
    sim2.context.setPositions(inpcrd.positions)
    e0 = sim2.context.getState(getEnergy=True).getPotentialEnergy()
    sim2.minimizeEnergy(maxIterations=5000)
    e1 = sim2.context.getState(getEnergy=True).getPotentialEnergy()
    log(f"  E: {e0} → {e1}", lf)

    sys.removeForce(bb_fidx)
    log("  BB restraints removed", lf)

    # Save minimized
    mp = sim2.context.getState(getPositions=True).getPositions()
    pdb_path = od / f"{a.name}_minimized.pdb"
    with open(pdb_path, 'w') as f:
        app.PDBFile.writeFile(top, mp, f)
    log(f"  Saved: {pdb_path}", lf)

    # Check key distances
    mp_ang = np.array([list(p.value_in_unit(unit.angstrom)) for p in mp])
    d_iso = float(np.linalg.norm(mp_ang[3395] - mp_ang[5964]))
    d_k315 = float(np.linalg.norm(mp_ang[11306] - mp_ang[5964]))
    d_cc = float(np.linalg.norm(mp_ang[1190] - mp_ang[5971]))
    log(f"  Isopeptide={d_iso:.2f}Å  K315→Ub={d_k315:.0f}Å  CC={d_cc:.0f}Å", lf)

    # Production
    log(f"Phase 2: Production ({a.prod_ns}ns)...", lf)
    sim3 = mk_sim(top, sys, plat, props)
    sim3.context.setPositions(mp)
    sim3.context.setVelocitiesToTemperature(300 * unit.kelvin, a.seed)

    sim3.reporters.append(app.DCDReporter(str(od / f"{a.name}.dcd"), 1000))
    sim3.reporters.append(app.StateDataReporter(str(od / f"{a.name}_state.log"), 10000,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True))

    n_steps = int(a.prod_ns * 1e6 / 2.0)
    chunk = int(10e6 / 2.0)
    for start in range(0, n_steps, chunk):
        steps = min(chunk, n_steps - start)
        sim3.step(steps)
        ns_done = (start + steps) * 2.0 / 1e6
        sim3.saveCheckpoint(str(od / f"{a.name}_{ns_done:.0f}ns.chk"))
        log(f"  {ns_done:.0f}ns checkpoint", lf)

    log(f"DONE: {a.prod_ns}ns", lf)


if __name__ == '__main__':
    main()
