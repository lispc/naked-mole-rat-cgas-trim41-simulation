#!/usr/bin/env python3
"""
Steered MD → Umbrella Sampling pipeline for K315 NZ ↔ Ub G76 C distance.

Phase 1 — Steered MD: pull K315 NZ from current distance to ~5 Å at ~1 Å/ns.
  Saves snapshots every ~1 Å as window starting structures.

Phase 2 — Umbrella Sampling: for each window, run 10ns with harmonic restraint
  at fixed distance. Uses --phase 2 --window N to run individual windows (parallel).

Usage:
  # Phase 1: Steered MD
  python run_steered_us.py --name us_WT --prmtop ... --inpcrd ... --start-pdb ... \
      --outdir ... --gpu 0 --phase 1

  # Phase 2: Umbrella sampling (one window per GPU invocation)
  python run_steered_us.py --name us_WT_win05 --prmtop ... \
      --outdir ... --gpu 0 --phase 2 --window 5 --target 5.0

FULL quaternary key indices (hardcoded):
  K315 NZ = 11306, Ub G76 C = 5964
  E2 K85 NZ = 3395 (isopeptide)
  RING C-term CA = 1190, SPRY N-term CA = 5971 (CC linker)
  RING+E2+Ub+SPRY = atoms 0-9398, cGAS = atoms 9399-14781
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


def build_system(prmtop_path):
    """Build OpenMM system with all restraints matching FULL quaternary MD."""
    prmtop = app.AmberPrmtopFile(prmtop_path)
    top = prmtop.topology
    sys = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds)

    # Isopeptide restraint
    iso = mm.CustomBondForce("0.5*k_iso*(r-r0)^2")
    iso.addGlobalParameter("k_iso", 5000 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    iso.addGlobalParameter("r0", 0.135 * unit.nanometer)
    iso.addBond(3395, 5964, [])
    sys.addForce(iso)

    # CC linker
    cc = mm.CustomBondForce("0.5*k_cc*max(0,r-cc_up)^2 + 0.5*k_cc*max(0,cc_lo-r)^2")
    cc.addGlobalParameter("k_cc", 100 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    cc.addGlobalParameter("cc_lo", 8.0 * unit.nanometer)
    cc.addGlobalParameter("cc_up", 12.0 * unit.nanometer)
    cc.addBond(1190, 5971, [])
    sys.addForce(cc)

    # COM restraint
    masses = [sys.getParticleMass(i).value_in_unit(unit.amu) for i in range(sys.getNumParticles())]
    com = mm.CustomCentroidBondForce(2,
        "0.5*k_cm*max(0,distance(g1,g2)-cm_up)^2 + 0.5*k_cm*max(0,cm_lo-distance(g1,g2))^2")
    com.addGlobalParameter("k_cm", 50 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    com.addGlobalParameter("cm_lo", 3.0 * unit.nanometer)
    com.addGlobalParameter("cm_up", 6.0 * unit.nanometer)
    g1 = com.addGroup(list(range(0, 9399)), [masses[i] for i in range(0, 9399)])
    g2 = com.addGroup(list(range(9399, 14782)), [masses[i] for i in range(9399, 14782)])
    com.addBond([g1, g2], [])
    sys.addForce(com)

    return sys, top, prmtop


def phase1_steered(args):
    """Steered MD: pull K315 NZ → Ub G76 C from current distance to ~5 Å."""
    od = Path(args.outdir)
    od.mkdir(parents=True, exist_ok=True)
    lf = od / f"{args.name}.log"
    log(f"=== Steered MD: {args.name} ===", lf)

    sys, top, prmtop_obj = build_system(args.prmtop)

    # Load starting positions from DCD (has all atoms including water)
    import MDAnalysis as mda
    u = mda.Universe(args.prmtop, args.start_dcd)
    u.trajectory[-1]
    # Convert MDAnalysis positions (Å) to OpenMM Quantity (nm)
    pos_ang = u.atoms.positions  # shape (N, 3) in Å
    pos_nm_list = [mm.Vec3(float(p[0])/10.0, float(p[1])/10.0, float(p[2])/10.0) for p in pos_ang]
    start_pos = unit.Quantity(pos_nm_list, unit.nanometer)
    log(f"Loaded {len(start_pos)} atom positions from DCD last frame", lf)

    # Current distance
    k3_pos = start_pos[11306].value_in_unit(unit.nanometer)
    ub_pos = start_pos[5964].value_in_unit(unit.nanometer)
    current_dist = float(np.linalg.norm(np.array(k3_pos) - np.array(ub_pos))) * 10  # nm → Å
    log(f"Current K315→Ub: {current_dist:.1f} Å", lf)

    # Steered MD parameters
    TARGET_DIST = 0.5  # nm (5 Å)
    PULL_RATE = 1.0  # Å/ns = 0.1 nm/ns = 0.0001 nm/ps
    PULL_RATE_NM_PS = 0.0001  # nm per ps

    # Compute steered MD duration
    start_nm = current_dist / 10.0
    pull_distance_nm = start_nm - TARGET_DIST
    pull_time_ps = pull_distance_nm / PULL_RATE_NM_PS
    pull_steps = int(pull_time_ps / 0.002)  # 2 fs timestep
    pull_ns = pull_time_ps / 1000
    log(f"Pulling from {start_nm*10:.0f}→{TARGET_DIST*10:.0f} Å, {pull_ns:.0f} ns ({pull_steps} steps)", lf)

    # Add steered restraint: moving harmonic between K315 NZ and Ub G76 C
    # We use a CustomCompoundBondForce with a time-dependent target
    # Actually, simpler: use CustomExternalForce on a virtual site
    # Simplest: use a restraint in the pulling phase and manually update it

    # Approach: use CustomBondForce with a parameter we update periodically
    pull_force = mm.CustomBondForce("0.5*k_pull*(r-r_pull)^2")
    pull_force.addGlobalParameter("k_pull", 500 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    pull_force.addGlobalParameter("r_pull", start_nm * unit.nanometer)
    pull_force.addBond(11306, 5964, [])
    sys.addForce(pull_force)

    plat = mm.Platform.getPlatformByName('CUDA')
    props = {'CudaDeviceIndex': args.gpu, 'Precision': 'mixed'}

    sim = mk_sim(top, sys, plat, props)
    sim.context.setPositions(start_pos)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    # Reporters
    sim.reporters.append(app.DCDReporter(str(od / f"{args.name}_steered.dcd"), 500))
    sim.reporters.append(app.StateDataReporter(str(od / f"{args.name}_steered_state.log"), 5000,
        step=True, time=True, potentialEnergy=True, temperature=True))

    # Chunked pulling with distance checkpoints
    chunk_steps = 5000  # 10 ps per chunk
    n_chunks = pull_steps // chunk_steps
    windows_saved = set()

    for i in range(n_chunks):
        # Update target distance
        elapsed_ps = (i + 1) * chunk_steps * 0.002
        new_target_nm = max(TARGET_DIST, start_nm - elapsed_ps * PULL_RATE_NM_PS)
        sim.context.setParameter("r_pull", new_target_nm * unit.nanometer)

        sim.step(chunk_steps)

        # Check current K315 distance
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions()
        current_k3 = float(np.linalg.norm(
            np.array(pos[11306].value_in_unit(unit.nanometer))
            - np.array(pos[5964].value_in_unit(unit.nanometer)))) * 10

        # Save window starting structure at integer Å distances
        dist_bin = int(round(current_k3))
        if dist_bin not in windows_saved and 5 <= dist_bin <= 22:
            windows_saved.add(dist_bin)
            win_pdb = od / f"{args.name}_window_{dist_bin:02d}A.pdb"
            with open(win_pdb, 'w') as f:
                app.PDBFile.writeFile(top, pos, f)
            log(f"  Window {dist_bin} Å saved: {win_pdb}", lf)

        if (i + 1) % 100 == 0:
            log(f"  Chunk {i+1}/{n_chunks}: target={new_target_nm*10:.1f}Å  actual={current_k3:.1f}Å", lf)

    log(f"Steered MD complete. Windows saved: {sorted(windows_saved)}", lf)
    return sorted(windows_saved)


def phase2_umbrella(args):
    """Umbrella sampling: one window with fixed harmonic restraint."""
    od = Path(args.outdir)
    od.mkdir(parents=True, exist_ok=True)
    lf = od / f"{args.name}.log"
    log(f"=== US Window: {args.name} target={args.target} Å ===", lf)

    sys, top, prmtop_obj = build_system(args.prmtop)

    # Load window starting position
    start_pdb = app.PDBFile(f"{args.outdir}/{args.name.replace(f'_win{args.window:02d}', '')}_window_{args.window:02d}A.pdb")
    start_pos = start_pdb.positions

    target_nm = args.target / 10.0

    # Umbrella restraint: fixed harmonic
    us_force = mm.CustomBondForce("0.5*k_us*(r-r_us)^2")
    us_force.addGlobalParameter("k_us", 500 * unit.kilojoules_per_mole / unit.nanometer ** 2)
    us_force.addGlobalParameter("r_us", target_nm * unit.nanometer)
    us_force.addBond(11306, 5964, [])
    sys.addForce(us_force)

    plat = mm.Platform.getPlatformByName('CUDA')
    props = {'CudaDeviceIndex': args.gpu, 'Precision': 'mixed'}

    sim = mk_sim(top, sys, plat, props)
    sim.context.setPositions(start_pos)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    sim.reporters.append(app.DCDReporter(str(od / f"{args.name}.dcd"), 500))
    sim.reporters.append(app.StateDataReporter(str(od / f"{args.name}_state.log"), 10000,
        step=True, time=True, potentialEnergy=True, temperature=True))

    prod_steps = int(args.prod_ns * 1e6 / 2.0)
    log(f"US production: {args.prod_ns} ns ({prod_steps} steps)", lf)

    chunk = int(5e6 / 2.0)  # 5ns per chunk
    for start in range(0, prod_steps, chunk):
        steps = min(chunk, prod_steps - start)
        sim.step(steps)
        ns_done = (start + steps) * 2.0 / 1e6
        log(f"  {ns_done:.0f}ns", lf)

    log(f"US window {args.window} complete", lf)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--prmtop', required=True)
    p.add_argument('--inpcrd', default=None)
    p.add_argument('--start-dcd', default=None, help='Starting DCD trajectory for steered MD')
    p.add_argument('--outdir', required=True)
    p.add_argument('--gpu', default='0')
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--phase', type=int, required=True, choices=[1, 2])
    p.add_argument('--window', type=int, default=0, help='Window number (phase 2)')
    p.add_argument('--target', type=float, default=0, help='Target distance in Å (phase 2)')
    p.add_argument('--prod-ns', type=float, default=10, help='US production length (ns)')
    a = p.parse_args()

    if a.phase == 1:
        phase1_steered(a)
    else:
        phase2_umbrella(a)


if __name__ == '__main__':
    main()
