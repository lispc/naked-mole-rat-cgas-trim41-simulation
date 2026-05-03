#!/usr/bin/env python3
"""Production MD for quaternary MVP (50 ns test)."""
import sys
import argparse
from pathlib import Path
from datetime import datetime
from openmm import app
import openmm as mm
from openmm import unit

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")

parser = argparse.ArgumentParser()
parser.add_argument('--prmtop', default=str(BASE / 'data/structures/quaternary_mvp/quaternary_mvp.prmtop'))
parser.add_argument('--pdb', default=str(BASE / 'data/structures/quaternary_mvp/quaternary_mvp_minimized.pdb'))
parser.add_argument('--name', default='quaternary_mvp_rep1')
parser.add_argument('--outdir', default=str(BASE / 'data/md_runs/quaternary_mvp/rep1'))
parser.add_argument('--prod-ns', type=float, default=50)
parser.add_argument('--platform', default='CUDA')
parser.add_argument('--gpu', default='2')
parser.add_argument('--seed', type=int, default=20260503)
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

logfile = outdir / f"{args.name}.log"
chkfile = outdir / f"{args.name}.chk"
dcdfile = outdir / f"{args.name}.dcd"

def log(msg):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line)
    with open(logfile, 'a') as f:
        f.write(line + '\n')

log(f"Starting {args.name} on GPU {args.gpu}")
log(f"Production: {args.prod_ns} ns")

prmtop = app.AmberPrmtopFile(args.prmtop)
pdb = app.PDBFile(args.pdb)

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds,
)

integrator = mm.LangevinMiddleIntegrator(
    300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond
)
integrator.setRandomNumberSeed(args.seed)

platform = mm.Platform.getPlatformByName(args.platform)
properties = {'CudaDeviceIndex': args.gpu, 'CudaPrecision': 'mixed'}

simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# Derive inpcrd path from prmtop path (same basename, .inpcrd extension)
inpcrd_path = str(Path(args.prmtop).with_suffix('.inpcrd'))
inpcrd = app.AmberInpcrdFile(inpcrd_path)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# reporters
dcd_reporter = app.DCDReporter(str(dcdfile), 5000)  # 10 ps interval
state_reporter = app.StateDataReporter(
    str(outdir / f"{args.name}_state.log"), 5000,
    step=True, time=True, potentialEnergy=True, temperature=True,
    speed=True, remainingTime=True, totalSteps=int(args.prod_ns * 1000 / 0.002)
)

simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)

log(f"Initial energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")

# Short NVT equilibration (100 ps)
log("NVT equilibration 100 ps...")
simulation.step(50000)

# NPT production
n_steps = int(args.prod_ns * 1000 / 0.002)
log(f"Production run: {n_steps} steps ({args.prod_ns} ns)...")
simulation.step(n_steps)

# Save final state
simulation.saveCheckpoint(str(chkfile))
final_state = simulation.context.getState(getPositions=True)
with open(outdir / f"{args.name}_final.pdb", 'w') as f:
    app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)

log("Done!")
