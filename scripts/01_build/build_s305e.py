#!/usr/bin/env python3
"""
Build S305E mutant system from Hsap_WT structure.

S305E is a phosphomimetic mutation: Ser -> Glu at residue 305.
This introduces a -1 charge to mimic the electrostatic effect of phosphorylation.

Usage:
    python build_s305e_system.py \
        --pdb data/md_runs/Hsap_WT/Hsap_WT_amber.pdb \
        --name Hsap_WT_S305E \
        --outdir data/md_runs/Hsap_WT_S305E
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def get_conda_bin():
    """Get the conda environment bin directory from sys.executable."""
    return Path(sys.executable).parent


def run_cmd(cmd, cwd=None, env=None):
    """Run shell command and print output."""
    print(f"  $ {cmd}")
    import os
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, env=run_env)
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"Command failed: {cmd}")
    return result.stdout


def build_s305e(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print(f"Building S305E system: {args.name}")
    print(f"Input PDB: {args.pdb}")
    print(f"Output dir: {outdir}")
    print("=" * 60)
    
    # Step 1: Apply S305E mutation manually
    print("\n[Step 1] Applying S305E mutation (SER B 324 -> GLU)...")
    
    mutated_pdb = outdir / f"{args.name}_mutated.pdb"
    run_cmd(
        f"{sys.executable} scripts/mutate_s305e.py "
        f"--input {args.pdb} "
        f"--output {mutated_pdb} "
        f"--chain B --resseq 324 --old SER --new GLU"
    )
    
    # Verify mutation
    print("\n  Verifying mutation...")
    result = run_cmd(f"grep 'GLU B 324' {mutated_pdb} | head -3")
    print(f"  {result.strip()}")
    
    # Step 2: pdb4amber --reduce
    print("\n[Step 2] Running pdb4amber with reduce...")
    amber_pdb = outdir / f"{args.name}_amber.pdb"
    conda_bin = get_conda_bin()
    run_cmd(
        f"{conda_bin}/pdb4amber -i {mutated_pdb} -o {amber_pdb} --reduce",
        env={"PATH": f"{conda_bin}:{os.environ.get('PATH', '')}"}
    )
    
    # Step 3: tleap build system
    print("\n[Step 3] Running tleap...")
    tleap_in = outdir / f"tleap_{args.name}.in"
    prmtop = outdir / f"{args.name}.prmtop"
    rst7 = outdir / f"{args.name}.rst7"
    solvated_pdb = outdir / f"{args.name}_solvated.pdb"
    
    tleap_script = f"""source leaprc.protein.ff19SB
source leaprc.water.opc

complex = loadPdb {amber_pdb}
check complex
solvateOct complex OPC 12.0
addIonsRand complex Na+ 0
addIonsRand complex Cl- 0
saveAmberParm complex {prmtop} {rst7}
savePdb complex {solvated_pdb}
charge complex
quit
"""
    
    with open(tleap_in, 'w') as f:
        f.write(tleap_script)
    
    tleap_log = outdir / "tleap.log"
    run_cmd(f"{conda_bin}/tleap -f {tleap_in} > {tleap_log} 2>&1")
    
    # Parse tleap results
    with open(tleap_log) as f:
        log = f.read()
    
    # Find charge and atom count
    for line in log.split('\n'):
        if 'Total atoms' in line or 'Total unperturbed' in line:
            print(f"  {line.strip()}")
        if 'Charge on' in line:
            print(f"  {line.strip()}")
    
    # Step 4: Energy minimization (crucial for mutated systems)
    print("\n[Step 4] Energy minimization with OpenMM...")
    minimized_pdb = outdir / f"{args.name}_minimized.pdb"
    
    from openmm import app, unit, CustomExternalForce, LangevinIntegrator, LocalEnergyMinimizer, Platform
    
    prmtop_obj = app.AmberPrmtopFile(str(prmtop))
    pdb_obj = app.PDBFile(str(solvated_pdb))
    system = prmtop_obj.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*app.unit.nanometer,
        constraints=app.HBonds
    )
    
    # Add weak position restraints to protein backbone
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", 100.0)  # kJ/mol/nm^2
    for atom in pdb_obj.topology.atoms():
        if atom.residue.chain.id in ['A', 'B'] and atom.name in ['N', 'CA', 'C', 'O']:
            force.addParticle(atom.index, pdb_obj.positions[atom.index])
    system.addForce(force)
    
    integrator = LangevinIntegrator(300*app.unit.kelvin, 1.0/app.unit.picosecond, 0.002*app.unit.picosecond)
    platform = Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(prmtop_obj.topology, system, integrator, platform)
    simulation.context.setPositions(pdb_obj.positions)
    
    # Initial energy
    state = simulation.context.getState(getEnergy=True)
    init_e = state.getPotentialEnergy().value_in_unit(app.unit.kilojoule_per_mole)
    print(f"  Initial energy: {init_e:.1f} kJ/mol")
    
    # Minimize
    print("  Minimizing...")
    LocalEnergyMinimizer.minimize(simulation.context, maxIterations=5000)
    
    # Final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_e = state.getPotentialEnergy().value_in_unit(app.unit.kilojoule_per_mole)
    print(f"  Final energy: {final_e:.1f} kJ/mol")
    print(f"  Delta: {final_e - init_e:.1f} kJ/mol")
    
    # Save
    with open(minimized_pdb, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    print("  Minimization complete")
    
    print(f"\n{'=' * 60}")
    print("Build complete!")
    print(f"  prmtop:  {prmtop}")
    print(f"  rst7:    {rst7}")
    print(f"  minimized PDB: {minimized_pdb}")
    print(f"{'=' * 60}")
    
    return {
        'prmtop': prmtop,
        'rst7': rst7,
        'minimized': minimized_pdb,
        'solvated': solvated_pdb,
    }


def main():
    parser = argparse.ArgumentParser(description='Build S305E mutant system')
    parser.add_argument('--pdb', default='data/md_runs/Hsap_WT/Hsap_WT_amber.pdb',
                        help='Input protein-only PDB (pdb4amber-processed)')
    parser.add_argument('--name', default='Hsap_WT_S305E',
                        help='System name')
    parser.add_argument('--outdir', default='data/md_runs/Hsap_WT_S305E',
                        help='Output directory')
    args = parser.parse_args()
    
    build_s305e(args)


if __name__ == '__main__':
    main()
