#!/usr/bin/env python3
"""Build MD system with configurable solvation box and conservative EM.

Improvements over v1:
- Supports truncated octahedron (--oct) for ~30% water savings
- Configurable buffer distance (--buffer)
- Conservative EM with position restraints on protein backbone
- Better handling of Rosetta poses (no PDBFixer needed for side chains)

Usage:
    python build_system_v2.py --pdb pose.pdb --name Hsap_WT --outdir data/md_runs/Hsap_WT \
        --oct --buffer 10.0 --platform CUDA
"""
import argparse
import os
import subprocess
import sys
import shutil


def prepare_pdb_rosetta(pdb_in, pdb_out):
    """Clean Rosetta PDB for Amber: strip H atoms, add TER between chains."""
    with open(pdb_in) as f:
        lines = f.readlines()
    
    with open(pdb_out, 'w') as out:
        prev_chain = None
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21]
                # Strip H atoms
                elem = line[76:78].strip()
                if elem == 'H':
                    continue
                # Add TER when chain changes
                if prev_chain is not None and chain != prev_chain:
                    out.write("TER\n")
                out.write(line)
                prev_chain = chain
            elif line.startswith('TER'):
                out.write(line)
        out.write("TER\n")
        out.write("END\n")
    
    # Count
    with open(pdb_out) as f:
        n_atoms = sum(1 for l in f if l.startswith('ATOM') or l.startswith('HETATM'))
    print(f"  Cleaned PDB: {n_atoms} non-H atoms")
    return pdb_out


def run_tleap(pdb_path, outdir, name, buffer=10.0, use_oct=False):
    """Run Amber tleap to build solvated system."""
    
    prmtop = os.path.join(outdir, f"{name}.prmtop")
    rst7 = os.path.join(outdir, f"{name}.rst7")
    
    if use_oct:
        solvate_cmd = f"solvateOct complex OPC {buffer}"
        box_type = "truncated octahedron"
    else:
        solvate_cmd = f"solvateBox complex OPC {buffer}"
        box_type = "rectangular"
    
    leap_script = f"""source leaprc.protein.ff19SB
source leaprc.water.opc
complex = loadPdb {pdb_path}
check complex
{solvate_cmd}
addIonsRand complex Na+ 0
addIonsRand complex Cl- 0
saveAmberParm complex {prmtop} {rst7}
charge complex
quit
"""
    
    leap_input = os.path.join(outdir, "tleap.in")
    with open(leap_input, 'w') as f:
        f.write(leap_script)
    
    print(f"\n[tleap] Building {box_type} box, buffer={buffer}Å")
    
    tleap_path = shutil.which("tleap") or os.path.expanduser("~/miniforge3/envs/cgas-md/bin/tleap")
    cmd = [tleap_path, "-f", leap_input]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    log_path = os.path.join(outdir, "tleap.log")
    with open(log_path, 'w') as f:
        f.write(result.stdout)
        f.write(result.stderr)
    
    if result.returncode != 0:
        print(f"  ❌ tleap failed")
        return None, None
    
    if not os.path.exists(prmtop) or not os.path.exists(rst7):
        print(f"  ❌ Output files missing")
        return None, None
    
    print(f"  ✅ tleap complete")
    return prmtop, rst7


def conservative_em(prmtop, rst7, outdir, name, max_steps=5000, platform_name='CUDA'):
    """EM with position restraints on protein backbone to prevent unfolding."""
    from openmm import app
    import openmm as mm
    
    print(f"\n[EM] Conservative minimization ({max_steps} steps, backbone restraints)...")
    
    prmtop_obj = app.AmberPrmtopFile(prmtop)
    inpcrd = app.AmberInpcrdFile(rst7)
    
    system = prmtop_obj.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*mm.unit.nanometer,
        constraints=None,
    )
    
    # Add harmonic restraints to backbone atoms (CA, C, N, O)
    force = mm.CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", 100.0 * mm.unit.kilojoule_per_mole / mm.unit.nanometer**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    
    for atom in prmtop_obj.topology.atoms():
        if atom.name in ('CA', 'C', 'N', 'O') and atom.residue.name not in ('WAT', 'Na+', 'Cl-'):
            idx = atom.index
            pos = inpcrd.positions[idx]
            force.addParticle(idx, [pos.x, pos.y, pos.z])
    
    system.addForce(force)
    print(f"  Restrained {force.getNumParticles()} backbone atoms")
    
    integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
    
    platform = mm.Platform.getPlatformByName(platform_name)
    props = {'CudaPrecision': 'mixed'} if platform_name == 'CUDA' else {}
    context = mm.Context(system, integrator, platform, props)
    context.setPositions(inpcrd.positions)
    
    e_init = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
    print(f"  Initial: {e_init:.1f} kJ/mol")
    
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=max_steps)
    
    state = context.getState(getEnergy=True, getPositions=True)
    e_final = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
    print(f"  Final:   {e_final:.1f} kJ/mol")
    print(f"  ΔE:      {e_final - e_init:.1f}")
    
    minimized_pdb = os.path.join(outdir, f"{name}_minimized.pdb")
    app.PDBFile.writeFile(prmtop_obj.topology, state.getPositions(), open(minimized_pdb, 'w'))
    
    del context
    print(f"  ✅ Saved: {minimized_pdb}")
    return minimized_pdb


def validate_system(prmtop, rst7):
    """Quick validation: load and report atom counts."""
    from parmed.amber import AmberParm
    p = AmberParm(prmtop)
    n_wat = sum(1 for r in p.residues if r.name == 'WAT')
    n_atoms = len(p.atoms)
    print(f"\n[System] {n_atoms} atoms, {n_wat} waters")
    return n_atoms, n_wat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--buffer', type=float, default=10.0)
    parser.add_argument('--oct', action='store_true', help='Use truncated octahedron')
    parser.add_argument('--em-steps', type=int, default=5000)
    parser.add_argument('--platform', default='CUDA')
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"{'='*60}")
    print(f"  MD System Builder v2")
    print(f"{'='*60}")
    print(f"  Input:  {args.pdb}")
    print(f"  Name:   {args.name}")
    print(f"  Box:    {'octahedron' if args.oct else 'rectangular'}, {args.buffer}Å buffer")
    print()
    
    # Step 1: Clean PDB
    clean_pdb = os.path.join(args.outdir, f"{args.name}_clean.pdb")
    prepare_pdb_rosetta(args.pdb, clean_pdb)
    
    # Step 2: tleap
    prmtop, rst7 = run_tleap(clean_pdb, args.outdir, args.name, 
                             buffer=args.buffer, use_oct=args.oct)
    if prmtop is None:
        sys.exit(1)
    
    # Step 3: Validate
    validate_system(prmtop, rst7)
    
    # Step 4: Conservative EM
    minimized = conservative_em(prmtop, rst7, args.outdir, args.name, 
                                max_steps=args.em_steps, platform_name=args.platform)
    
    # Step 5: Validate minimized structure dimensions
    from MDAnalysis import Universe
    u = Universe(minimized)
    ca = u.select_atoms('protein and name CA')
    dims = ca.positions.max(axis=0) - ca.positions.min(axis=0)
    print(f"\n[Minimized] Protein dimensions: {dims[0]:.1f} x {dims[1]:.1f} x {dims[2]:.1f} Å")
    
    print(f"\n{'='*60}")
    print(f"  ✅ System ready!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
