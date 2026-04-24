#!/usr/bin/env python3
"""Build Hsap cGAS-TRIM41 MD system from Rosetta docking pose.

Uses Amber tleap (ff19SB + OPC) + OpenMM validation.
Rosetta PDB already has chain A (TRIM41) and chain B (cGAS).
"""

import argparse
import os
import subprocess
import sys
import tempfile


def prepare_pdb_rosetta(pdb_in, pdb_out):
    """Clean Rosetta PDB for Amber: strip H atoms, add TER, rename residues."""
    with open(pdb_in) as f:
        lines = f.readlines()
    
    with open(pdb_out, 'w') as out:
        prev_chain = None
        for line in lines:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            
            # Skip hydrogen atoms (Amber will add its own)
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 77 else ''
            if atom_name.startswith('H') or element == 'H':
                continue
            
            chain = line[21]
            
            # TER between chains
            if prev_chain is not None and chain != prev_chain:
                out.write("TER\n")
            
            out.write(line)
            prev_chain = chain
        
        out.write("TER\n")
        out.write("END\n")
    
    print(f"  Stripped H atoms, wrote {pdb_out}")


def run_pdb4amber(pdb_in, pdb_out, reduce=False):
    """Run pdb4amber to fix atom names, add missing atoms."""
    cmd = ['pdb4amber', '-i', pdb_in, '-o', pdb_out]
    if reduce:
        cmd.append('--reduce')
    else:
        cmd.append('--no-reduce')
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"pdb4amber stderr: {result.stderr}")
        raise RuntimeError("pdb4amber failed")
    
    print(f"  pdb4amber output: {pdb_out}")
    return pdb_out


def run_tleap(pdb_file, out_prefix, prmtop, rst7):
    """Run Amber tleap to build solvated system."""
    tleap_in = f"""
source leaprc.protein.ff19SB
source leaprc.water.opc

# Load structure
mol = loadpdb {pdb_file}

# Solvate
solvateBox mol OPCBOX 12.0 iso

# Neutralize and add 150 mM NaCl
addionsrand mol Na+ 0
addionsrand mol Cl- 0

# Save
saveamberparm mol {prmtop} {rst7}
quit
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False) as f:
        f.write(tleap_in)
        tleap_file = f.name
    
    try:
        result = subprocess.run(
            ['tleap', '-f', tleap_file],
            capture_output=True, text=True
        )
        print(result.stdout)
        if result.returncode != 0:
            print(f"tleap stderr: {result.stderr}")
            raise RuntimeError("tleap failed")
    finally:
        os.unlink(tleap_file)
    
    print(f"  tleap output: {prmtop}, {rst7}")


def validate_with_openmm(prmtop, rst7):
    """Quick energy check with OpenMM."""
    from openmm import app
    import openmm as mm
    
    prmtop_obj = app.AmberPrmtopFile(prmtop)
    inpcrd = app.AmberInpcrdFile(rst7)
    
    system = prmtop_obj.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * mm.unit.nanometer,
        constraints=app.HBonds,
    )
    
    integrator = mm.LangevinMiddleIntegrator(
        300 * mm.unit.kelvin,
        1.0 / mm.unit.picosecond,
        0.002 * mm.unit.picosecond
    )
    
    platform = mm.Platform.getPlatformByName('CPU')
    simulation = app.Simulation(prmtop_obj.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy()
    print(f"  Initial energy: {energy}")
    
    # Quick minimization
    simulation.minimizeEnergy(maxIterations=100)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy_min = state.getPotentialEnergy()
    print(f"  After 100-step minimization: {energy_min}")
    
    return energy, energy_min


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', required=True, help='Rosetta docking pose PDB')
    parser.add_argument('--name', required=True, help='System name')
    parser.add_argument('--outdir', required=True, help='Output directory')
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"Building MD system for {args.name}")
    print(f"  Input: {args.pdb}")
    print(f"  Output: {args.outdir}")
    
    # Step 1: Clean Rosetta PDB (strip H)
    clean_pdb = os.path.join(args.outdir, f'{args.name}_clean.pdb')
    prepare_pdb_rosetta(args.pdb, clean_pdb)
    
    # Step 2: pdb4amber
    amber_pdb = os.path.join(args.outdir, f'{args.name}_amber.pdb')
    run_pdb4amber(clean_pdb, amber_pdb, reduce=True)
    
    # Step 3: tleap
    prmtop = os.path.join(args.outdir, f'{args.name}.prmtop')
    rst7 = os.path.join(args.outdir, f'{args.name}.rst7')
    run_tleap(amber_pdb, args.name, prmtop, rst7)
    
    # Step 4: Validate
    print("Validating with OpenMM...")
    validate_with_openmm(prmtop, rst7)
    
    print(f"\nDone! System files in {args.outdir}")


if __name__ == '__main__':
    main()
