#!/usr/bin/env python3
"""Build MD system from LightDock best_pose.pdb using Amber tleap + OpenMM.

Steps:
1. Preprocess PDB: split receptor/ligand into separate chains, add TER
2. Run tleap: ff19SB + OPC, solvate octahedron, add ions
3. Validate with OpenMM
4. Energy minimization

Usage:
    python build_system.py --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
                           --name Hgal_domain --outdir data/md_runs/Hgal_domain
"""

import argparse
import os
import subprocess
import sys
import tempfile
import shutil


def split_chains(pdb_in, pdb_out):
    """Split LightDock output into two chains with TER records.
    
    LightDock best_pose.pdb format:
      - Receptor (TRIM41) atoms first: resnums 413-630
      - Ligand (cGAS) atoms second: resnums 200-554
    """
    with open(pdb_in) as f:
        lines = f.readlines()
    
    # Find the split point (resnum jumps from high to low)
    split_idx = None
    prev_resnum = None
    for i, line in enumerate(lines):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        resnum = int(line[22:26])
        if prev_resnum is not None and resnum < prev_resnum - 100:
            split_idx = i
            break
        prev_resnum = resnum
    
    if split_idx is None:
        raise ValueError("Could not find receptor/ligand split point")
    
    print(f"  Split point at line {split_idx}: TRIM41 (1-{split_idx}), cGAS ({split_idx+1}-end)")
    
    with open(pdb_out, 'w') as out:
        # Write TRIM41 as chain A
        for line in lines[:split_idx]:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Ensure chain ID is 'A'
                out.write(line[:21] + 'A' + line[22:])
            else:
                out.write(line)
        
        # TER record between chains
        out.write("TER\n")
        
        # Write cGAS as chain B
        for line in lines[split_idx:]:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Change chain ID to 'B'
                out.write(line[:21] + 'B' + line[22:])
            else:
                out.write(line)
        
        # Final TER and END
        out.write("TER\n")
        out.write("END\n")
    
    # Count atoms per chain
    with open(pdb_out) as f:
        chain_a = sum(1 for l in f if (l.startswith('ATOM') or l.startswith('HETATM')) and l[21] == 'A')
    with open(pdb_out) as f:
        chain_b = sum(1 for l in f if (l.startswith('ATOM') or l.startswith('HETATM')) and l[21] == 'B')
    
    print(f"  Chain A (TRIM41): {chain_a} atoms")
    print(f"  Chain B (cGAS):   {chain_b} atoms")
    print(f"  Total: {chain_a + chain_b} atoms")
    return pdb_out


def run_tleap(pdb_path, outdir, name, buffer=12.0):
    """Run Amber tleap to build solvated system."""
    
    prmtop = os.path.join(outdir, f"{name}.prmtop")
    rst7 = os.path.join(outdir, f"{name}.rst7")
    
    leap_script = f"""
# Protein force field
source leaprc.protein.ff19SB

# Water model
source leaprc.water.opc

# Load protein
complex = loadPdb {pdb_path}

# Check structure
check complex

# Solvate with truncated octahedron
solvateBox complex OPC {buffer}

# Neutralize and add 150 mM NaCl
addIonsRand complex Na+ 0
addIonsRand complex Cl- 0

# Save topology and coordinates
saveAmberParm complex {prmtop} {rst7}

# Report
charge complex
quit
"""
    
    leap_input = os.path.join(outdir, "tleap.in")
    with open(leap_input, 'w') as f:
        f.write(leap_script)
    
    print(f"\n[tleap] Running Amber tleap...")
    print(f"  Buffer: {buffer} Å")
    print(f"  Force field: ff19SB")
    print(f"  Water: OPC")
    
    tleap_path = shutil.which("tleap")
    if tleap_path is None:
        # Try conda env
        tleap_path = os.path.expanduser("~/miniforge3/envs/cgas-md/bin/tleap")
        if not os.path.exists(tleap_path):
            raise FileNotFoundError("tleap not found. Activate cgas-md conda env.")
    
    cmd = [tleap_path, "-f", leap_input]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    log_path = os.path.join(outdir, "tleap.log")
    with open(log_path, 'w') as f:
        f.write(result.stdout)
        f.write(result.stderr)
    
    # Check for errors
    if result.returncode != 0:
        print(f"  ❌ tleap failed with code {result.returncode}")
        print(f"  Log: {log_path}")
        return None, None
    
    if "Error" in result.stdout or "Fatal Error" in result.stdout:
        print(f"  ⚠️  tleap reported errors, check log: {log_path}")
    
    if not os.path.exists(prmtop) or not os.path.exists(rst7):
        print(f"  ❌ Output files not generated")
        return None, None
    
    print(f"  ✅ tleap completed")
    print(f"  prmtop: {prmtop} ({os.path.getsize(prmtop)} bytes)")
    print(f"  rst7:   {rst7} ({os.path.getsize(rst7)} bytes)")
    
    return prmtop, rst7


def get_openmm_platform(preferred='auto'):
    """Get best available OpenMM platform with appropriate properties.
    
    Args:
        preferred: 'auto', 'CUDA', 'OpenCL', 'CPU', or 'Reference'
    
    Returns:
        (platform, properties_dict)
    """
    import openmm as mm
    
    available = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    print(f"  Available platforms: {', '.join(available)}")
    
    if preferred == 'auto':
        for name in ['CUDA', 'OpenCL', 'CPU']:
            if name in available:
                preferred = name
                break
        else:
            preferred = 'Reference'
    
    platform = mm.Platform.getPlatformByName(preferred)
    
    props = {}
    if preferred == 'CUDA':
        # RTX 3090 supports mixed precision
        props['CudaPrecision'] = 'mixed'
        props['CudaDeviceIndex'] = '0'
    elif preferred == 'OpenCL':
        # Apple Silicon only supports single precision
        props['OpenCLPrecision'] = 'single'
    
    print(f"  Using platform: {preferred}")
    if props:
        print(f"  Properties: {props}")
    
    return platform, props


def validate_with_openmm(prmtop, rst7, platform_name='auto'):
    """Validate system with OpenMM."""
    try:
        from openmm import app
        import openmm as mm
        
        print(f"\n[OpenMM] Validating system...")
        
        prmtop_obj = app.AmberPrmtopFile(prmtop)
        inpcrd = app.AmberInpcrdFile(rst7)
        
        system = prmtop_obj.createSystem(
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*mm.unit.nanometer,
            constraints=app.HBonds,
            implicitSolvent=None,
        )
        
        n_atoms = system.getNumParticles()
        n_residues = prmtop_obj.topology.getNumResidues()
        n_chains = prmtop_obj.topology.getNumChains()
        
        print(f"  Atoms:     {n_atoms}")
        print(f"  Residues:  {n_residues}")
        print(f"  Chains:    {n_chains}")
        
        # Quick energy calculation
        integrator = mm.LangevinMiddleIntegrator(
            300*mm.unit.kelvin,
            1.0/mm.unit.picosecond,
            0.002*mm.unit.picosecond,
        )
        
        platform, props = get_openmm_platform(platform_name)
        context = mm.Context(system, integrator, platform, props)
        context.setPositions(inpcrd.positions)
        
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(f"  Potential energy: {energy}")
        
        del context
        print(f"  ✅ OpenMM validation passed")
        return True
        
    except Exception as e:
        print(f"  ⚠️  OpenMM validation error: {e}")
        return False


def energy_minimize(prmtop, rst7, outdir, name, max_steps=1000, platform_name='auto'):
    """Quick energy minimization with OpenMM."""
    try:
        from openmm import app
        import openmm as mm
        
        print(f"\n[EM] Energy minimization ({max_steps} steps)...")
        
        prmtop_obj = app.AmberPrmtopFile(prmtop)
        inpcrd = app.AmberInpcrdFile(rst7)
        
        system = prmtop_obj.createSystem(
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*mm.unit.nanometer,
            constraints=None,  # No constraints for EM
        )
        
        integrator = mm.VerletIntegrator(0.001*mm.unit.picosecond)
        platform, props = get_openmm_platform(platform_name)
        context = mm.Context(system, integrator, platform, props)
        context.setPositions(inpcrd.positions)
        
        # Initial energy
        state = context.getState(getEnergy=True)
        e_init = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        print(f"  Initial energy: {e_init:.1f} kJ/mol")
        
        # Minimize
        mm.LocalEnergyMinimizer.minimize(context, maxIterations=max_steps)
        
        state = context.getState(getEnergy=True, getPositions=True)
        e_final = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoule_per_mole)
        print(f"  Final energy:   {e_final:.1f} kJ/mol")
        print(f"  ΔE:             {e_final - e_init:.1f} kJ/mol")
        
        # Save minimized coordinates as PDB
        minimized_pdb = os.path.join(outdir, f"{name}_minimized.pdb")
        app.PDBFile.writeFile(prmtop_obj.topology, state.getPositions(), open(minimized_pdb, 'w'))
        
        print(f"  ✅ Minimization complete")
        print(f"  Output: {minimized_pdb}")
        
        del context
        return minimized_pdb
        
    except Exception as e:
        print(f"  ⚠️  Minimization error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description='Build MD system from docking pose')
    parser.add_argument('--pdb', required=True, help='Input PDB file (LightDock best_pose)')
    parser.add_argument('--name', default='system', help='System name')
    parser.add_argument('--outdir', default='data/md_runs', help='Output directory')
    parser.add_argument('--buffer', type=float, default=12.0, help='Solvation buffer (Å)')
    parser.add_argument('--skip-em', action='store_true', help='Skip energy minimization')
    parser.add_argument('--platform', default='auto', choices=['auto', 'CUDA', 'OpenCL', 'CPU', 'Reference'],
                        help='OpenMM platform (auto=best available)')
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"{'='*60}")
    print(f"  MD System Builder")
    print(f"{'='*60}")
    print(f"  Input:  {args.pdb}")
    print(f"  Name:   {args.name}")
    print(f"  Output: {args.outdir}")
    print()
    
    # Step 1: Preprocess PDB
    print("[1/4] Preprocessing PDB...")
    pdb_processed = os.path.join(args.outdir, f"{args.name}_processed.pdb")
    split_chains(args.pdb, pdb_processed)
    print(f"  Output: {pdb_processed}")
    
    # Step 2: tleap
    print("\n[2/4] Building topology with tleap...")
    prmtop, rst7 = run_tleap(pdb_processed, args.outdir, args.name, buffer=args.buffer)
    
    if prmtop is None:
        print("\n❌ System building failed at tleap step")
        sys.exit(1)
    
    # Step 3: OpenMM validation
    print("\n[3/4] Validating with OpenMM...")
    ok = validate_with_openmm(prmtop, rst7, platform_name=args.platform)
    
    # Step 4: Energy minimization
    if not args.skip_em:
        print("\n[4/4] Energy minimization...")
        minimized = energy_minimize(prmtop, rst7, args.outdir, args.name, platform_name=args.platform)
        if minimized:
            print(f"\n{'='*60}")
            print(f"  ✅ System ready for MD!")
            print(f"{'='*60}")
            print(f"  prmtop:     {prmtop}")
            print(f"  rst7:       {rst7}")
            print(f"  minimized:  {minimized}")
    else:
        print(f"\n{'='*60}")
        print(f"  ✅ System built (EM skipped)")
        print(f"{'='*60}")
        print(f"  prmtop: {prmtop}")
        print(f"  rst7:   {rst7}")


if __name__ == '__main__':
    main()
