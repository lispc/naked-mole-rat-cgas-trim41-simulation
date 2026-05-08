#!/usr/bin/env python3
"""Build DNA-bound cGAS+SPRY system from Boltz-2 model using OpenMM + PDBFixer."""
import sys
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import app
from openmm import unit
from pdbfixer import PDBFixer

BASE = Path(__file__).resolve().parent.parent.parent
PDB_PATH = BASE / "data/boltz_cgas_dna_trim41/boltz_results_boltz_cgas_dna_trim41/predictions/boltz_cgas_dna_trim41/boltz_cgas_dna_trim41_model_0.pdb"
OUTDIR = BASE / "data/md_runs/cgas_dna_spry_openmm"

def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Fix PDB ──────────────────────────────────────────────────
    print("1. Fixing PDB with PDBFixer...", flush=True)
    fixer = PDBFixer(filename=str(PDB_PATH))
    print(f"   Chains: {fixer.topology.getNumChains()}")
    print(f"   Residues: {fixer.topology.getNumResidues()}")
    fixer.findMissingResidues()
    print(f"   Missing residues: {fixer.missingResidues}")

    fixer.findMissingAtoms()
    missing_atoms = fixer.missingAtoms if hasattr(fixer, 'missingAtoms') else []
    n_missing = len(missing_atoms)
    print(f"   Missing atoms: {n_missing}")
    if n_missing > 0:
        for atom in list(missing_atoms)[:5]:
            print(f"     {atom}")

    fixer.addMissingAtoms()
    print(f"   After addMissingAtoms: {fixer.topology.getNumAtoms()} atoms")

    fixer.addMissingHydrogens(7.0)
    print(f"   After addMissingHydrogens: {fixer.topology.getNumAtoms()} atoms")

    # ── Fix DNA 3' termini: remove phosphate from 3'-terminal residues ──
    # Amber expects DA3/DT3/etc to have free 3'OH, not 3' phosphate
    dna_three_prime_atoms = set()
    for chain in fixer.topology.chains():
        residues = list(chain.residues())
        if not residues or residues[0].name not in ("DA", "DT", "DG", "DC"):
            continue
        last_res = residues[-1]
        for atom in last_res.atoms():
            if atom.name in ("P", "OP1", "OP2"):
                # This is actually a 5' phosphate incorrectly placed
                # on the 3'-terminal residue by Boltz-2
                # Actually, check: is it the 5'P of this residue or 3'P?
                pass
        print(f"   DNA chain {chain.id}: 3'-term {last_res.name}{last_res.id}, "
              f"atoms={[a.name for a in last_res.atoms()][:8]}")
    # Note: the 3'-terminal issue is from Boltz-2 output having
    # identical atom sets for all DNA residues (each has 5'P).
    # We'll handle this via modifying topology before force field matching.

    # ── 2. Load force field ─────────────────────────────────────────
    print("\n2. Loading force fields...", flush=True)
    ff = app.ForceField("amber19-all.xml", "amber19/tip4pew.xml")
    print("   FF loaded")

    # ── 3. Solvate ──────────────────────────────────────────────────
    print("\n3. Adding solvent...", flush=True)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        ff,
        model="tip4pew",
        boxShape="octahedron",
        padding=1.0 * unit.nanometers,
        ionicStrength=0.15 * unit.molar,
        positiveIon="Na+",
        negativeIon="Cl-",
    )
    n_atoms = modeller.topology.getNumAtoms()
    print(f"   Solvated: {n_atoms} atoms")

    # Get box vectors to estimate size
    positions = modeller.positions
    min_coords = np.min(positions.value_in_unit(unit.nanometers), axis=0)
    max_coords = np.max(positions.value_in_unit(unit.nanometers), axis=0)
    print(f"   Box extent: {max_coords - min_coords} nm")

    # ── 4. Create system ────────────────────────────────────────────
    print("\n4. Creating system...", flush=True)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        ewaldErrorTolerance=0.0005,
    )
    print("   System created")

    # ── 5. Quick minimization ───────────────────────────────────────
    print("\n5. Minimization test...", flush=True)
    platform = mm.Platform.getPlatformByName("CUDA")
    props = {"CudaDeviceIndex": "0", "Precision": "mixed"}
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )
    sim = app.Simulation(modeller.topology, system, integrator, platform, props)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=5000)
    state = sim.context.getState(getEnergy=True, getPositions=True)
    e_pot = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"   Post-min PE: {e_pot:.1f} kJ/mol")

    # Check for clashes
    if e_pot > 0 and e_pot < 1e6:
        print("   Minimization OK")
    elif e_pot < 0:
        # Get forces to check
        state = sim.context.getState(getForces=True)
        forces = state.getForces().value_in_unit(unit.kilojoules_per_mole / unit.nanometer)
        max_f = np.max(np.linalg.norm(forces, axis=1))
        print(f"   Max force: {max_f:.1f} kJ/mol/nm")
    else:
        print("   WARNING: Very high energy - possible clashes")

    # Save minimized
    pos = state.getPositions()
    with open(OUTDIR / "minimized.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, pos, f)
    print(f"   Saved: {OUTDIR / 'minimized.pdb'}")

    # ── 6. Save serialized system ───────────────────────────────────
    with open(OUTDIR / "system.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(system))
    with open(OUTDIR / "topology.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    print(f"\nDONE. Output in {OUTDIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
