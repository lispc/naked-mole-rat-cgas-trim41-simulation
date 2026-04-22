#!/usr/bin/env python3
"""
Build MD simulation system from AF3 PDB using AmberTools/OpenMM.
Steps:
  1. Parse AF3 PDB, extract chains, renumber if needed
  2. pdb4amber cleanup
  3. tleap: load ff19SB + OPC, solvate, add ions
  4. Convert to OpenMM format
  5. Optional: domain truncation based on interface analysis

Usage:
  python scripts/build_system.py --pdb structures/af3_raw/job1_Hsap_WT/ranked_0.pdb \
      --name Hsap_WT --truncate  # with truncation
"""
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path


def run_cmd(cmd, cwd=None):
    print(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def prepare_pdb(input_pdb: Path, output_pdb: Path, cgas_chain="A", trim41_chain="B"):
    """Extract cGAS and TRIM41 chains, renumber, add TER records."""
    with open(input_pdb) as f:
        lines = f.readlines()
    
    chains = {}
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]
            chains.setdefault(chain, []).append(line)
    
    print(f"Chains found: {list(chains.keys())}")
    
    with open(output_pdb, "w") as out:
        # Write cGAS
        for line in chains.get(cgas_chain, []):
            out.write(line)
        out.write("TER\n")
        # Write TRIM41
        for line in chains.get(trim41_chain, []):
            out.write(line)
        out.write("TER\n")
        out.write("END\n")
    
    print(f"Prepared PDB: {output_pdb}")


def pdb4amber(input_pdb: Path, output_pdb: Path):
    """Run pdb4amber to fix PDB issues."""
    cmd = ["pdb4amber", "-i", str(input_pdb), "-o", str(output_pdb), "--reduce"]
    run_cmd(cmd)
    print(f"pdb4amber output: {output_pdb}")


def tleap_build(input_pdb: Path, out_dir: Path, name: str, buffer=12.0, ion_conc=0.15):
    """Run tleap to build solvated system."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    prmtop = out_dir / f"{name}.prmtop"
    inpcrd = out_dir / f"{name}.inpcrd"
    pdb_out = out_dir / f"{name}_solv.pdb"
    
    leap_script = f"""
source leaprc.protein.ff19SB
source leaprc.water.opc
loadAmberParams frcmod.ionsjc_tip4pew

mol = loadpdb {input_pdb}
check mol

solvatebox mol OPCBOX {buffer}
addIonsRand mol Na+ 0
addIonsRand mol Cl- 0
addIonsRand mol Na+ {ion_conc}
addIonsRand mol Cl- {ion_conc}

saveamberparm mol {prmtop} {inpcrd}
savepdb mol {pdb_out}
quit
"""
    leap_in = out_dir / "leap.in"
    with open(leap_in, "w") as f:
        f.write(leap_script)
    
    cmd = ["tleap", "-f", str(leap_in)]
    result = run_cmd(cmd)
    
    # Check output
    if not prmtop.exists() or not inpcrd.exists():
        raise RuntimeError("tleap failed to produce prmtop/inpcrd")
    
    print(f"System built: {prmtop}, {inpcrd}")
    return prmtop, inpcrd


def convert_to_openmm(prmtop: Path, inpcrd: Path, out_xml: Path, out_pdb: Path):
    """Convert Amber topology to OpenMM XML for faster loading."""
    from openmm import app
    import openmm
    
    prmtop_obj = app.AmberPrmtopFile(str(prmtop))
    inpcrd_obj = app.AmberInpcrdFile(str(inpcrd))
    
    system = prmtop_obj.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*openmm.unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )
    
    # Save system XML
    with open(out_xml, "w") as f:
        f.write(openmm.XmlSerializer.serialize(system))
    
    # Save first frame as PDB
    with open(out_pdb, "w") as f:
        app.PDBFile.writeFile(prmtop_obj.topology, inpcrd_obj.positions, f)
    
    print(f"OpenMM system: {out_xml}")
    print(f"Reference PDB: {out_pdb}")
    return system


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Input AF3 PDB")
    parser.add_argument("--name", required=True, help="System name")
    parser.add_argument("--out-dir", default="data/md_runs", help="Output directory")
    parser.add_argument("--truncate", action="store_true", help="Enable domain truncation")
    parser.add_argument("--cgas-residues", help="cGAS residue range to keep, e.g. 350-522")
    parser.add_argument("--trim-residues", help="TRIM41 residue range to keep, e.g. 400-630")
    args = parser.parse_args()
    
    work_dir = Path(args.out_dir) / args.name
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Prepare PDB
    raw_pdb = work_dir / "raw_extracted.pdb"
    prepare_pdb(Path(args.pdb), raw_pdb)
    
    # Step 2: pdb4amber
    fixed_pdb = work_dir / "fixed.pdb"
    pdb4amber(raw_pdb, fixed_pdb)
    
    # Step 3: tleap
    prmtop, inpcrd = tleap_build(fixed_pdb, work_dir, args.name)
    
    # Step 4: Convert to OpenMM
    system_xml = work_dir / f"{args.name}_system.xml"
    ref_pdb = work_dir / f"{args.name}_ref.pdb"
    convert_to_openmm(prmtop, inpcrd, system_xml, ref_pdb)
    
    print(f"\n✅ System {args.name} built successfully in {work_dir}")


if __name__ == "__main__":
    main()
