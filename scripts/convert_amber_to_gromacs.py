#!/usr/bin/env python3
"""
Convert AMBER prmtop + NPT trajectory last frame to GROMACS top + gro.
Uses MDAnalysis for coordinate extraction (preserves EPW + box),
parmed for topology conversion.
"""
import numpy as np
from pathlib import Path
import MDAnalysis as mda
import parmed as pmd

SYSTEMS = {
    "Hsap_WT": {
        "prmtop": "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        "npt_dcd": "data/md_runs/Hsap_WT/rep{rep}/Hsap_WT_rep{rep}_npt.dcd",
        "out_dir": "data/md_runs_gmx/Hsap_WT/rep{rep}",
    },
    "Hsap_4mut": {
        "prmtop": "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
        "npt_dcd": "data/md_runs/Hsap_4mut/rep{rep}/Hsap_4mut_rep{rep}_npt.dcd",
        "out_dir": "data/md_runs_gmx/Hsap_4mut/rep{rep}",
    },
}

def convert_topology(prmtop_path, dcd_path, out_top, out_gro):
    """Convert AMBER topology to GROMACS using parmed + MDAnalysis coordinates."""
    
    # Load topology with parmed
    amber = pmd.load_file(prmtop_path)
    
    # Load trajectory with MDAnalysis and get last frame
    u = mda.Universe(prmtop_path, dcd_path)
    last_ts = u.trajectory[-1]
    
    # Set coordinates (MDAnalysis uses Angstrom, parmed uses Angstrom)
    for i, atom in enumerate(amber.atoms):
        atom.xx = u.atoms[i].position[0]
        atom.xy = u.atoms[i].position[1]
        atom.xz = u.atoms[i].position[2]
    
    # Set box: MDAnalysis box is [a, b, c, alpha, beta, gamma] in Angstrom, degrees
    # parmed expects [a, b, c, alpha, beta, gamma] in Angstrom, degrees
    box = last_ts.dimensions
    amber.box = box
    
    print(f"  Atoms: {len(amber.atoms)}, Box: {box}")
    
    # Save GROMACS topology and coordinates
    amber.save(out_top, overwrite=True)
    amber.save(out_gro, overwrite=True)
    return True

def main():
    for system_name, paths in SYSTEMS.items():
        for rep in [1, 2]:
            print(f"\n=== {system_name} rep{rep} ===")
            
            prmtop = paths["prmtop"]
            npt_dcd = paths["npt_dcd"].format(rep=rep)
            out_dir = Path(paths["out_dir"].format(rep=rep))
            out_dir.mkdir(parents=True, exist_ok=True)
            
            out_top = out_dir / f"{system_name}_rep{rep}.top"
            out_gro = out_dir / f"{system_name}_rep{rep}.gro"
            
            if not Path(prmtop).exists():
                print(f"Missing prmtop: {prmtop}")
                continue
            if not Path(npt_dcd).exists():
                print(f"Missing NPT dcd: {npt_dcd}")
                continue
            
            print("Converting to GROMACS...")
            ok = convert_topology(prmtop, npt_dcd, str(out_top), str(out_gro))
            if ok:
                print(f"Saved: {out_top}, {out_gro}")
            else:
                print("Conversion failed")

if __name__ == "__main__":
    main()
