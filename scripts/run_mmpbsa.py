#!/usr/bin/env python3
"""
Run MM-GBSA / MM-PBSA on MD trajectory using MMPBSA.py (AmberTools).

Usage:
  python scripts/run_mmpbsa.py \
      --topology data/md_runs/Hsap_WT/Hsap_WT.prmtop \
      --trajectory data/md_runs/Hsap_WT/production/Hsap_WT_rep1.dcd \
      --receptor-mask ":1-522" \
      --ligand-mask ":523-1152" \
      --name Hsap_WT_rep1 \
      --method gb
"""
import argparse
import subprocess
import tempfile
from pathlib import Path


def write_mmpbsa_input(out_path, method="gb", decomp=True):
    """Write MMPBSA.py input file."""
    if method == "gb":
        gb_section = """
&gb
  igb=5, saltcon=0.150,
  molsurf=0, sprob=1.4,
/
"""
    elif method == "pb":
        gb_section = """
&pb
  istrng=0.150, fillratio=4.0,
/
"""
    else:
        raise ValueError(f"Unknown method: {method}")
    
    decomp_section = """
&decomp
  idecomp=2, dec_verbose=1,
  per_residue=1,
/
""" if decomp else ""

    content = f"""&general
  startframe=1, endframe=0, interval=20,
  keep_files=0, strip_mask=":WAT,Cl-,Na+",
  use_sander=0,
/
{gb_section}{decomp_section}"""
    
    with open(out_path, "w") as f:
        f.write(content)


def run_mmpbsa(topology, trajectory, receptor_mask, ligand_mask, name, out_dir, method="gb"):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert DCD to NetCDF if needed (MMPBSA.py prefers nc or mdcrd)
    traj_stem = Path(trajectory).stem
    nc_path = out_dir / f"{traj_stem}.nc"
    
    # Use cpptraj to convert
    cpptraj_in = out_dir / "convert.in"
    with open(cpptraj_in, "w") as f:
        f.write(f"parm {topology}\n")
        f.write(f"trajin {trajectory}\n")
        f.write(f"trajout {nc_path} netcdf\n")
        f.write("run\n")
    
    print("Converting trajectory to NetCDF...")
    subprocess.run(["cpptraj", "-i", str(cpptraj_in)], check=True, capture_output=True)
    
    # Write MMPBSA input
    mmpbsa_in = out_dir / f"{name}_mmpbsa.in"
    write_mmpbsa_input(mmpbsa_in, method=method)
    
    # Run MMPBSA.py
    cmd = [
        "MMPBSA.py",
        "-i", str(mmpbsa_in),
        "-o", str(out_dir / f"{name}_results.dat"),
        "-do", str(out_dir / f"{name}_decomp.dat"),
        "-cp", topology,
        "-y", str(nc_path),
        "-yr", receptor_mask,
        "-yl", ligand_mask,
    ]
    
    print(f"Running MMPBSA.py ({method})...")
    print(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"STDERR:\n{result.stderr}")
        raise RuntimeError("MMPBSA.py failed")
    
    print(f"Results: {out_dir / f'{name}_results.dat'}")
    print(f"Decomposition: {out_dir / f'{name}_decomp.dat'}")
    return out_dir / f"{name}_results.dat"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--receptor-mask", required=True, help="Amber mask for receptor (cGAS)")
    parser.add_argument("--ligand-mask", required=True, help="Amber mask for ligand (TRIM41)")
    parser.add_argument("--name", required=True)
    parser.add_argument("--out-dir", default="data/analysis/mmpbsa")
    parser.add_argument("--method", choices=["gb", "pb"], default="gb")
    args = parser.parse_args()
    
    run_mmpbsa(
        topology=args.topology,
        trajectory=args.trajectory,
        receptor_mask=args.receptor_mask,
        ligand_mask=args.ligand_mask,
        name=args.name,
        out_dir=args.out_dir,
        method=args.method,
    )


if __name__ == "__main__":
    main()
