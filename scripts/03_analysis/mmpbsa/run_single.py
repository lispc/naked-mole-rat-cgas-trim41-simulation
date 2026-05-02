#!/usr/bin/env python3
"""
Run MM-GBSA on MD trajectory using MMPBSA.py (AmberTools).

Usage:
  python scripts/run_mmpbsa.py \
      --topology data/md_runs/Hsap_WT/Hsap_WT.prmtop \
      --trajectory data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd \
      --receptor-mask ":1-218" \
      --ligand-mask ":219-1152" \
      --name Hsap_WT_rep1 \
      --method gb \
      --interval 20 \
      --outdir data/analysis/mmpbsa

Notes:
  - DCD is converted to NetCDF via cpptraj for MMPBSA.py compatibility.
  - --interval N: sample every Nth frame from trajectory.
  - For multi-replica, run separately then average results.
"""
import argparse
import os
import subprocess
import sys
import json
from pathlib import Path


def find_amber_bin(name):
    """Find Amber binary in conda env or PATH."""
    conda_bin = Path.home() / "miniforge3/envs/cgas-md/bin" / name
    if conda_bin.exists():
        return str(conda_bin)
    return shutil.which(name) or name


def write_mmpbsa_input(out_path, method="gb", decomp=True, interval=20,
                         receptor_mask=None, ligand_mask=None):
    """Write MMPBSA.py input file."""
    if method == "gb":
        gb_section = """
&gb
  igb=5, saltcon=0.150,
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
  print_res="all",
/
""" if decomp else ""

    # Build general section
    mask_line = ""
    if receptor_mask:
        mask_line += f"\n  receptor_mask=\"{receptor_mask}\","
    if ligand_mask:
        mask_line += f"\n  ligand_mask=\"{ligand_mask}\","

    content = f"""&general
  startframe=1, endframe=99999, interval={interval},
  keep_files=0, strip_mask=":WAT,Cl-,Na+",{mask_line}
  use_sander=0,
/
{gb_section}{decomp_section}"""

    with open(out_path, "w") as f:
        f.write(content)


def convert_dcd_to_nc(topology, trajectory, out_nc):
    """Use cpptraj to convert DCD to NetCDF."""
    cpptraj = find_amber_bin("cpptraj")
    traj_stem = Path(trajectory).stem
    
    cpptraj_in = out_nc.parent / f"{traj_stem}_convert.in"
    with open(cpptraj_in, "w") as f:
        f.write(f"parm {topology}\n")
        f.write(f"trajin {trajectory}\n")
        f.write(f"trajout {out_nc} netcdf\n")
        f.write("run\n")

    print(f"Converting {Path(trajectory).name} → NetCDF...")
    result = subprocess.run(
        [cpptraj, "-i", str(cpptraj_in)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"cpptraj stderr:\n{result.stderr[-500:]}")
        raise RuntimeError("cpptraj conversion failed")
    
    if not out_nc.exists():
        raise RuntimeError(f"Output NetCDF not created: {out_nc}")
    
    print(f"  ✅ Converted: {out_nc.name}")
    return out_nc


def parse_mmpbsa_results(dat_file):
    """Parse MMPBSA.py results .dat file for summary energies."""
    results = {}
    with open(dat_file) as f:
        lines = f.readlines()
    
    # Find DELTA G section
    in_delta = False
    for line in lines:
        if "DELTA G" in line or "DELTA TOTAL" in line:
            in_delta = True
        if in_delta and line.strip() and not line.startswith("-"):
            parts = line.split()
            if len(parts) >= 3:
                key = parts[0]
                try:
                    val = float(parts[1])
                    results[key] = val
                except ValueError:
                    pass
    return results


def run_mmpbsa(topology, trajectory, receptor_mask, ligand_mask, 
               name, out_dir, method="gb", interval=20, decomp=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Convert DCD to NetCDF
    nc_path = out_dir / f"{Path(trajectory).stem}.nc"
    if not nc_path.exists():
        convert_dcd_to_nc(topology, trajectory, nc_path)
    else:
        print(f"Using existing NetCDF: {nc_path.name}")

    # Write MMPBSA input (with masks)
    mmpbsa_in = out_dir / f"{name}_mmpbsa.in"
    write_mmpbsa_input(mmpbsa_in, method=method, decomp=decomp, interval=interval,
                       receptor_mask=receptor_mask, ligand_mask=ligand_mask)

    # Run MMPBSA.py
    mmpbsa_py = find_amber_bin("MMPBSA.py")
    cmd = [
        sys.executable, mmpbsa_py,
        "-i", str(mmpbsa_in),
        "-o", str(out_dir / f"{name}_results.dat"),
        "-do", str(out_dir / f"{name}_decomp.dat"),
        "-cp", topology,
        "-y", str(nc_path),
    ]

    print(f"\nRunning MMPBSA.py ({method.upper()}, interval={interval})...")
    print(f"$ {' '.join(cmd[-8:])}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    log_path = out_dir / f"{name}_mmpbsa.log"
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    if result.returncode != 0:
        print(f"MMPBSA.py failed. Log: {log_path}")
        print(f"STDERR tail:\n{result.stderr[-1000:]}")
        raise RuntimeError("MMPBSA.py failed")

    results_dat = out_dir / f"{name}_results.dat"
    decomp_dat = out_dir / f"{name}_decomp.dat"
    
    print(f"\n✅ Results: {results_dat.name}")
    print(f"   Decomp:  {decomp_dat.name}")

    # Parse summary
    if results_dat.exists():
        summary = parse_mmpbsa_results(results_dat)
        print(f"\n   ΔG binding components:")
        for k, v in summary.items():
            print(f"     {k}: {v:.2f} kcal/mol")
        
        # Save JSON
        json_path = out_dir / f"{name}_summary.json"
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"   JSON:    {json_path.name}")

    return results_dat, decomp_dat


def main():
    parser = argparse.ArgumentParser(description="Run MM-GBSA on MD trajectory")
    parser.add_argument("--topology", required=True, help="Amber prmtop file")
    parser.add_argument("--trajectory", required=True, help="MD trajectory (DCD)")
    parser.add_argument("--receptor-mask", required=True, 
                        help="Amber mask for receptor (e.g., ':1-218')")
    parser.add_argument("--ligand-mask", required=True,
                        help="Amber mask for ligand (e.g., ':219-1152')")
    parser.add_argument("--name", required=True, help="Run identifier")
    parser.add_argument("--out-dir", default="data/analysis/mmpbsa")
    parser.add_argument("--method", choices=["gb", "pb"], default="gb")
    parser.add_argument("--interval", type=int, default=20,
                        help="Frame sampling interval (default: every 20th frame)")
    parser.add_argument("--no-decomp", action="store_true",
                        help="Skip per-residue decomposition")
    args = parser.parse_args()

    run_mmpbsa(
        topology=args.topology,
        trajectory=args.trajectory,
        receptor_mask=args.receptor_mask,
        ligand_mask=args.ligand_mask,
        name=args.name,
        out_dir=args.out_dir,
        method=args.method,
        interval=args.interval,
        decomp=not args.no_decomp,
    )


if __name__ == "__main__":
    import shutil
    main()
