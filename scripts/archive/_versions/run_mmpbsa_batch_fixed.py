#!/usr/bin/env python3
"""
Fixed batch MM-GBSA runner.

Uses original prmtop + strip_mask + receptor/ligand masks (same as run_mmpbsa.py).
Runs sequentially to avoid temp file conflicts.

Usage:
    python scripts/run_mmpbsa_batch_fixed.py
"""
import os
import subprocess
import sys
from pathlib import Path

AMBERHOME = Path.home() / "miniforge3/envs/cgas-md"
OUTDIR = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/analysis/mmpbsa")


def run_one(name, topology, trajectory, nc_path, receptor_mask, ligand_mask):
    """Run MMPBSA.py for one replica using original prmtop + masks."""
    outdir = OUTDIR / name
    outdir.mkdir(parents=True, exist_ok=True)
    workdir = outdir / "_work"
    workdir.mkdir(parents=True, exist_ok=True)

    # Write MMPBSA input with masks
    inp = workdir / "mmpbsa.in"
    with open(inp, 'w') as f:
        f.write(f"""&general
  startframe=1, endframe=99999, interval=50,
  keep_files=0, strip_mask=":WAT,Cl-,Na+",
  receptor_mask="{receptor_mask}",
  ligand_mask="{ligand_mask}",
  use_sander=0,
/
&gb
  igb=5, saltcon=0.150,
/
""")

    mmpbsa_py = str(AMBERHOME / "bin/MMPBSA.py")
    cmd = [
        sys.executable, mmpbsa_py,
        "-i", str(inp),
        "-o", str(outdir / "results.dat"),
        "-cp", str(topology),
        "-y", str(nc_path),
    ]

    log = outdir / "mmpbsa.log"
    print(f"\n{'='*60}")
    print(f"[START] {name}")
    print(f"{'='*60}")

    env = {**os.environ, 'AMBERHOME': str(AMBERHOME)}
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd=str(workdir), env=env)

    with open(log, 'w') as f:
        f.write(result.stdout)
        f.write(result.stderr)

    results_file = outdir / "results.dat"
    if results_file.exists():
        content = results_file.read_text()
        if "DELTA TOTAL" in content:
            # Extract value
            for line in content.splitlines():
                if "DELTA TOTAL" in line:
                    print(f"[DONE] {name}: {line.strip()}")
                    return (name, 0)

    print(f"[FAIL] {name}")
    print(f"  Log: {log}")
    print(f"  STDERR tail: {result.stderr[-500:]}")
    return (name, 1)


def main():
    jobs = []

    # Hsap_WT rep2, rep3
    for rep in [2, 3]:
        name = f"Hsap_WT_rep{rep}"
        topology = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT/Hsap_WT.prmtop"
        nc = OUTDIR / f"Hsap_WT_rep{rep}_prot.nc"
        jobs.append((name, topology, None, nc, ":1-218", ":219-541"))

    # S305-phos rep1, rep2, rep3
    for rep in [1, 2, 3]:
        name = f"Hsap_WT_S305phos_rep{rep}"
        topology = "/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation/data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop"
        nc = OUTDIR / f"Hsap_WT_S305phos_rep{rep}_prot.nc"
        jobs.append((name, topology, None, nc, ":1-218", ":219-541"))

    print(f"Running {len(jobs)} MM-GBSA jobs sequentially...")
    print(f"Estimated time: ~{len(jobs) * 30} min ({len(jobs) * 0.5}h)")

    success = 0
    for job in jobs:
        name, code = run_one(*job)
        success += (code == 0)

    print(f"\n{'='*60}")
    print(f"Batch complete: {success}/{len(jobs)} succeeded")
    print(f"{'='*60}")

    for name, _, _, _, _, _ in jobs:
        results_file = OUTDIR / name / "results.dat"
        if results_file.exists() and "DELTA TOTAL" in results_file.read_text():
            for line in results_file.read_text().splitlines():
                if "DELTA TOTAL" in line:
                    print(f"  ✅ {name}: {line.strip()}")
                    break
        else:
            print(f"  ❌ {name}: failed")


if __name__ == '__main__':
    main()
