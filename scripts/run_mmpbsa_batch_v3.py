#!/usr/bin/env python3
"""
Batch MM-GBSA runner v3 — uses protein-only prmtop + protein-only trajectory.

This avoids the "bad atom type: EP" error from OPC water extra points
by running MMPBSA.py on stripped protein systems.

Usage:
    python run_mmpbsa_batch_v3.py

Requires AMBERHOME env var (or auto-detects from sander path).
"""

import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
MD_DIR   = BASE_DIR / "data" / "md_runs"
OUT_DIR  = BASE_DIR / "data" / "analysis" / "mmpbsa"

JOBS = [
    # (name, protein_prmtop, protein_nc)
    ("Hsap_WT_rep2",
     MD_DIR / "Hsap_WT" / "Hsap_WT_protein.prmtop",
     OUT_DIR / "Hsap_WT_rep2_prot.nc"),
    ("Hsap_WT_rep3",
     MD_DIR / "Hsap_WT" / "Hsap_WT_protein.prmtop",
     OUT_DIR / "Hsap_WT_rep3_prot.nc"),
    ("Hsap_WT_S305phos_rep1",
     MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos_protein.prmtop",
     OUT_DIR / "Hsap_WT_S305phos_rep1_prot.nc"),
    ("Hsap_WT_S305phos_rep2",
     MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos_protein.prmtop",
     OUT_DIR / "Hsap_WT_S305phos_rep2_prot.nc"),
    ("Hsap_WT_S305phos_rep3",
     MD_DIR / "Hsap_WT_S305phos" / "Hsap_WT_S305phos_protein.prmtop",
     OUT_DIR / "Hsap_WT_S305phos_rep3_prot.nc"),
]

RECEPTOR_MASK = ":1-218"
LIGAND_MASK   = ":219-541"
INTERVAL      = 50   # analyse every 50th frame
# ---------------------------------------------------------------------------


def find_amber_bin(name: str) -> Path:
    """Find an AMBER executable."""
    for p in [Path(os.environ.get("AMBERHOME", "")), Path(sys.executable).parent]:
        candidate = p / "bin" / name
        if candidate.exists():
            return candidate
    candidate = Path(sys.executable).parent / name
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"Cannot find {name}")


def run_single(job_name: str, prmtop: Path, nc: Path, out_dir: Path) -> dict:
    """Run MMPBSA.py for a single replica."""
    mmpbsa_py = find_amber_bin("MMPBSA.py")

    # Write MMPBSA input
    mmpbsa_in = out_dir / f"{job_name}_mmpbsa.in"
    with open(mmpbsa_in, "w") as f:
        f.write(f"""&general
startframe=1, endframe=99999, interval={INTERVAL},
keep_files=0,
receptor_mask="{RECEPTOR_MASK}",
ligand_mask="{LIGAND_MASK}",
/
&gb
igb=5, saltcon=0.150,
/
&decomp
idecomp=2, dec_verbose=1,
print_res="all",
/
""")

    results_dat = out_dir / f"{job_name}_results.dat"
    decomp_dat  = out_dir / f"{job_name}_decomp.dat"
    log_file    = out_dir / f"{job_name}_mmpbsa.log"

    cmd = [
        sys.executable, str(mmpbsa_py),
        "-i", str(mmpbsa_in),
        "-o", str(results_dat),
        "-do", str(decomp_dat),
        "-cp", str(prmtop),
        "-y", str(nc),
    ]

    env = {**dict(os.environ), "AMBERHOME": str(mmpbsa_py.parent.parent)}
    # Limit threads per job to avoid oversubscription when parallel
    env["OMP_NUM_THREADS"] = "4"

    print(f"[{job_name}] Starting MMPBSA.py ...")
    print(f"[{job_name}] {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(out_dir),
        env=env,
    )

    with open(log_file, "w") as f:
        f.write(result.stdout)
        if result.stderr:
            f.write("\n--- STDERR ---\n")
            f.write(result.stderr)

    if result.returncode == 0:
        print(f"[{job_name}] ✅ MMPBSA.py finished successfully")
        return {"name": job_name, "ok": True, "log": log_file}
    else:
        print(f"[{job_name}] ❌ MMPBSA.py failed (exit {result.returncode})")
        print(f"[{job_name}] Log: {log_file}")
        # Print tail of stderr for quick diagnosis
        err_tail = result.stderr[-800:] if len(result.stderr) > 800 else result.stderr
        if err_tail:
            print(f"[{job_name}] STDERR tail:\n{err_tail}")
        return {"name": job_name, "ok": False, "log": log_file, "stderr": result.stderr}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Ensure AMBERHOME is discoverable
    try:
        sander = find_amber_bin("sander")
        os.environ.setdefault("AMBERHOME", str(sander.parent.parent))
    except FileNotFoundError:
        pass

    print(f"Output directory: {OUT_DIR}")
    print(f"Jobs to run: {len(JOBS)}")
    print(f"Interval: every {INTERVAL} frames")
    print(f"Receptor mask: {RECEPTOR_MASK}")
    print(f"Ligand mask: {LIGAND_MASK}")
    print("-" * 60)

    # Run serially — MMPBSA.py/sander already parallelises internally.
    # Running multiple jobs in parallel causes oversubscription on a 64-core
    # machine and can actually slow things down.
    results = []
    for name, prmtop, nc in JOBS:
        if not prmtop.exists():
            print(f"[{name}] Skipping — prmtop not found: {prmtop}")
            results.append({"name": name, "ok": False, "reason": "prmtop missing"})
            continue
        if not nc.exists():
            print(f"[{name}] Skipping — trajectory not found: {nc}")
            results.append({"name": name, "ok": False, "reason": "nc missing"})
            continue
        res = run_single(name, prmtop, nc, OUT_DIR)
        results.append(res)
        print("-" * 60)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    ok_count = sum(1 for r in results if r.get("ok"))
    fail_count = len(results) - ok_count
    for r in results:
        status = "✅ PASS" if r.get("ok") else "❌ FAIL"
        print(f"  {status}  {r['name']}")
    print(f"\nTotal: {ok_count}/{len(results)} succeeded")

    if fail_count:
        sys.exit(1)


if __name__ == "__main__":
    main()
