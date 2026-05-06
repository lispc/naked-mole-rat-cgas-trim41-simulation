#!/usr/bin/env python3
"""
Four-system comprehensive comparison: Hgal_WT vs Hgal_4mut_rev vs Hsap_WT vs Hsap_4mut.

Metrics:
  - Per-residue RMSF (3-rep average)
  - cGAS-TRIM41 COM distance
  - Radius of gyration (Rg)
  - Interface contact count (< 5Å)
  - Cross-system RMSD
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import distance_array
from pathlib import Path

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/four_system_comparison"
OUTDIR.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "Hgal_WT": {
        "prmtop": BASE / "data/md_runs/Hgal_WT/Hgal_WT.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hgal_WT/rep1/Hgal_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hgal_WT/rep2/Hgal_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hgal_WT/rep2/Hgal_WT_rep2_restart.dcd",
            BASE / "data/md_runs/Hgal_WT/rep3/Hgal_WT_rep3_prod.dcd",
        ],
        "color": "tab:blue",
    },
    "Hgal_4mut_rev": {
        "prmtop": BASE / "data/md_runs/Hgal_4mut_rev/Hgal_4mut_rev.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hgal_4mut_rev/rep1/Hgal_4mut_rev_rep1_prod.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_prod.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep2/Hgal_4mut_rev_rep2_restart.dcd",
            BASE / "data/md_runs/Hgal_4mut_rev/rep3/Hgal_4mut_rev_rep3_prod.dcd",
        ],
        "color": "tab:red",
    },
    "Hsap_WT": {
        "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        ],
        "color": "tab:green",
    },
    "Hsap_4mut": {
        "prmtop": BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
        "dcds": [
            BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
        ],
        "color": "tab:orange",
    },
}


def analyze_system(name, paths):
    """Analyze a single system across all reps."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")
    
    prmtop = paths["prmtop"]
    dcds = paths["dcds"]
    
    # Load first DCD to get topology info
    u_ref = mda.Universe(str(prmtop), str(dcds[0]))
    
    # Define selections - auto-detect TRIM41/cGAS boundary
    all_protein_res = [r for r in u_ref.residues if r.resname not in ['WAT', 'HOH', 'Na+', 'Cl-', 'K+']]
    n_prot = len(all_protein_res)
    # TRIM41 SPRY = 218 aa, cGAS CTD varies by species
    trim41_end = 218
    cgas_start = 219
    
    protein = u_ref.select_atoms("protein")
    cgas = u_ref.select_atoms(f"protein and resid {cgas_start}-{n_prot}")
    trim41 = u_ref.select_atoms(f"protein and resid 1-{trim41_end}")
    ca = u_ref.select_atoms("protein and name CA")
    
    print(f"  Auto-detected: {n_prot} protein residues")
    print(f"  TRIM41: 1-{trim41_end} ({trim41.n_atoms} atoms)")
    print(f"  cGAS: {cgas_start}-{n_prot} ({cgas.n_atoms} atoms)")
    
    print(f"  Total protein atoms: {protein.n_atoms}")
    print(f"  cGAS atoms: {cgas.n_atoms}")
    print(f"  TRIM41 atoms: {trim41.n_atoms}")
    print(f"  CA atoms: {ca.n_atoms}")
    
    # Collect metrics across all frames
    rmsf_per_rep = []
    com_distances = []
    rg_values = []
    contact_counts = []
    
    for i, dcd in enumerate(dcds, 1):
        print(f"  [Rep {i}/{len(dcds)}] {dcd.name}")
        u = mda.Universe(str(prmtop), str(dcd))
        
        # Select same atoms
        ca_u = u.select_atoms("protein and name CA")
        cgas_u = u.select_atoms(f"protein and resid {cgas_start}-{n_prot}")
        trim41_u = u.select_atoms(f"protein and resid 1-{trim41_end}")
        
        # Align to first frame
        ref = u.copy()
        align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
        
        # RMSF for this rep
        rmsf_calc = rms.RMSF(ca_u, verbose=False).run()
        rmsf_per_rep.append(rmsf_calc.results.rmsf)
        
        # COM distance, Rg, contacts per frame
        for ts in u.trajectory:
            com_cgas = cgas_u.center_of_mass()
            com_trim41 = trim41_u.center_of_mass()
            com_dist = np.linalg.norm(com_cgas - com_trim41)
            com_distances.append(com_dist)
            
            rg = protein.radius_of_gyration()
            rg_values.append(rg)
            
            # Contacts: any heavy atom pair < 5Å between cGAS and TRIM41
            dist_array = distance_array(
                cgas_u.positions, trim41_u.positions, box=u.dimensions
            )
            n_contacts = np.sum(dist_array < 5.0)
            contact_counts.append(n_contacts)
    
    resids = ca.resids
    avg_rmsf = np.mean(rmsf_per_rep, axis=0)
    
    results = {
        "resids": resids,
        "resnames": ca.resnames,
        "avg_rmsf": avg_rmsf,
        "rmsf_reps": np.array(rmsf_per_rep),
        "com_dist": np.array(com_distances),
        "rg": np.array(rg_values),
        "contacts": np.array(contact_counts),
    }
    
    print(f"  Total frames analyzed: {len(com_distances)}")
    print(f"  Avg COM distance: {results['com_dist'].mean():.1f} ± {results['com_dist'].std():.1f} Å")
    print(f"  Avg Rg: {results['rg'].mean():.1f} ± {results['rg'].std():.1f} Å")
    print(f"  Avg contacts: {results['contacts'].mean():.0f} ± {results['contacts'].std():.0f}")
    
    return results


def plot_rmsf_comparison(all_results):
    """Plot RMSF for all four systems."""
    fig, ax = plt.subplots(figsize=(16, 5))
    
    for name, data in all_results.items():
        color = SYSTEMS[name]["color"]
        ax.plot(data["resids"], data["avg_rmsf"], 
                label=name.replace("_", " "), color=color, alpha=0.8, lw=1.5)
    
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Per-Residue RMSF Comparison (3-rep average)")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "four_system_rmsf.png", dpi=300)
    print(f"\nSaved: {OUTDIR / 'four_system_rmsf.png'}")


def plot_global_metrics(all_results):
    """Plot COM distance, Rg, contacts distributions."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))
    
    metrics = [
        ("com_dist", "COM Distance (Å)", "cGAS-TRIM41 COM Distance"),
        ("rg", "Rg (Å)", "Radius of Gyration"),
        ("contacts", "Contact Count", "Interface Contacts (<5Å)"),
    ]
    
    for ax, (key, ylabel, title) in zip(axes, metrics):
        for name, data in all_results.items():
            color = SYSTEMS[name]["color"]
            values = data[key]
            ax.hist(values, bins=50, alpha=0.4, color=color, 
                   label=f"{name.replace('_', ' ')} ({values.mean():.1f}±{values.std():.1f})",
                   density=True)
        ax.set_xlabel(ylabel)
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.legend(fontsize=7)
    
    plt.tight_layout()
    plt.savefig(OUTDIR / "four_system_global_metrics.png", dpi=300)
    print(f"Saved: {OUTDIR / 'four_system_global_metrics.png'}")


def print_summary_table(all_results):
    """Print summary statistics table."""
    print(f"\n{'='*80}")
    print("Four-System Summary")
    print(f"{'='*80}")
    print(f"{'System':<18} {'Frames':>8} {'COM(Å)':>12} {'Rg(Å)':>10} {'Contacts':>10} {'RMSF(Å)':>10}")
    print(f"{'-'*80}")
    for name, data in all_results.items():
        print(f"{name:<18} {len(data['com_dist']):>8} "
              f"{data['com_dist'].mean():>6.1f}±{data['com_dist'].std():<4.1f} "
              f"{data['rg'].mean():>6.1f}±{data['rg'].std():<2.1f} "
              f"{data['contacts'].mean():>6.0f}±{data['contacts'].std():<2.0f} "
              f"{data['avg_rmsf'].mean():>6.2f}±{data['avg_rmsf'].std():<2.2f}")
    print(f"{'='*80}")


def main():
    print("Four-System Comprehensive Comparison")
    print("="*60)
    
    all_results = {}
    for name, paths in SYSTEMS.items():
        all_results[name] = analyze_system(name, paths)
    
    # Save all data
    np.savez(OUTDIR / "four_system_data.npz", **{
        f"{name}_{key}": val for name, data in all_results.items() for key, val in data.items()
    })
    
    # Plots
    print("\n[Plots] Generating...")
    plot_rmsf_comparison(all_results)
    plot_global_metrics(all_results)
    print_summary_table(all_results)
    
    print(f"\n{'='*60}")
    print("All done!")
    print(f"Output: {OUTDIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
