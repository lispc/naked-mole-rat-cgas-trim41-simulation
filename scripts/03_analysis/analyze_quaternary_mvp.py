#!/usr/bin/env python3
"""
Analyze quaternary MVP trajectories (E2~Ub-TRIM41-cGAS).

Metrics:
- K315 NZ -> Ub G76 C distance (primary catalytic geometry)
- K315 NZ -> E2 K85 NZ distance (alternative)
- E2 K85 NZ -> Ub G76 C distance (E2~Ub closed conformation)
- K315 side-chain SASA
- RING-cGAS interface contacts
- RMSD of RING dimer

Expected values (from experimental design):
                WT      4mut
K315->Ub G76    ~12 Å   ~7 Å
E2~Ub closed    ~30%    ~60%
K315 SASA       Low     High
"""
import sys
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.distances import distance_array

BASE = Path("/home/scroll/personal/naked-mole-rat-cgas-trim41-simulation")
OUTDIR = BASE / "data/analysis/quaternary_mvp"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Residue numbers in prmtop (from clean_renum.txt)
RES_K315 = 962       # cGAS K315 (chain C)
RES_UB_G76 = 421     # Ub G76 (chain U)
RES_E2_K85 = 273     # E2 K85 (chain E)

# Atom selections
SEL_K315_NZ = f"resid {RES_K315} and name NZ"
SEL_UB_G76_C = f"resid {RES_UB_G76} and name C"
SEL_E2_K85_NZ = f"resid {RES_E2_K85} and name NZ"
SEL_K315_SC = f"resid {RES_K315} and (name CB CG CD CE NZ)"
SEL_RING = "chainID R S and name CA"
SEL_CGAS = "chainID C and name CA"
SEL_SPRY = "chainID P and name CA"


def compute_sasa_atomgroup(ag, probe_radius=1.4):
    """Compute SASA for an atomgroup using Shrake-Rupley algorithm."""
    try:
        from MDAnalysis.analysis.hole2 import sr
        # Not available in all MDA versions
        return None
    except ImportError:
        pass

    # Simple Shrake-Rupley implementation
    coords = ag.positions
    radii = np.array([1.7 if a.name.startswith('C') else 1.5 if a.name.startswith('N') else 1.2 for a in ag])
    
    # Golden spiral sphere sampling
    n_points = 960
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    
    sphere = np.zeros((n_points, 3))
    sphere[:, 0] = np.sin(phi) * np.cos(theta)
    sphere[:, 1] = np.sin(phi) * np.sin(theta)
    sphere[:, 2] = np.cos(phi)
    
    total_sasa = 0.0
    for i, (c, r) in enumerate(zip(coords, radii)):
        probe_r = r + probe_radius
        points = c + probe_r * sphere
        # Check which points are exposed
        other_coords = np.delete(coords, i, axis=0)
        other_radii = np.delete(radii, i, axis=0)
        
        exposed = 0
        for p in points:
            dists = np.linalg.norm(other_coords - p, axis=1)
            if np.all(dists > other_radii + probe_radius):
                exposed += 1
        
        total_sasa += 4 * np.pi * probe_r**2 * (exposed / n_points)
    
    return total_sasa


def analyze_trajectory(prmtop, dcd, name):
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")

    u = mda.Universe(str(prmtop), str(dcd))
    n_frames = len(u.trajectory)
    print(f"Frames: {n_frames}")

    # Selections
    k315_nz = u.select_atoms(SEL_K315_NZ)
    ub_g76_c = u.select_atoms(SEL_UB_G76_C)
    e2_k85_nz = u.select_atoms(SEL_E2_K85_NZ)
    k315_sc = u.select_atoms(SEL_K315_SC)
    ring_ca = u.select_atoms(SEL_RING)
    cgas_ca = u.select_atoms(SEL_CGAS)
    spry_ca = u.select_atoms(SEL_SPRY)

    # Reference for RMSD
    ref = mda.Universe(str(prmtop), str(dcd))
    ref_ring_ca = ref.select_atoms(SEL_RING)

    # Arrays
    time_ns = np.zeros(n_frames)
    d_k315_ub = np.zeros(n_frames)
    d_k315_e2 = np.zeros(n_frames)
    d_e2_ub = np.zeros(n_frames)
    k315_sasa = np.zeros(n_frames)
    ring_rmsd = np.zeros(n_frames)
    ring_cgas_contacts = np.zeros(n_frames)
    spry_cgas_contacts = np.zeros(n_frames)

    for i, ts in enumerate(u.trajectory):
        time_ns[i] = ts.time / 1000.0  # ps -> ns

        # Distances
        d_k315_ub[i] = np.linalg.norm(k315_nz.center_of_geometry() - ub_g76_c.center_of_geometry())
        d_k315_e2[i] = np.linalg.norm(k315_nz.center_of_geometry() - e2_k85_nz.center_of_geometry())
        d_e2_ub[i] = np.linalg.norm(e2_k85_nz.center_of_geometry() - ub_g76_c.center_of_geometry())

        # RMSD
        ring_rmsd[i] = rms.rmsd(ring_ca.positions, ref_ring_ca.positions, superposition=False)

        # Contacts (heavy atoms within 5 Å)
        ring_heavy = u.select_atoms("chainID R S and not name H*")
        cgas_heavy = u.select_atoms("chainID C and not name H*")
        spry_heavy = u.select_atoms("chainID P and not name H*")

        dist_ring_cgas = distance_array(ring_heavy.positions, cgas_heavy.positions)
        ring_cgas_contacts[i] = np.sum(dist_ring_cgas < 5.0)

        dist_spry_cgas = distance_array(spry_heavy.positions, cgas_heavy.positions)
        spry_cgas_contacts[i] = np.sum(dist_spry_cgas < 5.0)

        if i % 100 == 0:
            print(f"  Frame {i}/{n_frames}: K315->Ub={d_k315_ub[i]:.1f}Å, E2~Ub={d_e2_ub[i]:.1f}Å")

    results = {
        'name': name,
        'time_ns': time_ns,
        'd_k315_ub': d_k315_ub,
        'd_k315_e2': d_k315_e2,
        'd_e2_ub': d_e2_ub,
        'ring_rmsd': ring_rmsd,
        'ring_cgas_contacts': ring_cgas_contacts,
        'spry_cgas_contacts': spry_cgas_contacts,
    }

    # Summary stats
    print(f"\nSummary ({name}):")
    print(f"  K315->Ub G76:  {np.mean(d_k315_ub):.2f} ± {np.std(d_k315_ub):.2f} Å")
    print(f"  K315->E2 K85:  {np.mean(d_k315_e2):.2f} ± {np.std(d_k315_e2):.2f} Å")
    print(f"  E2 K85->Ub G76:{np.mean(d_e2_ub):.2f} ± {np.std(d_e2_ub):.2f} Å")
    print(f"  RING RMSD:     {np.mean(ring_rmsd):.2f} ± {np.std(ring_rmsd):.2f} Å")
    print(f"  RING-cGAS contacts: {np.mean(ring_cgas_contacts):.0f} ± {np.std(ring_cgas_contacts):.0f}")
    print(f"  SPRY-cGAS contacts: {np.mean(spry_cgas_contacts):.0f} ± {np.std(spry_cgas_contacts):.0f}")

    # Closed fraction: E2~Ub distance < 8 Å (generous threshold for closed state)
    closed_frac = np.mean(d_e2_ub < 8.0) * 100
    print(f"  E2~Ub closed fraction (<8Å): {closed_frac:.1f}%")
    results['closed_fraction'] = closed_frac

    return results


def plot_comparison(wt, mut):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # K315->Ub distance
    ax = axes[0, 0]
    ax.plot(wt['time_ns'], wt['d_k315_ub'], label='WT', alpha=0.7)
    ax.plot(mut['time_ns'], mut['d_k315_ub'], label='4mut', alpha=0.7)
    ax.axhline(7, color='green', linestyle='--', alpha=0.5, label='Target (~7Å)')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('K315 NZ → Ub G76 C (Å)')
    ax.set_title('Primary Catalytic Distance')
    ax.legend()

    # E2~Ub distance
    ax = axes[0, 1]
    ax.plot(wt['time_ns'], wt['d_e2_ub'], label='WT', alpha=0.7)
    ax.plot(mut['time_ns'], mut['d_e2_ub'], label='4mut', alpha=0.7)
    ax.axhline(8, color='red', linestyle='--', alpha=0.5, label='Closed threshold')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('E2 K85 NZ → Ub G76 C (Å)')
    ax.set_title('E2~Ub Conformation')
    ax.legend()

    # Interface contacts
    ax = axes[1, 0]
    ax.plot(wt['time_ns'], wt['ring_cgas_contacts'], label='WT RING-cGAS', alpha=0.7)
    ax.plot(mut['time_ns'], mut['ring_cgas_contacts'], label='4mut RING-cGAS', alpha=0.7)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Heavy-atom contacts (<5Å)')
    ax.set_title('RING-cGAS Interface Contacts')
    ax.legend()

    # RMSD
    ax = axes[1, 1]
    ax.plot(wt['time_ns'], wt['ring_rmsd'], label='WT', alpha=0.7)
    ax.plot(mut['time_ns'], mut['ring_rmsd'], label='4mut', alpha=0.7)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RING CA RMSD (Å)')
    ax.set_title('RING Domain Stability')
    ax.legend()

    plt.tight_layout()
    plt.savefig(OUTDIR / 'quaternary_mvp_comparison.png', dpi=150)
    plt.close()
    print(f"Saved comparison plot to {OUTDIR / 'quaternary_mvp_comparison.png'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wt-prmtop', default=str(BASE / 'data/structures/quaternary_mvp/quaternary_mvp.prmtop'))
    parser.add_argument('--wt-dcd', default=str(BASE / 'data/md_runs/quaternary_mvp/WT_rep1/quaternary_mvp_WT_rep1.dcd'))
    parser.add_argument('--mut-prmtop', default=str(BASE / 'data/structures/quaternary_mvp/quaternary_mvp_4mut.prmtop'))
    parser.add_argument('--mut-dcd', default=str(BASE / 'data/md_runs/quaternary_mvp/4mut_rep1/quaternary_mvp_4mut_rep1.dcd'))
    args = parser.parse_args()

    wt = analyze_trajectory(args.wt_prmtop, args.wt_dcd, 'WT')
    mut = analyze_trajectory(args.mut_prmtop, args.mut_dcd, '4mut')

    plot_comparison(wt, mut)

    # Save summary table
    with open(OUTDIR / 'quaternary_mvp_summary.txt', 'w') as f:
        f.write("Quaternary MVP Analysis Summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"{'Metric':<30} {'WT':>12} {'4mut':>12}\n")
        f.write("-"*60 + "\n")
        f.write(f"{'K315->Ub G76 (Å)':<30} {np.mean(wt['d_k315_ub']):>12.2f} {np.mean(mut['d_k315_ub']):>12.2f}\n")
        f.write(f"{'K315->E2 K85 (Å)':<30} {np.mean(wt['d_k315_e2']):>12.2f} {np.mean(mut['d_k315_e2']):>12.2f}\n")
        f.write(f"{'E2~Ub closed fraction (%)':<30} {wt['closed_fraction']:>12.1f} {mut['closed_fraction']:>12.1f}\n")
        f.write(f"{'RING RMSD (Å)':<30} {np.mean(wt['ring_rmsd']):>12.2f} {np.mean(mut['ring_rmsd']):>12.2f}\n")
        f.write(f"{'RING-cGAS contacts':<30} {np.mean(wt['ring_cgas_contacts']):>12.0f} {np.mean(mut['ring_cgas_contacts']):>12.0f}\n")
        f.write(f"{'SPRY-cGAS contacts':<30} {np.mean(wt['spry_cgas_contacts']):>12.0f} {np.mean(mut['spry_cgas_contacts']):>12.0f}\n")

    print(f"\nSaved summary to {OUTDIR / 'quaternary_mvp_summary.txt'}")


if __name__ == '__main__':
    main()
