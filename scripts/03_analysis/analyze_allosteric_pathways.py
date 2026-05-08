#!/usr/bin/env python3
"""
P1-5: Allosteric pathway identification from 4mut sites to N-terminal interface.

Efficient approach:
  - Uses existing DCCM data (pre-computed from binary MD)
  - Computes CA-level contact maps (CA-CA < 8 Å for >50% frames)
  - Builds weighted graph: edge weight = -log(|correlation|)
  - Finds shortest paths from 4mut sites to N-term interface
  - Compares WT vs 4mut pathway architectures
"""
import json
import warnings
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import distance_array

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

BASE = Path(__file__).resolve().parent.parent.parent
OUTDIR = BASE / "data/analysis/allosteric_pathways"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Config ──────────────────────────────────────────────────────────────
SYSTEMS = {
    "WT": {
        "prmtop": BASE / "data/md_runs/Hsap_WT/Hsap_WT.prmtop",
        "trajs": [
            BASE / "data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep2/Hsap_WT_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_WT/rep3/Hsap_WT_rep3_prod.dcd",
        ],
    },
    "4mut": {
        "prmtop": BASE / "data/md_runs/Hsap_4mut/Hsap_4mut.prmtop",
        "trajs": [
            BASE / "data/md_runs/Hsap_4mut/rep1/Hsap_4mut_rep1_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep2/Hsap_4mut_rep2_prod.dcd",
            BASE / "data/md_runs/Hsap_4mut/rep3/Hsap_4mut_rep3_prod.dcd",
        ],
    },
}

# cGAS full-length → PDB resid mapping: PDB = 219 + (FL_res - 200)
def cgas(res):
    return 219 + (res - 200)

MUT_SITES = {
    "D431S": cgas(431),  # PDB 450
    "K479E": cgas(479),  # PDB 498
    "L495Y": cgas(495),  # PDB 514
    "K498T": cgas(498),  # PDB 517
}
IFACE_TARGETS = [cgas(r) for r in range(211, 220)]  # PDB 230-238

CONTACT_CUTOFF = 8.0   # Å, CA-CA distance
CONTACT_FRAC = 0.50    # min fraction for persistent contact
STRIDE = 20            # analyze every Nth frame (100ps × 20 = 2ns)


def get_ca_resids(u):
    """Get CA residue IDs for all protein residues."""
    ca = u.select_atoms("protein and name CA")
    return ca.resids


def compute_contact_map(u, resids, stride=STRIDE):
    """Compute CA-CA contact fraction matrix."""
    ca = u.select_atoms("protein and name CA")
    n_res = len(resids)
    n_frames = len(u.trajectory)
    indices = range(0, n_frames, stride)
    n_sample = len(indices)

    contact_count = np.zeros((n_res, n_res), dtype=np.int32)
    for s, fi in enumerate(indices):
        u.trajectory[fi]
        dist_mat = distance_array(ca.positions, ca.positions)
        contact_count += (dist_mat < CONTACT_CUTOFF).astype(np.int32)

    return contact_count.astype(np.float32) / n_sample


def compute_dccm_from_positions(u, resids, stride=STRIDE):
    """Compute DCCM from CA positions (sampled every stride frames)."""
    ca = u.select_atoms("protein and name CA")
    n_res = len(resids)
    n_frames = len(u.trajectory)
    indices = list(range(0, n_frames, stride))
    n_sample = len(indices)

    # Collect CA positions
    ca_traj = np.zeros((n_sample, n_res, 3), dtype=np.float32)
    for s, fi in enumerate(indices):
        u.trajectory[fi]
        ca_traj[s] = ca.positions

    # Center
    ca_centered = ca_traj - ca_traj.mean(axis=0, keepdims=True)

    # Compute DCCM efficiently: C_ij = <dr_i · dr_j> / sqrt(<dr_i²> <dr_j²>)
    # dr is (n_frames, n_res, 3)
    # cov = mean over frames of dot product
    cov = np.einsum("fid,fjd->ij", ca_centered, ca_centered) / n_sample
    diag = np.sqrt(np.diag(cov))
    # Avoid division by zero
    diag[diag < 1e-10] = 1e-10
    dccm = cov / np.outer(diag, diag)
    return dccm


def build_graph(contact_frac, resids, dccm):
    """Build weighted graph from contacts and correlations."""
    G = nx.Graph()
    n_res = len(resids)
    resid_to_idx = {r: i for i, r in enumerate(resids)}

    for i, r in enumerate(resids):
        G.add_node(r)

    for i in range(n_res):
        for j in range(i + 1, n_res):
            if contact_frac[i, j] < CONTACT_FRAC:
                continue
            c = abs(dccm[i, j])
            weight = -np.log(max(c, 0.001))
            G.add_edge(resids[i], resids[j], weight=weight)

    return G


def find_allosteric_paths(G, sources, targets):
    """Find shortest paths from each source to each target."""
    paths = []
    for src_name, src_res in sources.items():
        if src_res not in G:
            continue
        for tgt_res in targets:
            if tgt_res not in G:
                continue
            try:
                p = nx.shortest_path(G, source=src_res, target=tgt_res, weight="weight")
                pl = nx.shortest_path_length(G, source=src_res, target=tgt_res, weight="weight")
                paths.append((src_name, src_res, tgt_res, p, pl))
            except nx.NetworkXNoPath:
                pass
    return paths


def analyze_rep(u, name):
    """Analyze one replica."""
    print(f"    [{name}] Aligning...", flush=True)
    u.trajectory[0]
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()

    resids = get_ca_resids(u)
    print(f"    [{name}] {len(resids)} CA residues, {len(u.trajectory)} frames")

    contact_frac = compute_contact_map(u, resids)
    print(f"    [{name}] Contacts >50%: {(contact_frac > CONTACT_FRAC).sum() // 2} pairs")

    dccm = compute_dccm_from_positions(u, resids)
    print(f"    [{name}] DCCM computed")

    G = build_graph(contact_frac, resids, dccm)
    print(f"    [{name}] Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    paths = find_allosteric_paths(G, MUT_SITES, IFACE_TARGETS)
    print(f"    [{name}] Found {len(paths)} paths")

    return paths, G, resids


def main():
    print("=" * 60)
    print("P1-5: Allosteric Pathway Identification")
    print("=" * 60)
    print(f"\nSources: {MUT_SITES}")
    print(f"Targets: {IFACE_TARGETS}")

    all_data = {}

    for sys_name, cfg in SYSTEMS.items():
        print(f"\n{'─'*40}\nSystem: {sys_name}\n{'─'*40}")
        all_paths = []
        all_graphs = []

        for rep_idx, traj_path in enumerate(cfg["trajs"]):
            u = mda.Universe(str(cfg["prmtop"]), str(traj_path))
            paths, G, resids = analyze_rep(u, f"{sys_name}_rep{rep_idx+1}")
            all_paths.extend(paths)
            all_graphs.append((G, resids))

        all_data[sys_name] = {"paths": all_paths, "graphs": all_graphs}

    # ── Report ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("RESULTS")
    print(f"{'=' * 60}")

    for sys_name in ["WT", "4mut"]:
        data = all_data[sys_name]
        path_counter = Counter()
        all_res = set()
        for _, src, tgt, p, pl in data["paths"]:
            for r in p[1:-1]:
                path_counter[r] += 1
            all_res.update(p)

        sources_set = set(MUT_SITES.values())
        targets_set = set(IFACE_TARGETS)
        hubs = {r: c for r, c in path_counter.items()
                if r not in sources_set and r not in targets_set}

        print(f"\n{sys_name}:")
        print(f"  Total paths: {len(data['paths'])}")
        print(f"  Unique residues in all paths: {len(all_res)}")
        if data["paths"]:
            lengths = [p[4] for p in data["paths"]]
            print(f"  Path length: mean={np.mean(lengths):.2f}, median={np.median(lengths):.2f}")

        top_hubs = sorted(hubs.items(), key=lambda x: -x[1])[:15]
        print(f"  Top 15 hub residues:")
        for res, count in top_hubs:
            chain = "cGAS" if res >= 219 else "TRIM41"
            print(f"    PDB {res:4d} ({chain}): {count:3d} paths")

    # ── Comparison ──────────────────────────────────────────────────────
    wt_path_counter = Counter()
    for _, src, tgt, p, pl in all_data["WT"]["paths"]:
        for r in p[1:-1]:
            wt_path_counter[r] += 1
    mut_path_counter = Counter()
    for _, src, tgt, p, pl in all_data["4mut"]["paths"]:
        for r in p[1:-1]:
            mut_path_counter[r] += 1

    sources_set = set(MUT_SITES.values())
    targets_set = set(IFACE_TARGETS)

    wt_hubs = {r for r, c in sorted(
        {r: cnt for r, cnt in wt_path_counter.items()
         if r not in sources_set and r not in targets_set
        }.items(), key=lambda x: -x[1])[:30]}
    mut_hubs = {r for r, c in sorted(
        {r: cnt for r, cnt in mut_path_counter.items()
         if r not in sources_set and r not in targets_set
        }.items(), key=lambda x: -x[1])[:30]}

    shared = wt_hubs & mut_hubs
    print(f"\nHub comparison (top 30):")
    print(f"  Shared: {len(shared)}")
    print(f"  WT-only: {len(wt_hubs - mut_hubs)}")
    print(f"  4mut-only: {len(mut_hubs - wt_hubs)}")
    if wt_hubs - mut_hubs:
        print(f"  WT-only hubs: {sorted(wt_hubs - mut_hubs)}")
    if mut_hubs - wt_hubs:
        print(f"  4mut-only hubs: {sorted(mut_hubs - wt_hubs)}")

    # ── WT vs 4mut path architecture ────────────────────────────────────
    print(f"\nPer-mutation consensus pathways:")
    for mut_name, mut_res in MUT_SITES.items():
        for sys_name in ["WT", "4mut"]:
            site_paths = [(p, pl) for n, s, t, p, pl in all_data[sys_name]["paths"] if n == mut_name]
            if not site_paths:
                print(f"  {mut_name} ({sys_name}): no paths")
                continue
            # Count intermediate residues
            pc = Counter()
            for p, pl in site_paths:
                for r in p[1:-1]:
                    pc[r] += 1
            n = len(site_paths)
            consensus = [(r, c) for r, c in pc.most_common(10) if c >= n * 0.5]
            if consensus:
                print(f"  {mut_name} ({sys_name}, {n} paths): "
                      f"consensus={[(r, f'{c}/{n}') for r, c in consensus[:5]]}")
            else:
                print(f"  {mut_name} ({sys_name}, {n} paths): no consensus (>50%)")

    # ── Plot: hub comparison ────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, sys_name in enumerate(["WT", "4mut"]):
        ax = axes[ax_idx]
        pc = Counter()
        for _, src, tgt, p, pl in all_data[sys_name]["paths"]:
            for r in p[1:-1]:
                pc[r] += 1
        sources_set = set(MUT_SITES.values())
        targets_set = set(IFACE_TARGETS)
        filtered = {r: c for r, c in pc.items()
                     if r not in sources_set and r not in targets_set and c > 0}

        if filtered:
            res_list = sorted(filtered.keys())
            counts = [filtered[r] for r in res_list]
            colors = ["steelblue" if r < 219 else "coral" for r in res_list]
            ax.bar(range(len(res_list)), counts, color=colors, alpha=0.8)
            ax.set_xlabel("Residue (sorted by PDB resid)")
            ax.set_ylabel("Path frequency")
            ax.set_title(f"{sys_name}: Pathway hub residues")
            # Label top 5
            for res, count in sorted(filtered.items(), key=lambda x: -x[1])[:5]:
                idx = res_list.index(res)
                ax.annotate(f"{res}", (idx, count), fontsize=6, ha="center",
                            textcoords="offset points", xytext=(0, 5))
        else:
            ax.set_title(f"{sys_name}: No paths found")

    fig.tight_layout()
    fig.savefig(OUTDIR / "pathway_hubs.png", dpi=300)
    plt.close(fig)
    print(f"\nSaved: {OUTDIR / 'pathway_hubs.png'}")

    # ── Export JSON ─────────────────────────────────────────────────────
    export = {}
    for sys_name in ["WT", "4mut"]:
        pc = Counter()
        for _, src, tgt, p, pl in all_data[sys_name]["paths"]:
            for r in p[1:-1]:
                pc[r] += 1
        sources_set = set(MUT_SITES.values())
        targets_set = set(IFACE_TARGETS)
        hubs = {str(r): c for r, c in pc.most_common(50)
                if r not in sources_set and r not in targets_set}
        export[sys_name] = {
            "n_paths": len(all_data[sys_name]["paths"]),
            "hub_residues": hubs,
        }
    with open(OUTDIR / "pathway_results.json", "w") as f:
        json.dump(export, f, indent=2)
    print(f"Saved: {OUTDIR / 'pathway_results.json'}")

    print(f"\n{'=' * 60}\nDone.\n{'=' * 60}")


if __name__ == "__main__":
    main()
