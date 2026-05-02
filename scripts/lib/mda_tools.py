"""Common MDAnalysis helper functions."""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import distance_array


def load_universe(prmtop, trajectory):
    """Load MDAnalysis Universe, guessing elements if missing."""
    u = mda.Universe(str(prmtop), str(trajectory))
    if not hasattr(u.atoms, 'elements') or u.atoms.elements is None or len(u.atoms.elements) == 0:
        u.atoms.guess_elements()
    return u


def align_trajectory(u, ref, selection="protein and name CA"):
    """In-memory alignment of trajectory to reference."""
    aligner = align.AlignTraj(u, ref, select=selection, in_memory=True)
    aligner.run()
    return aligner


def get_com_distance(u, sel1, sel2):
    """Return COM distance between two selections for current frame."""
    a1 = u.select_atoms(sel1)
    a2 = u.select_atoms(sel2)
    return np.linalg.norm(a1.center_of_mass() - a2.center_of_mass())


def count_hbonds_simple(u, sel1, sel2, distance=3.5):
    """
    Count heavy-atom donor-acceptor pairs between two selections.
    Simplified criterion: N/O distance < cutoff (no angle check).
    Returns array of counts per frame.
    """
    heavy1 = u.select_atoms(f"({sel1}) and (name N* or name O*)")
    heavy2 = u.select_atoms(f"({sel2}) and (name N* or name O*)")
    n_frames = len(u.trajectory)
    counts = np.zeros(n_frames, dtype=int)
    for i, ts in enumerate(u.trajectory):
        dist_mat = distance_array(heavy1.positions, heavy2.positions)
        counts[i] = np.sum(dist_mat < distance) // 2
    return counts


def get_interface_residues(u, sel1, sel2, cutoff=5.0):
    """Return sets of interface residue IDs within cutoff (Å)."""
    a1 = u.select_atoms(sel1)
    a2 = u.select_atoms(sel2)
    u.trajectory[0]
    dist_mat = distance_array(a1.positions, a2.positions)
    close = np.argwhere(dist_mat < cutoff)
    resids1 = set(a1.resids[close[:, 0]])
    resids2 = set(a2.resids[close[:, 1]])
    return resids1, resids2
