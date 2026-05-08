#!/usr/bin/env python3
"""Convert Boltz-2 CIF predictions to PDB and check K315/K347 positions."""
import gemmi
import numpy as np
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent.parent
CIF_DIR = BASE / "data/boltz_cgas_dna_trim41/boltz_results_boltz_cgas_dna_trim41/predictions/boltz_cgas_dna_trim41"

for model_idx in range(5):
    cif_path = CIF_DIR / f"boltz_cgas_dna_trim41_model_{model_idx}.cif"
    pdb_path = CIF_DIR / f"boltz_cgas_dna_trim41_model_{model_idx}.pdb"
    if not cif_path.exists():
        print(f"Model {model_idx}: CIF not found")
        continue
    try:
        st = gemmi.read_structure(str(cif_path))
        st.write_pdb(str(pdb_path))
        print(f"Model {model_idx}: {len(st[0])} chains")
        for ch in st[0]:
            residues = [r for r in ch if r.seqid.num > 0]
            if residues:
                print(f"  Chain {ch.name}: {len(residues)} residues")
        # Check for K315 and K347 (need to infer chain assignments)
        # For now just report the conversion
    except Exception as e:
        print(f"Model {model_idx}: conversion failed - {e}")

print("Done.")
