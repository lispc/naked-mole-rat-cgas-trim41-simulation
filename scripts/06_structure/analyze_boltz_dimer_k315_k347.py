import gemmi
import numpy as np
from pathlib import Path

CIF_DIR = Path("data/boltz_cgas_dna_dimer_trim41/boltz_results_boltz_cgas_dna_dimer_trim41/predictions/boltz_cgas_dna_dimer_trim41")

def get_ca_pos(chain, res_idx_1based):
    for atom in chain[res_idx_1based - 1]:
        if atom.name == "CA":
            return np.array([atom.pos.x, atom.pos.y, atom.pos.z])
    return None

def get_nz_pos(chain, res_idx_1based):
    for atom in chain[res_idx_1based - 1]:
        if atom.name == "NZ":
            return np.array([atom.pos.x, atom.pos.y, atom.pos.z])
    return None

def chain_heavy_atoms(chain):
    pos = []
    for res in chain:
        for atom in res:
            if atom.element.name != "H":
                pos.append([atom.pos.x, atom.pos.y, atom.pos.z])
    return np.array(pos)

print("=" * 95)
print("Boltz-2 cGAS-DNA Dimer + TRIM41 SPRY: K315(pos225) / K@293 Analysis")
print("=" * 95)
print("NOTE: Human K347 maps to CIF pos 257 = LEU (not LYS in this construct).")
print("      CIF pos 293 is a construct-specific LYS not corresponding to human K347.")
print("      Human equivalent of CIF pos 293 ≈ K383, but YAML seq differs from Hsap.")

results = []

for model_idx in range(5):
    cif_path = CIF_DIR / f"boltz_cgas_dna_dimer_trim41_model_{model_idx}.cif"
    st = gemmi.read_structure(str(cif_path))
    
    cgas1 = st[0]["A"]
    cgas2 = st[0]["B"]
    spry = st[0]["C"]
    
    spry_heavy = chain_heavy_atoms(spry)
    spry_ca = []
    for res in spry:
        for atom in res:
            if atom.name == "CA":
                spry_ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
    spry_com = np.array(spry_ca).mean(axis=0)
    
    print(f"\n--- Model {model_idx} ---")
    
    model_data = {"model": model_idx}
    
    for chain_label, cgas in [("cgas1 (A)", cgas1), ("cgas2 (B)", cgas2)]:
        k315_nz = get_nz_pos(cgas, 225)      # K315 (SASKMLSKFRK motif)
        k293_nz = get_nz_pos(cgas, 293)      # K@293 (construct-specific; NOT human K347)
        k315_ca = get_ca_pos(cgas, 225)
        k293_ca = get_ca_pos(cgas, 293)

        k315_d_com = np.linalg.norm(k315_nz - spry_com) if k315_nz is not None else None
        k315_d_min = np.min(np.linalg.norm(spry_heavy - k315_nz, axis=1)) if k315_nz is not None else None
        k293_d_com = np.linalg.norm(k293_nz - spry_com) if k293_nz is not None else None
        k293_d_min = np.min(np.linalg.norm(spry_heavy - k293_nz, axis=1)) if k293_nz is not None else None
        k315_k293_ca = np.linalg.norm(k315_ca - k293_ca) if (k315_ca is not None and k293_ca is not None) else None

        print(f"\n  {chain_label}:")
        print(f"    K315 (pos 225):          SPRY COM={k315_d_com:.1f}A, min={k315_d_min:.1f}A")
        print(f"    K@293 (NOT human K347):  SPRY COM={k293_d_com:.1f}A, min={k293_d_min:.1f}A")
        print(f"    K315-K@293 CA distance: {k315_k293_ca:.1f}A")

        # All LYS distances
        lys_info = []
        for i, res in enumerate(cgas):
            if res.name == "LYS":
                nz = get_nz_pos(cgas, i+1)
                if nz is not None:
                    d_com = np.linalg.norm(nz - spry_com)
                    d_min = np.min(np.linalg.norm(spry_heavy - nz, axis=1))
                    lys_info.append((i+1, d_com, d_min))
        lys_info.sort(key=lambda x: x[1])

        print(f"    Top 5 nearest LYS to SPRY COM:")
        for resid, d_com, d_min in lys_info[:5]:
            marker = ""
            if resid == 225: marker = " <- K315"
            elif resid == 293: marker = " <- K@293 (NOT K347)"
            print(f"      resid {resid:3d}: {d_com:5.1f}A COM, {d_min:5.1f}A min{marker}")

        prefix = "cgas1" if "cgas1" in chain_label else "cgas2"
        model_data[f"{prefix}_k315_com"] = k315_d_com
        model_data[f"{prefix}_k315_min"] = k315_d_min
        model_data[f"{prefix}_k293_com"] = k293_d_com
        model_data[f"{prefix}_k293_min"] = k293_d_min
        model_data[f"{prefix}_k315_k293_ca"] = k315_k293_ca

    # Inter-chain K@293 distance
    cgas1_k293_nz = get_nz_pos(cgas1, 293)
    cgas2_k293_nz = get_nz_pos(cgas2, 293)

    print(f"\n  Cross-chain (SPRY -> K@293 on OTHER chain):")
    if cgas2_k293_nz is not None:
        d_com = np.linalg.norm(cgas2_k293_nz - spry_com)
        d_min = np.min(np.linalg.norm(spry_heavy - cgas2_k293_nz, axis=1))
        print(f"    cgas2 K@293 -> SPRY: COM={d_com:.1f}A, min={d_min:.1f}A")
        model_data["cross_k293_cgas2_min"] = d_min
    if cgas1_k293_nz is not None:
        d_com = np.linalg.norm(cgas1_k293_nz - spry_com)
        d_min = np.min(np.linalg.norm(spry_heavy - cgas1_k293_nz, axis=1))
        print(f"    cgas1 K@293 -> SPRY: COM={d_com:.1f}A, min={d_min:.1f}A")
        model_data["cross_k293_cgas1_min"] = d_min

    # Inter-chain K@293 distance
    if cgas1_k293_nz is not None and cgas2_k293_nz is not None:
        d_inter = np.linalg.norm(cgas1_k293_nz - cgas2_k293_nz)
        print(f"    K@293(cgas1) <-> K@293(cgas2) NZ: {d_inter:.1f}A")
        model_data["k293_inter_nz"] = d_inter
    
    results.append(model_data)

# Summary table
print("\n" + "=" * 95)
print("SUMMARY TABLE (NOTE: K@293 is NOT human K347; human K347=Leu in this construct)")
print("=" * 95)
print(f"{'Model':>6} | {'K315 min':>10} | {'K@293 min':>10} | {'K315-K@293':>10} | {'cross K@293':>10} | {'K@293-K@293':>10}")
print("       |   min (A)  |   min (A)  |   CA (A)   |   min (A)  |   NZ (A)   ")
print("-" * 95)
for r in results:
    k315_min = r.get("cgas1_k315_min", 0)
    k293_min = r.get("cgas1_k293_min", 0)
    k315_k293 = r.get("cgas1_k315_k293_ca", 0)
    cross = min(r.get("cross_k293_cgas1_min", 999), r.get("cross_k293_cgas2_min", 999))
    k293_inter = r.get("k293_inter_nz", 0)
    print(f"{r['model']:>6} | {k315_min:>10.1f} | {k293_min:>10.1f} | {k315_k293:>10.1f} | {cross:>10.1f} | {k293_inter:>10.1f}")

print("\nNOTES:")
print("  - K315 (CIF pos 225) = human K315 (verified by SASKMLSKFRK motif)")
print("  - K@293 (CIF pos 293) is a construct-specific LYS; NOT human K347")
print("  - Human K347 maps to CIF pos 257, which is LEU in this YAML construct")
print("  - Human equivalent of CIF pos 293 = pos 383, but YAML seq diverges from Hsap")
print("  - Catalytic range for ubiquitination: NZ -> Ub G76 C < 8-10 A typically")
print("=" * 95)
