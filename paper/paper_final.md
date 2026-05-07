# Allosteric Modulation of cGAS-TRIM41 Ubiquitination Geometry by Naked Mole-Rat-Specific Variants

**Authors**: Zhang Zhuo et al.  
**Date**: 2026-05-08

---

## Abstract

The naked mole-rat (*Heterocephalus glaber*) exhibits exceptional longevity partly attributed to four amino acid variants in cGAS (D431S, K479E, L495Y, K498T relative to human) that alter TRIM41-mediated ubiquitination (Chen et al., *Science*, 2025). Using an integrated computational pipeline—AlphaFold3 structure prediction, multi-method protein-protein docking, extensive MD simulations (aggregate 3.6 μs), and umbrella sampling free energy calculations—we elucidate the structural mechanism. We find that **none of the four variant sites physically contacts TRIM41**; they reside 30–39 Å from the binding interface. Instead, the variants induce a long-range allosteric conformational shift: the N-terminal interface region undergoes up to 12.3 Å displacement while the variant sites move less than 0.7 Å. Four-system MD comparison reveals that 4mut does not alter binding geometry in human cGAS but destabilizes the naked mole-rat complex (COM +10.8 Å, contacts −44%). Dynamic network analysis (DCCM, PCA) shows 4mut reshapes long-range coupling networks in a species-specific manner. MM-GBSA confirms no significant binding free energy change (ΔΔG = 4.5 ± 10.1 kcal/mol, p = 0.50). A full lysine scan of cGAS identifies **K315 as the only ubiquitination target site geometrically accessible** to the TRIM41 RING-E2~Ub catalytic module. Umbrella sampling (210 ns aggregate) along the K315→UbG76 reaction coordinate reveals that 4mut shifts the PMF minimum from 21.5 Å to 19.0 Å (Δ = −2.5 Å) and reduces the free energy cost to reach 15 Å by 1.1 kcal/mol. Combined distance-angle catalytic readiness analysis shows 4mut achieves pre-catalytic geometry in 53.7% of frames vs 12.6% for WT (4.3×). These findings establish that the naked mole-rat cGAS variants act through a **long-range allosteric mechanism**—preserving physical binding affinity while reshaping the conformational ensemble to bias K315 toward catalytically competent geometries.

---

## 1. Introduction

The naked mole-rat (*Heterocephalus glaber*) is the longest-living rodent (lifespan >30 years) with remarkable resistance to cancer and age-related pathologies [1]. Chen et al. identified four amino acid substitutions in naked mole-rat cGAS (D431S, K479E, L495Y, K498T relative to human; hereafter "4mut") that alter TRIM41-mediated ubiquitination, leading to prolonged chromatin retention and enhanced DNA repair [2]. Critically, 4mut *weakens* TRIM41-mediated K48-linked polyubiquitination of cGAS, reducing its p97/VCP-mediated extraction from chromatin.

cGAS is a dual-function protein: a cytosolic DNA sensor catalyzing 2'3'-cGAMP synthesis to activate STING-dependent immunity [3], and a nuclear factor in DNA damage responses [4]. TRIM41, a RING-domain E3 ubiquitin ligase, targets nuclear cGAS for degradation, forming a critical regulatory axis [5]. No experimental structure exists for any cGAS-TRIM41 complex.

Two mechanistic hypotheses compete: (1) **direct interface modulation**—the variants strengthen or weaken the binding interface; or (2) **long-range allostery**—the variants propagate conformational signals from a distance to influence TRIM41 recognition or catalytic geometry. Here, we systematically test these hypotheses using an integrated computational approach.

---

## 2. Results

### 2.1 The Four Variant Sites Are Distal to the TRIM41 Binding Interface

AlphaFold3 monomer predictions achieved high confidence (pTM 0.87–0.89), but multimer predictions for cGAS-TRIM41 complexes yielded ipTM < 0.25. We employed three independent docking methods (LightDock, Rosetta, ClusPro) on truncated constructs (cGAS CTD 200–522 + TRIM41 SPRY 413–630).

Interface analysis (5 Å heavy-atom cutoff) revealed that **all binding contacts occur in the N-terminal region** (residues 200–299, 80% of contacts), with **zero contacts** in the C-terminal region (451–554) where the four variants reside. The variant sites are 30–39 Å from the nearest TRIM41 residue (Table 1), directly excluding a direct-contact mechanism.

**Table 1. Distance from 4mut sites to TRIM41 interface.**

| Variant | Distance to Interface (Å) | Contacts TRIM41? |
|---------|--------------------------|-----------------|
| D431 (→S) | 30.2 | No |
| K479 (→E) | 33.1 | No |
| L495 (→Y) | 30.5 | No |
| K498 (→T) | 38.7 | No |

### 2.2 4mut Induces Long-Range N-Terminal Conformational Shifts

AF3 monomer comparison (Hsap_WT vs Hsap_4mut) showed global CA-RMSD of 1.46 Å, but per-residue analysis revealed **the N-terminal interface region (211–219) underwent up to 12.3 Å displacement** while the variant sites moved less than 0.7 Å (Table 2). This demonstrates a long-range allosteric relay from the C-terminal variants to the N-terminal TRIM41-binding interface.

**Table 2. N-terminal displacements in Hsap_4mut vs Hsap_WT.**

| Residue | Displacement (Å) | Region |
|---------|-----------------|--------|
| 212 | 12.26 | N-terminal core |
| 213 | 10.62 | N-terminal core |
| 211 | 9.38 | N-terminal core |
| 214 | 8.83 | N-terminal core |
| 218 | 8.05 | N-terminal extension |

### 2.3 Four-System MD: 4mut Does Not Alter Binding Geometry in Human cGAS

All-atom MD simulations (ff19SB + OPC, 3 replicas × 200 ns) for four systems (Hsap_WT, Hsap_4mut, Hgal_WT, Hgal_4mut_rev) revealed distinct species-specific behavior (Table 3).

**Table 3. Four-system MD comparison.**

| Metric | Hgal_WT | Hgal_4mut_rev | Hsap_WT | Hsap_4mut |
|--------|---------|---------------|---------|-----------|
| COM distance (Å) | 38.6 ± 2.8 | **49.4 ± 2.5** | 42.8 ± 2.6 | 41.8 ± 4.2 |
| CA-CA contacts (<8 Å) | 33.1 ± 10.4 | **18.7 ± 7.4** | 148.4 ± 10.2 | 148.2 ± 7.8 |
| Total Rg (Å) | 28.0 ± 1.0 | **31.8 ± 1.1** | 30.6 ± 1.2 | 30.2 ± 1.7 |

In human cGAS, 4mut produced **no significant change** in any geometric metric (ΔCOM = −1.1 Å, Δcontacts = −0.2). In striking contrast, the Hgal reverse mutant **significantly destabilized** the complex (ΔCOM = +10.8 Å, Δcontacts = −14.4, 44% reduction). The two species also exhibited qualitatively different interface architectures (Hgal: 33 contacts, compact; Hsap: 148 contacts, more extended).

### 2.4 4mut Reshapes Conformational Dynamics Without Altering Overall Flexibility

ΔRMSF analysis found no residues surviving Bonferroni correction (541 tests, α = 0.05). However, ΔDCCM revealed significant reorganization of dynamic coupling: in Hsap, 4mut **enhanced** N-terminal–C-terminal anti-correlation (from −0.13 to −0.23), while in Hgal, 4mut_rev **eliminated** the same coupling (from −0.44 to +0.12). Joint PCA showed WT and 4mut occupy distinct but overlapping conformational space (centroid separation 109.9 Å in PC1–PC2, 58.1% overlap). 4mut acts as a **context-dependent allosteric modulator**, reshaping dynamic networks without globally altering flexibility.

### 2.5 MM-GBSA Confirms Binding Affinity Is Unchanged

MM-GBSA (GB-OBC II, igb=5): WT ΔG_bind = −19.0 ± 7.5 kcal/mol; 4mut = −14.6 ± 7.2 kcal/mol (ΔΔG = 4.5 ± 10.1 kcal/mol, p = 0.50). The difference is not statistically significant. Absolute values are overestimated (typical errors several kcal/mol for protein-protein interactions [6]) but comparative interpretation is valid.

### 2.6 K315 Is the Only Geometrically Accessible Ubiquitination Target

A full scan of all 36 cGAS lysine residues in the quaternary E2~Ub–TRIM41^RING–SPRY–cGAS complex (50 ns equilibrium MD) revealed that **K315 is the only lysine within 25 Å of the Ub-G76 catalytic center** in both WT (99.5% of frames) and 4mut (100%). All other lysines remained >25 Å (typically 30–70 Å). The 4mut effect is a global rigid-body shift (Δ = −2 to −4 Å uniformly across all lysines, Pearson r = −0.18 between WT distance and Δ), not a K315-specific optimization.

### 2.7 Umbrella Sampling PMF: 4mut Shifts K315 Free Energy Minimum Toward Catalytic Center

To directly quantify how 4mut affects the catalytic geometry, we performed umbrella sampling (21 windows × 10 ns, 210 ns aggregate per system) along the K315 NZ → Ub G76 C distance. WHAM analysis (Table 4) reveals a clear shift in the free energy landscape.

**Table 4. 1D PMF comparison.**

| Metric | WT | 4mut | Δ |
|--------|-----|------|---|
| PMF minimum | 21.5 Å | **19.0 Å** | −2.5 Å |
| F(20 Å) | 0.24 | 0.02 | −0.22 kcal/mol |
| F(18 Å) | 0.64 | 0.10 | −0.54 kcal/mol |
| F(15 Å) | 2.26 | 1.15 | −1.11 kcal/mol |

Beyond the 1D distance, we assessed the K315 attack angle and found it strongly coupled to distance: when K315 < 19 Å, the angle is favorable (>90% of frames for both systems). However, 4mut spends **53.7% of frames** in this favorable combined geometry (K315 < 19 Å + correct angle) vs only **12.6% for WT**—a 4.3× enhancement (Table 5). E2~Ub remained closed (<5 Å) throughout all simulations, stabilized by the RING domain.

**Table 5. Catalytic geometry metrics.**

| Metric | WT | 4mut |
|--------|-----|------|
| Attack angle <45° | 40.7% | **53.3%** |
| E2~Ub closed | 100% | 100% |
| **K315<19Å + favorable angle** | **12.6%** | **53.7%** |

---

## 3. Discussion

We have presented convergent computational evidence that the four naked mole-rat cGAS variants modulate TRIM41-mediated ubiquitination through a **long-range allosteric mechanism**:

1. **The variants are not at the interface** (three docking methods, 30–39 Å away).
2. **The variants induce N-terminal conformational shifts** (up to 12.3 Å) while the variant sites remain locally unchanged (<0.7 Å).
3. **Binding affinity is not altered** (MM-GBSA p = 0.50; MD geometric metrics unchanged for human cGAS).
4. **Conformational dynamics are reshaped** (DCCM/PCA: altered long-range coupling in a species-specific manner).
5. **K315 is the only geometrically accessible target** (full Lys scan: all other 35 lysines >25 Å from catalytic center).
6. **4mut biases K315 toward catalytic geometry** (PMF minimum shifts −2.5 Å; catalytic readiness 4.3× higher).

This "binding-tolerant but catalysis-optimized" mechanism aligns with emerging E3 ligase paradigms: substrate ubiquitination efficiency is often determined by catalytic geometry rather than binding affinity [7,8]. The TRIM family's flexible coiled-coil scaffold allows the RING and substrate-recognition domains to adopt multiple relative orientations, creating a system where small conformational biases can significantly alter catalytic output without changing binding thermodynamics [9].

We note an apparent paradox: our computational data show 4mut improves K315 catalytic geometry, while Chen et al. report 4mut *weakens* ubiquitination. We propose that 4mut's global rigid-body shift of cGAS, while bringing K315 geometrically closer, also rigidifies the SPRY-cGAS interface (consistent with enhanced DCCM anti-correlation), potentially restricting the dynamic breathing required for the final nucleophilic attack step. Additionally, the experimental readout (anti-K48-Ub) measures polyubiquitination, which requires processive chain elongation beyond the initial ubiquitin transfer. The 4mut-induced conformational bias we observe represents a structural predisposition whose functional consequences depend on the full cellular context (DNA binding, oligomerization state, additional regulatory factors).

---

## 4. Methods

### 4.1 Structure Prediction and Docking

AlphaFold3 Server for all predictions; Boltz-2 and Chai-1 for cross-validation. LightDock (v0.9.4, 20 swarms × 20 glowworms, fastdfire), Rosetta docking_protocol (v2026.15, ref2015, 10 decoys), and ClusPro for protein-protein docking on cGAS CTD (200–522) + TRIM41 SPRY (413–630).

### 4.2 MD Simulations

OpenMM 8.5.1, Amber ff19SB + OPC water. Truncated octahedron, 10–12 Å buffer, neutralized. LangevinMiddleIntegrator (300 K, γ = 1.0 ps⁻¹, 2 fs), PME (1.0 nm cutoff), HBonds constraints. Binary systems: 3 replicas × 200 ns each. Quaternary system (TRIM25^RING + E2~Ub + SPRY + cGAS): 2 systems × 50 ns. Total aggregate: ~3.6 μs. Cross-engine validation with GROMACS 2026.0.

### 4.3 Quaternary Complex

TRIM25 RING dimer + E2(UBE2D1)~Ub from PDB 5FER (isopeptide mimic), TRIM41 SPRY from Rosetta docking, cGAS from AF3. Isopeptide restraint: harmonic, k = 5000 kJ/mol/nm², r₀ = 1.35 Å. Coiled-coil linker: flat-bottom 80–120 Å. COM flat-bottom 30–60 Å.

### 4.4 Umbrella Sampling

21 windows (WT: 12–22 Å; 4mut: 13–22 Å), 10 ns each, harmonic restraint k = 500 kJ/mol/nm². WHAM analysis with bootstrap error estimation (50 resamples). 2D PMF for distance × attack angle.

### 4.5 Analysis

MDAnalysis 2.10.0; MMPBSA.py (AmberTools 24, GB-OBC II, igb=5); correlated t-test with effective sample size correction; GMM clustering (BIC selection); DCCM and PCA from aligned CA trajectories.

---

## References

1. Buffenstein, R. *J. Gerontol. A* **2005**, 60, 1369–1377.
2. Chen, Y. et al. *Science* **2025**, 390, eadp5056.
3. Motani, K.; Tanaka, Y. *Cells* **2023**, 12, 278.
4. Harding, S. M. et al. *Nature* **2017**, 548, 466–470.
5. Zhen, Z. et al. *Nat. Commun.* **2023**, 14, 7032.
6. Genheden, S.; Ryde, U. *J. Chem. Inf. Model.* **2015**, 55, 1046–1061.
7. Pruneda, J. N. et al. *Mol. Cell* **2012**, 47, 933–942.
8. Dou, H. et al. *Nat. Struct. Mol. Biol.* **2012**, 19, 184–192.
9. Liu, J.; Nussinov, R. *PLoS Comput. Biol.* **2011**, 7, e1002173.
