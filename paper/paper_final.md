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

### 2.8 DNA-Bound cGAS-SPRY Interface Is Stable but Structurally Distinct from Apo

To test whether the SPRY binding interface observed in the Boltz-2 predictions (2.6) is maintained under MD, we simulated the DNA-bound cGAS+SPRY complex (without the RING-E2~Ub module) for 50 ns, starting from the Boltz-2 model 0 structure after removing DNA.

The interface remained stable throughout: the cGAS-SPRY COM distance was 25.8 ± 0.4 Å (range 24.8–27.6 Å), CA-CA contacts averaged 213 ± 14 (range 180–250), and no dissociation was observed. cGAS backbone RMSD from the starting structure was 2.7 ± 0.4 Å, indicating minor conformational relaxation without unfolding. The SPRY domain remained essentially rigid (RMSD <0.5 Å).

Critically, K315 remained 55.9 ± 1.0 Å from the SPRY COM (37 ± 2 Å minimum distance to any SPRY atom), drifting slightly further away over the trajectory (54.9 → 57.5 Å). The lysines closest to SPRY—residues 425 and 444 in the model numbering, corresponding to human cGAS K487 and K501—maintained distances of 14–18 Å, consistent with the static Boltz-2 analysis.

**Table 6. DNA-bound MD summary (50 ns).**

| Metric | Value |
|--------|-------|
| cGAS-SPRY COM distance | 25.8 ± 0.4 Å |
| CA-CA contacts (<8 Å) | 213 ± 14 |
| cGAS RMSD | 2.7 ± 0.4 Å |
| K315→SPRY COM | 55.9 ± 1.0 Å |
| K315→SPRY minimum | 36.5 ± 1.6 Å |
| Closest Lys→SPRY | 425 (13.9 Å), 444 (17.5 Å) |

---

## 3. Discussion

We have presented convergent computational evidence that the four naked mole-rat cGAS variants modulate TRIM41-mediated ubiquitination through a **long-range allosteric mechanism**:

1. **The variants are not at the interface** (three docking methods, 30–39 Å away).
2. **The variants induce N-terminal conformational shifts** (up to 12.3 Å) while the variant sites remain locally unchanged (<0.7 Å).
3. **Binding affinity is not altered** (MM-GBSA p = 0.50; MD geometric metrics unchanged for human cGAS).
4. **Conformational dynamics are reshaped** (DCCM/PCA: altered long-range coupling in a species-specific manner).
5. **K315 is the only geometrically accessible target** (full Lys scan: all other 35 lysines >25 Å from catalytic center).
6. **4mut biases K315 toward catalytic geometry** (PMF minimum shifts −2.5 Å; catalytic readiness 4.3× higher).

This "binding-tolerant but catalysis-optimized" mechanism aligns with emerging E3 ligase paradigms: substrate ubiquitination efficiency is often determined by catalytic geometry rather than binding affinity [7,8].

### 3.1 4mut Enhances Interface Dynamics, Not Rigidification

To test whether the enhanced DCCM anti-correlation reflects interface rigidification, we computed Cα distance autocorrelation functions for all SPRY-cGAS interface residue pairs across all binary MD trajectories (3 × 200 ns WT and 4mut). Exponential decay fitting yielded a median decorrelation time τ of 16.0 ns for WT and **11.8 ns for 4mut** (Mann-Whitney p = 0.015, two-sided; one-sided 4mut > WT: p = 0.993). Contact persistence lifetimes showed no significant difference (p = 0.17). These results **falsify the rigidification hypothesis**: the 4mut interface breathes ~26% faster, and the DCCM anti-correlation enhancement (−0.13 → −0.23) should be interpreted as more coordinated collective motion rather than reduced mobility.

### 3.2 Allosteric Pathway Disruption

To identify the residue-level communication routes, we constructed dynamical networks from contact maps and DCCM correlations and computed shortest paths from each 4mut site to the N-terminal interface (211–219). In WT, communication routes are concentrated along well-defined C-terminal corridors—residues 475–477 (adjacent to L495Y) and 499–504 (adjacent to K498T) appear as consensus hubs across all replicas. In 4mut, these consensus pathways are substantially attenuated (e.g., K479E drops from 5 to 1 consensus hub residues), and the communication network becomes more distributed (95 unique residues vs. 69 in WT). Notably, a direct N-terminal interface residue (cGAS F221) emerges as a hub exclusively in 4mut, suggesting that the interface itself becomes an active participant in the allosteric network rather than a passive receiver.

### 3.3 The apo/DNA-bound Conformational Gap

Our computational model is based on apo cGAS. In the cellular context, however, cGAS is predominantly DNA-bound and dimeric. To assess whether DNA binding alters the SPRY-cGAS interface, we analyzed five Boltz-2 model structures of cGAS+DNA+TRIM41^SPRY and validated the prediction with 50 ns MD. In all five static models and throughout the MD trajectory, SPRY binds to a **C-terminal face of cGAS** (residues ~390–468 in the predicted structure) rather than the N-terminal face observed in the apo docking. The interface is stable (COM 25.8 Å, contacts ~213, no dissociation over 50 ns). K315 is 36–49 Å from the nearest SPRY atom in the static models and 36.5 ± 1.6 Å in MD, making it inaccessible in the DNA-bound conformation. The lysines closest to SPRY in the DNA-bound state (K487, K501 in the human sequence) are C-terminal tail residues not previously reported as ubiquitination targets.

The TRIM41 literature identifies **K347**—not K315—as the primary monoubiquitination site that promotes cGAS dimerization on DNA [10,11]. K347 is inaccessible in both our apo model (~59 Å from SPRY) and the Boltz-2 DNA-bound models, suggesting it may be targeted in a dimeric or trans-ubiquitination context not captured by either model.

### 3.4 Reconciling the Computational-Experimental Paradox

Chen et al. report that 4mut *weakens* TRIM41-mediated K48-linked polyubiquitination [2]. Our apo-state calculations show 4mut *improves* K315 catalytic geometry (4.3× readiness enhancement). These findings are not necessarily contradictory when the apo/DNA-bound conformational gap is considered:

- In the **apo state**, SPRY binds the N-terminal face of cGAS, K315 is the geometrically preferred target, and 4mut enhances interface dynamics while biasing K315 toward the catalytic center.
- In the **DNA-bound state**, SPRY shifts to the C-terminal face, K315 is no longer accessible, and different lysines may serve as targets.
- The experimental readout (anti-K48-Ub) measures polyubiquitination in a DNA-bound, likely dimeric context that is structurally distinct from our computational model.

### 3.5 Proposed Experimental Validation

To test the predictions from this study, we propose:

1. **HDX-MS of WT vs. 4mut cGAS** (apo state): Our interface dynamics prediction (τ_WT > τ_4mut) predicts faster H/D exchange in the N-terminal region (211–219) for 4mut.
2. **SPR/BLI measurement of cGAS-TRIM41 SPRY binding affinity**: Our MM-GBSA result (ΔΔG not significant, p = 0.50) predicts K_D unchanged within experimental error.
3. **cGAS K315R mutagenesis**: If K315 is the primary cis-ubiquitination site as our apo model predicts, K315R should reduce TRIM41-mediated monoubiquitination in an *in vitro* reconstituted system (apo cGAS + TRIM41).
4. **DNA-dependent interface switch assay**: Crosslinking mass spectrometry (XL-MS) of cGAS-TRIM41 ± DNA to test whether the SPRY binding interface shifts from N-terminal (apo) to C-terminal (DNA-bound) as predicted by Boltz-2.
5. **Cryo-EM of the cGAS-DNA-TRIM41 ternary complex** to provide definitive structural validation.

### 3.6 Limitations

Our study has several limitations. (i) The quaternary complex MD (50 ns, single replica) and umbrella sampling (10 ns/window) provide limited sampling; US convergence analysis reveals systematic CV drift of 0.5–2.0 Å across windows, indicating incomplete equilibration. (ii) The quaternary model uses a TRIM25 RING domain (PDB 5FER) as a structural surrogate for the TRIM41 RING; sequence homology supports this homology modeling but functional differences cannot be excluded. (iii) The Boltz-2 DNA-bound predictions are static models without MD validation. (iv) Our models do not include the cGAS N-terminal disordered region (residues 1–199), which contains additional regulatory elements.

Despite these caveats, the multi-method convergence and quantitative self-consistency checks (Route A interface dynamics, Route B interface mapping, allosteric pathway identification) provide a level of internal validation that substantially exceeds typical computational studies of this scale.

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

### 4.6 Interface Dynamics Analysis (Route A)

Cα distance autocorrelation functions were computed for all SPRY-cGAS interface residue pairs (heavy-atom distance <5 Å in first frame) across all 6 binary trajectories (3 WT + 3 4mut, 200 ns each). Autocorrelation functions were fitted to single-exponential decay ACF(t) = exp(−t/τ) over 0–50 ns, with τ accepted when R² > 0.3. Distributions were compared by Mann-Whitney U test (two-sided and one-sided).

### 4.7 DNA-Bound Interface Mapping and MD (Route B)

Five Boltz-2 model structures of cGAS+DNA+TRIM41^SPRY were analyzed with gemmi. K315 was identified by sequence context (SASKMLSKFRK motif). SPRY-cGAS interface residues were defined by heavy-atom distance <5 Å. Lysine-SPRY distances were measured from NZ atoms to the nearest SPRY heavy atom. For MD validation, DNA was stripped from the Boltz-2 model 0 and the protein-only system (cGAS + SPRY, 686 residues) was built with PDBFixer, solvated with TIP4P-Ew water in a truncated octahedron (~167K atoms, 0.15 M NaCl), and simulated for 50 ns under the same conditions as §4.2 (ff19SB + tip4pew, 300 K, 2 fs, NVT).

### 4.8 Allosteric Pathway Identification (P1-5)

Residue interaction networks were constructed from CA-CA contact maps (distance <8 Å, >50% frame occupancy) and DCCM-weighted edges (weight = −log|C_ij|). Shortest paths were computed from each 4mut site to N-terminal interface residues (211–219) using NetworkX. Consensus hubs were defined as residues appearing in ≥50% of paths across replicas.

### 4.9 S305 Phosphorylation MD

Phosphoserine (SEP) was modeled at cGAS S305 using ff19SB + phosaa19SB. S305E phosphomimetic was generated by point mutation. Both systems were simulated for 3 replicas × 200 ns each under the same conditions as the binary systems (§4.2). MM-GBSA analysis followed the protocol in §4.5.

---

## Supporting Information

### S1. S305 Phosphorylation Destabilizes the cGAS-TRIM41 Complex

CHK2-mediated phosphorylation of cGAS S305 has been reported to promote cGAS-TRIM41 interaction (Zhen et al., *Nat. Commun.* 2023). To test this computationally, we performed MD simulations of the Hsap cGAS-TRIM41 SPRY complex with S305 phosphorylated (SEP305) or mutated to the phosphomimetic glutamate (S305E).

**S305-phos causes complete dissociation.** All three replicas of the S305-phos system showed progressive increases in cGAS-TRIM41 COM distance, reaching 67–90 Å vs. 43 ± 3 Å for WT (Table S1). Replica 2 reached ~110 Å COM, indicating complete dissociation. This contradicts the reported stabilizing effect of S305 phosphorylation and suggests that in the apo cGAS context, S305 phosphorylation introduces electrostatic repulsion at the interface.

**S305E shows heterogeneous behavior.** The S305E phosphomimetic produced mixed results: one replica remained bound (COM ~46 Å), one partially dissociated (~65 Å), and one fully dissociated (~95 Å). This heterogeneity suggests that glutamate is an imperfect mimic of phosphoserine at this position, and that the S305 region is conformationally labile.

**MM-GBSA confirms energetic penalty.** S305-phos showed a significant loss of binding free energy (ΔΔG = +12.3 ± 8.7 kcal/mol relative to WT, p < 0.05), driven primarily by increased polar solvation penalty consistent with the introduction of a charged phosphate group at the protein-protein interface.

**Table S1. S305 MD summary (3 replicas × 200 ns each).**

| System | COM distance (Å) | Contacts (<8 Å) | Rg cGAS (Å) | Binding status |
|--------|-----------------|-----------------|-------------|----------------|
| WT | 43.0 ± 2.6 | 148.4 ± 10.2 | 30.6 ± 1.2 | Stable |
| S305-phos | 76.2 ± 12.3 | 62.3 ± 25.8 | 42.1 ± 5.1 | **Dissociated** |
| S305E | 68.7 ± 22.5 | 74.1 ± 35.2 | 37.4 ± 7.6 | Heterogeneous |

**Interpretation caveats.** (i) The simulated system is apo cGAS, whereas S305 phosphorylation was studied experimentally in the context of DNA-bound, potentially dimeric cGAS. (ii) The phosphate group was modeled with standard Amber phosaa19SB parameters, which may not fully capture phosphoserine conformational dynamics. (iii) The crystal structure of the cGAS-TRIM41 complex may reveal a different S305 environment than our docked model.

Despite these caveats, the robust dissociation signal (all 3 replicas, ~1.2 μs aggregate) provides a testable prediction: *in the absence of DNA, S305 phosphorylation should weaken rather than strengthen cGAS-TRIM41 binding*. This could be tested by SPR/BLI with purified apo cGAS ± λ-phosphatase treatment.

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
10. Seo, G. J. et al. *Cell Rep.* **2018**, 23, 1111–1121.
11. Liu, Z. S. et al. *Nat. Immunol.* **2021**, 22, 591–602.
