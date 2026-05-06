# Quaternary MVP Structure Assessment & Rebuild Plan

**Date**: 2026-04-23 (updated 2026-04-23)
**Status**: 🟡 Chai-1 predictions completed; WT vs 4mut interface difference is marginal; quaternary rebuild deprioritized.

---

## Executive Summary

Static and trajectory analyses of the E2~Ub-TRIM41-cGAS quaternary MVP reveal **serious structural problems** in both WT and 4mut systems. The current MD trajectories (WT ~5.7 ns, 4mut ~4.6 ns) sample incorrect conformational spaces and **have been stopped** pending structural correction.

However, subsequent efforts to rebuild the 4mut cGAS-SPRY interface have yielded important insights:

1. **Rosetta FastRelax from WT AF3 pose failed** — score explodes to 317,306 when coordinate constraints are released, indicating the 4mut clashes cannot be resolved by local optimization.
2. **Chai-1 truncated predictions (cGAS CTD + SPRY) succeeded on GPU** for both WT and 4mut.
3. **Chai-1 WT vs 4mut comparison shows marginal interface differences** — SPRY RMSD ~3–8 Å, comparable to within-set sampling noise.
4. **Mutation sites (D431S, K479E, L495Y, K498T) are distant from the cGAS-SPRY interface** in both WT and 4mut Chai-1 models (~24–39 Å).

**Key implication**: These four mutations do not directly disrupt the cGAS-SPRY binding interface. Any functional effect is likely allosteric or indirect.

---

## 1. Root Cause Analysis

### 1.1 4mut cGAS-SPRY Pose Failure

We had a script (`prepare_4mut_input.py`) designed to build 4mut input from the WT pose by rigid Kabsch alignment. However, this produced a structure with **catastrophic clashes**:

| Metric | WT Pose | 4mut Input (post-Kabsch) |
|--------|---------|--------------------------|
| cGAS-SPRY COM distance | 42.18 Å | 20.70 Å |
| Min CA distance | 5.07 Å | **0.36 Å** ❌ |
| CA contacts (<8 Å) | 18 | **999** ❌ |

**Why it failed**: Kabsch alignment is rigid-body. The 4mut mutations (D431S, K479E, L495Y, K498T) induce local conformational changes that cannot be accommodated by rigid translation/rotation. The aligned 4mut cGAS simply penetrates the SPRY surface.

This clash-ridden input was fed into Rosetta docking, which could not recover a valid pose:

| Metric | WT Docking | 4mut Docking |
|--------|-----------|--------------|
| CAPRI rank | **3** ✅ | **0** ❌ |
| Fnat | **0.857** | **0.000** |
| I_sc | -14.99 | -9.97 |
| Irms | 0.98 Å | 12.55 Å |

The final 4mut quaternary complex uses this failed docking pose, resulting in a **biologically meaningless cGAS-SPRY interface**.

### 1.2 Coordinate Drift Between Raw PDB and inpcrd

| System | Raw PDB K315→Ub | inpcrd K315→Ub | Δ |
|--------|-----------------|----------------|---|
| WT | 16.60 Å | 26.25 Å | +9.7 Å ❌ |
| 4mut (shifted) | 20.19 Å | 28.39 Å | +8.2 Å ❌ |

The tleap/pdb4amber pipeline somehow displaced the protein coordinates between raw PDB and inpcrd. The mechanism is unclear (possible re-boxing, centering, or an unintended translate command), but the consequence is that **MD starts from a wrong geometry**.

### 1.3 E2~Ub Conformational Drift

Even though the static E2~Ub RMSD vs 5FER is 0.00 Å (directly extracted from 5FER), during MD:
- WT: E2→Ub distance grows from 3.7 Å to **7.2 ± 1.5 Å**
- 4mut: E2→Ub stays ~4.8 Å (more stable)
- Both systems show **E2 RMSD vs 5FER of 29–57 Å**, indicating a global opening of the E2~Ub conformation

This suggests the RING-E2~Ub interface in the model may not be sufficiently stabilizing the closed state.

---

## 2. Rebuild Attempts & Results

### 2.1 Rosetta FastRelax from WT AF3 Pose ❌ FAILED

**Approach**: Apply 4mut mutations (D431S, K479E, L495Y, K498T) to the WT AF3 pose, then run FastRelax with coordinate constraints.

**Protocol**:
- `MutateResidue` movers for 4 mutations
- `FastRelax` with `constrain_relax_to_start_coords=1`, `coord_constrain_sidechains=1`, `ramp_constraints=1`
- `ref2015` score function
- 5 structures requested

**Results**:
- Initial score: 5,537.84 (very high — indicates major problems)
- During relaxation cycles, score drops to −81.89 (good!) with coord constraints active
- **When `coord_cst_weight` is ramped to 0**, score explodes to **99,263.5 → 172,995 → 317,306**
- Final RMSD: 56.67 Å

**Conclusion**: The 4mut mutations create clashes that **cannot be resolved by local relaxation**. When constraints are released, the structure collapses. Rigid-body mutation + local optimization is insufficient.

### 2.2 Chai-1 Truncated Prediction ✅ SUCCESS

**Approach**: Predict cGAS CTD (residues 200–522, 323 aa) + TRIM41 SPRY (413–630, 217 aa) = 540 aa total. Run on GPU 3 (RTX 3090, 24 GB) with `--low-memory`.

**Parameters**:
- `--num-trunk-recycles 3`
- `--num-diffn-timesteps 200`
- `--num-diffn-samples 5`
- `--device cuda:3`

**Status**: Both WT and 4mut completed successfully (~30–50 min each).

#### 2.2.1 Chai-1 WT Truncated Results

| Model | Confidence | cGAS-SPRY CA contacts (<8 Å) | Min CA distance | Interface cGAS | Interface SPRY |
|-------|-----------|------------------------------|-----------------|----------------|----------------|
| 0 | 0.6445 | 27 | 3.91 Å | 12 | 12 |
| 1 | 0.6515 | 36 | 3.94 Å | 19 | 15 |
| 2 | 0.6317 | 32 | 4.92 Å | 17 | 14 |
| 3 | 0.6368 | 24 | 5.70 Å | 16 | 11 |
| 4 | 0.6329 | 34 | 3.67 Å | 21 | 17 |

#### 2.2.2 Chai-1 4mut Truncated Results

| Model | Confidence | cGAS-SPRY CA contacts (<8 Å) | Min CA distance | Interface cGAS | Interface SPRY |
|-------|-----------|------------------------------|-----------------|----------------|----------------|
| 0 | 0.6310 | 21 | 5.54 Å | 19 | 14 |
| 1 | 0.6297 | 29 | 5.34 Å | 21 | 18 |
| 2 | 0.6522 | 31 | 5.85 Å | 25 | 24 |
| 3 | 0.6627 | 33 | 5.46 Å | 24 | 20 |
| 4 | 0.6576 | 30 | 3.84 Å | 26 | 18 |

#### 2.2.3 Chai-1 WT vs 4mut Direct Comparison

**Aligned on cGAS CTD, comparing SPRY position**:

| Comparison | SPRY RMSD | Interpretation |
|-----------|-----------|----------------|
| WT1 vs 4m3 (best models) | **4.27 Å** | Moderate difference |
| WT internal range | 3.26 – 9.90 Å | Sampling noise |
| 4mut internal range | 4.16 – 12.63 Å | Sampling noise |
| WT vs 4mut full range | 2.97 – 12.48 Å | Comparable to noise |

**Mutation site distances to nearest SPRY CA** (in Chai-1 models):

| Site | WT dist | 4mut dist | Δ |
|------|---------|-----------|---|
| D431 | 24.20 Å | 29.12 Å | +4.92 Å |
| K479 | 25.52 Å | 34.13 Å | +8.61 Å |
| L495 | 32.08 Å | 33.88 Å | +1.80 Å |
| K498 | 28.20 Å | 38.84 Å | +10.63 Å |

**Interface residue overlap** (best models, CA < 8 Å):
- cGAS interface: WT=19, 4mut=6, overlap=3
- SPRY interface: WT=15, 4mut=6, overlap=4

#### 2.2.4 Key Conclusions from Chai-1

1. **Interface difference between WT and 4mut is marginal** — SPRY RMSD (3–8 Å) is comparable to Chai-1's own sampling noise.
2. **Mutation sites are far from the interface** in both WT and 4mut (~24–39 Å). They do not make direct contact with SPRY.
3. **Chai-1 predicts a different cGAS-SPRY interface than AF3** — when compared to the AF3 WT pose, Chai-1 WT SPRY RMSD = ~88 Å. This is a method-level difference, not a mutation effect.
4. **Therefore, fair comparison requires same-method (Chai-1 WT vs Chai-1 4mut)**, which shows minimal differences.

---

## 3. Static Structure Analysis Results (Original)

### WT (quaternary_mvp_raw.pdb)

**Catalytic geometry**:
- K315 NZ → Ub G76 C: **16.60 Å**
- K315 NZ → E2 K85 NZ: 15.00 Å
- E2 K85 NZ → Ub G76 C: 2.41 Å
- K315 side chain angle toward Ub: **6.2°** ✅

**Interface contacts (CA-CA < 8 Å)**:
- cGAS–SPRY: 18 ✅
- cGAS–E2: 2
- cGAS–Ub: 0
- RING1–E2: 18 ✅
- E2–Ub: 57 ✅

**Overall**:
- Total Rg: 35.25 Å
- E2~Ub closed (RMSD vs 5FER): 0.00 Å

### 4mut (quaternary_mvp_4mut_shifted.pdb)

**Catalytic geometry**:
- K315 NZ → Ub G76 C: **20.19 Å**
- K315 NZ → E2 K85 NZ: 18.10 Å
- E2 K85 NZ → Ub G76 C: 2.41 Å
- K315 side chain angle toward Ub: **118.2°** ❌ (points away)

**Interface contacts**:
- cGAS–SPRY: **0** ❌
- cGAS–RING1: **315** ❌ (cGAS embedded in RING1)
- cGAS–E2: **247** ❌
- cGAS–Ub: **159** ❌
- RING1–E2: 18 ✅
- E2–Ub: 57 ✅

**Overall**:
- Total Rg: 35.44 Å
- E2~Ub closed (RMSD vs 5FER): 0.00 Å

---

## 4. Trajectory Analysis Results (~5.7 ns WT / ~4.6 ns 4mut)

| Metric | WT_rep1 | 4mut_rep1 | Interpretation |
|--------|---------|-----------|----------------|
| K315→Ub (mean ± std) | **29.28 ± 1.29 Å** | **25.56 ± 1.29 Å** | 4mut ~4 Å closer, but both >> catalytic window |
| K315→Ub (min–max) | 25.16 – 31.59 Å | 20.91 – 29.20 Å | No significant approach trend |
| E2→Ub (mean ± std) | **7.24 ± 1.49 Å** | **4.75 ± 0.72 Å** | 4mut isopeptide mimic more stable |
| Total Rg | 37.57 ± 0.32 Å | 42.66 ± 1.00 Å | 4mut overall more expanded |
| E2 RMSD vs 5FER | ~57 Å | ~29 Å | Both deviate massively from closed crystal structure |

**Conclusion**: Neither trajectory samples catalytically relevant geometry. Continuing to 50 ns would not fix the fundamental structural errors.

---

## 5. Implications for the 4mut Allosteric Hypothesis

### Prior hypothesis (from AF3 / Boltz-1 analysis)

The original working hypothesis was that D431S, K479E, L495Y, K498T (humanizing mutations on naked mole rat cGAS, or "reverse" mutations on human cGAS) alter cGAS function through **allosteric modulation of the cGAS-SPRY interface**. The AF3-predicted WT pose placed D431 and L495 within ~7 Å of SPRY, suggesting they could influence interface stability.

### Chai-1 evidence

1. **Chai-1 predicts D431/L495 are ~24–32 Å from SPRY** in both WT and 4mut. This contradicts the AF3 pose where they are ~7 Å from SPRY.
2. **The difference between Chai-1 WT and Chai-1 4mut is marginal** (SPRY RMSD ~3–8 Å, within sampling noise).
3. **Therefore, Chai-1 data does NOT support the hypothesis that these 4 mutations directly remodel the cGAS-SPRY interface**.

### Reconciliation

| Aspect | AF3 / Boltz-1 | Chai-1 | Interpretation |
|--------|--------------|--------|----------------|
| D431/L495 proximity to SPRY | ~7 Å (interface-adjacent) | ~24–32 Å (distant) | Methodological disagreement on interface location |
| WT vs 4mut interface difference | Large (different poses) | Marginal (3–8 Å) | Chai-1: mutations do not change binding mode |
| Interface confidence | AF3 iptm high for WT | Chai-1 conf ~0.63–0.66 | Moderate confidence |

**Chai-1 data contradicts the "interface-remodeling" component** of the allosteric hypothesis. However, it **does not rule out allosteric effects entirely** — the mutations could still propagate conformational signals through cGAS CTD without changing the cGAS-SPRY interface geometry.

**Revised hypothesis**: The 4 mutations may exert their effects through:
1. **Altering cGAS DNA-binding affinity** (mutations are in the CTD, which forms the DNA-binding groove)
2. **Modulating cGAS catalytic activity** (allosteric effects on the active site, not SPRY interface)
3. **Changing cGAS oligomerization state** (CTD mutations can affect dimerization)
4. **Indirect SPRY effects** through long-range allostery (possible but not supported by Chai-1 interface data)

---

## 6. Revised Action Items

| # | Task | Status | Priority | Notes |
|---|------|--------|----------|-------|
| 1 | 🛑 Stop current quaternary MDs | ✅ Done | P0 | ~5.7 ns WT / ~4.6 ns 4mut collected |
| 2 | 📝 Document findings | ✅ Done | P0 | This document |
| 3 | 🔧 Rosetta FastRelax 4mut from WT AF3 pose | ❌ **Failed** | P1 | Score explodes to 317k; route abandoned |
| 4 | 🤖 Chai-1 truncated prediction (WT + 4mut) | ✅ **Done** | P1 | GPU ~30–50 min; 5 structures each |
| 5 | 📊 Chai-1 WT vs 4mut comparison | ✅ **Done** | P1 | Interface difference marginal |
| 6 | 🔍 **Fix `analyze_s305e.py` final-50ns nan bug** | ⏳ Pending | P2 | Paper data needed |
| 7 | 🚀 **Continue Hgal MD to 200 ns** | 🟢 **Running** | P0 | 5 replicas active; ~60–104 ns current |
| 8 | 📝 **Revise allosteric hypothesis** | ⏳ Pending | P1 | Chai-1 contradicts interface-remodeling model |

---

## 7. Lessons Learned

1. **Rigid-body alignment is insufficient for mutant modeling**. Always check for clashes after alignment; if min distance < 2 Å, do relaxation or repacking.
2. **Always verify inpcrd vs raw PDB coordinates before MD**. A simple distance check on key atoms (K315→Ub) would have caught the coordinate drift immediately.
3. **CAPRI rank and Fnat are essential quality filters**. Fnat = 0 means the docking pose has zero native contacts—biologically meaningless.
4. **Trajectory monitoring should start from frame 0**, not just after several nanoseconds. Frame-0 distance checks validate the initial structure.
5. **Same-method comparison is required for fair mutant assessment**. Comparing AF3 WT vs Chai-1 4mut conflates method differences with mutation effects.
6. **Chai-1 predictions have significant sampling variability**. SPRY RMSD 3–12 Å within the same sequence means small cross-group differences (3–8 Å) may not be biologically meaningful.

---

## Appendix: Key File Paths

```
# Chai-1 truncated outputs
structures/quaternary_mvp/chai_wt_truncated_output/pred.model_idx_*.cif
structures/quaternary_mvp/chai_4mut_truncated_output/pred.model_idx_*.cif

# WT AF3 pose
structures/af3_raw/job1_Hsap_WT/ranked_0_chain_A.pdb   # cGAS (522 aa)
structures/af3_raw/job1_Hsap_WT/ranked_0_chain_B.pdb   # TRIM41 (630 aa)

# 4mut sequences
sequences/Hsap_cGAS_4mut.fasta                         # D431S,K479E,L495Y,K498T
sequences/cgas_trim41_sequences.fasta

# Failed 4mut docking
structures/docking/rosetta/hsap_4mut_global.sc         # CAPRI rank 0

# Analysis scripts
scripts/03_analysis/analyze_quaternary_static.py
scripts/03_analysis/analyze_quaternary_trajectory.py

# Stopped MD outputs
data/md_runs/quaternary_mvp/WT_rep1/analysis.csv
data/md_runs/quaternary_mvp/4mut_rep1/analysis.csv
```
