# cGAS-TRIM41 结合分析：计算实验设计方案

> 基于论文：*A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging* (Chen et al., Science 2025)
>
> 目标：通过计算模拟，解释/验证裸鼹鼠 cGAS 4个氨基酸变异（S463, E511, Y527, T530）如何影响其与 TRIM41 的功能性相互作用，进而影响泛素化与 DNA 修复。

---

## 一、PDB 结构现状

| 蛋白 | 实验结构 | 备注 |
|------|----------|------|
| 人 cGAS | ✅ 有多个（如 `4LEZ`, `6CT9`） | 包括 DNA 结合态、二聚体 |
| 裸鼹鼠 cGAS | ❌ 无 | 需同源建模或 AF 预测 |
| TRIM41 | ❌ 无实验结构 | UniProt: `Q8WV44` |
| cGAS-TRIM41 复合物 | ❌ 无 | 2025 年 Science 才报道 |

**结论**：所有工作必须从预测结构开始。

---

## 二、方案分层

### 🔹 方案 A：快速筛选（1-2 周）

**目标**：快速判断 4 个突变是否落在 TRIM41 结合界面上，低成本。

#### 步骤
1. **结构预测**：
   - 工具：**ColabFold**（AlphaFold-Multimer v2/v3）或 **AlphaFold3**
   - 预测以下复合物：
     - 人 cGAS (WT) + TRIM41
     - 人 cGAS_4mut + TRIM41（人→裸鼹鼠）
     - 裸鼹鼠 cGAS (WT) + TRIM41
     - 裸鼹鼠 cGAS_4mut_rev + TRIM41（裸鼹鼠→人）

2. **界面质量初筛**：
   - 比较 **ipTM**（interface pTM）和 **pLDDT**
   - ipTM > 0.7 表示高置信度界面；< 0.5 可能不稳定
   - 检查 4 个突变位点是否直接接触 TRIM41

3. **快速能量估算**：
   - **FoldX** 或 **Rosetta ddG**：单点/组合突变扫描
   - 计算 ΔΔG_bind，若 > 1 kcal/mol 则认为显著影响结合

#### 预期输出
- 哪些突变位点直接参与 TRIM41 结合
- WT vs 突变体的结合强度趋势（定性/半定量）

---

### 🔹 方案 B：中等精度验证（1-2 个月）【推荐主力】

**目标**：半定量比较 WT 与突变体的结合自由能，支持论文结论。

#### 步骤
1. **结构预测**（同方案 A，用 AF-Multimer/AF3）

2. **MD 模拟**（关键改进）：
   - 力场：`Amber ff19SB` + `OPC` 水 或 `CHARMM36m`
   - 体系：cGAS-TRIM41 复合物 + 显式水 + 150 mM NaCl
   - 时长：**至少 500 ns - 1 μs**（蛋白-蛋白界面收敛慢）
   - ⚠️ **不要用短于 100 ns 的轨迹算结合能**

3. **结合自由能计算**：
   - **MM-PBSA**：`gmx_MMPBSA` 或 `Amber MMPBSA.py`
   - 更好的选择：**MM-GBSA**（对界面极性残基更敏感，计算更快）
   - 做 **能量分解**（per-residue），看 4 个突变位点各自贡献
   - ΔΔG_bind = G_bind(突变体) - G_bind(WT)

4. **增强采样**（可选但推荐）：
   - **aMD** / **GaMD**：观察界面开合、解离趋势
   - 或 **umbrella sampling**：计算解离路径的自由能面

#### 风险控制
- 若 MD 后 RMSD > 3 Å 且界面明显松散 → **预测结构本身不可靠，停止**
- 需做至少 **3 个独立重复**，统计显著性

---

### 🔹 方案 C：高精度计算（3-6 个月，仅发高水平计算论文时）

**目标**：定量计算精确 ΔΔG。

- **FEP（自由能微扰）**：`Amber TI` 或 `Schrodinger FEP+`
- 做 alchemical mutation（人→4mut 逐步突变）
- 需要多个 λ 窗口（5-10 个），每个 10-20 ns
- 误差可控制在 ~0.5 kcal/mol

⚠️ **风险**：TRIM41 无实验结构，FEP 建立在预测结构上的可靠性会打折扣。

---

## 三、推荐的具体流程

```
1. 序列分析
   └─> 4个位点是否在TRIM41已知底物识别motif附近？

2. 结构预测（AF-Multimer/AF3）
   ├─> 人cGAS + TRIM41
   ├─> 人cGAS_4mut + TRIM41
   ├─> 裸鼹鼠cGAS + TRIM41
   └─> 裸鼹鼠cGAS_4mut_rev + TRIM41

3. 界面分析（不用MD）
   ├─> 比较 ipTM / pLDDT
   ├─> 界面面积、氢键数、盐桥数
   └─> 4个位点是否直接接触TRIM41？

4. 短MD弛豫（100-200 ns）
   └─> 观察 RMSF：突变位点柔性是否增加？
       └─> 若4mut导致界面loop更flexible → 结合可能不稳定

5. 结合能计算（MM-PBSA/GBSA）
   └─> 至少3个独立重复，统计显著性

6. 对照实验
   ├─> 阴性对照：已知无相互作用的蛋白对
   └─> 阳性对照：文献中有Kd数据的TRIM41-底物
```

---

## 四、务实替代方案

若 AF-Multimer 预测的复合物置信度低（ipTM < 0.6）：

1. **只做单体 MD**：比较人 vs 裸鼹鼠 cGAS 表面构象差异，看 4 个突变是否改变了泛素化位点附近的局部结构。
2. **分子表面分析**：预测 TRIM41 底物识别表面，比较 cGAS 表面静电势/形状互补性。
3. **等实验结构**：建议合作方尝试冷冻电镜。

---

## 五、关键参数速查

| 项目 | 推荐 |
|------|------|
| 结构预测 | ColabFold (AF-Multimer v3) / AlphaFold3 |
| MD 力场 | Amber ff19SB + OPC / CHARMM36m |
| MD 时长 | ≥ 500 ns（结合能计算） |
| 结合能方法 | MM-GBSA（快）+ MM-PBSA（交叉验证） |
| 快速突变扫描 | FoldX / Rosetta ddG |
| 高精度 ΔΔG | FEP+（仅结构可靠时） |

---

## 六、项目当前进度（2026-04-30）

### 已完成

| 阶段 | 任务 | 结果 |
|------|------|------|
| 结构预测 | AlphaFold3 (4 systems) | ✅ 完成，详见 `docs/af3_report.md` |
| 分子对接 | ClusPro (Hsap_WT/Hsap_4mut/Hgal_WT/Hgal_4mut_rev) | ✅ 完成，详见 `docs/docking_report.md` |
| MD 准备 | OpenMM 体系构建 (prmtop/rst7) | ✅ 4 物种 × 3 reps |
| MD 生产 | OpenMM Hsap_WT 200ns × 3 reps | ✅ 完成 |
| MD 生产 | OpenMM Hgal_WT 200ns × 3 reps | ✅ 完成 |
| 磷酸化 | S305-phos 体系构建 + EM | ✅ 完成 |
| 磷酸化 | OpenMM S305-phos 200ns × 3 reps | 🔄 ~15.5ns/200ns，运行中 |
| GROMACS 验证 | 旧转换（CMAP bug）Hsap_WT/Hsap_4mut | ✅ 完成，但数据不可靠 |
| GROMACS 验证 | GROMACS 2026 native amber19sb.ff | 🔄 ~9.8ns/200ns，运行中 |

### 关键发现

1. **GROMACS 与 OpenMM 差异**: 旧 GROMACS（parmed 转换）RMSD 是 OpenMM 的 3-4×，根本原因是 CMAP 残基特异性丢失（14 types → 1 type）。已用 GROMACS 2026 原生 `amber19sb.ff` 修复，验证运行中。
2. **磷酸化位点**: S305（CHK2 靶点）是当前构建体（200-522）内唯一可模拟的磷酸化位点。
3. **S305-phos 构建**: 需先能量最小化（-215k → -1033k kJ/mol）再运行 MD，否则 heating 阶段 NaN。

### 下一步

| 优先级 | 任务 |
|--------|------|
| 🔴 高 | 等待 S305-phos 3× replica 完成（~12h） |
| 🔴 高 | 等待 GROMACS 2026 验证完成（~32h），对比 RMSD |
| 🟡 中 | 构建 S305E 体系 + 3× replica MD |
| 🟡 中 | 分析已完成数据：Hsap_WT vs Hgal_WT vs 4mut |
| 🟢 低 | 磷酸化 vs 磷酸模拟（S305-phos vs S305E）对比分析 |

---

*文档创建时间：2026-04-22*
*更新时间：2026-04-30*
