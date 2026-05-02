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
   - **MM-GBSA**：使用 `MMPBSA.py` (AmberTools)，`igb=5` (GB-OBC II)
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

2. 结构预测（AF3 + Boltz-2 交叉验证）
   ├─> 人cGAS + TRIM41（全长 + 截断）
   ├─> 人cGAS_4mut + TRIM41
   ├─> 裸鼹鼠cGAS + TRIM41
   ├─> 裸鼹鼠cGAS_4mut_rev + TRIM41
   └─> ⚠️ **当 ipTM < 0.6 时，不直接信任任一模型的复合物构象**（见下方教训）

3. 界面分析（不用MD）
   ├─> 比较 ipTM / pLDDT（AF3 vs Boltz-2）
   ├─> 界面面积、氢键数、盐桥数
   ├─> 4个位点是否直接接触TRIM41？
   └─> ⚠️ **若 AF3 与 Boltz-2 的界面接触网络完全不同 → 预测不可靠，需对接+MD**

4. 蛋白对接（当预测置信度低时）
   ├─> 工具：ClusPro（remote）或 SDOCK2.0（local）
   ├─> 输入：单体预测结构（可信部分）
   └─> 输出：1000 个对接构象 → 聚类 → 选择 top poses

5. 短MD弛豫 + 验证（100-200 ns）
   ├─> 观察 COM distance 和 Rg：复合物是否稳定？
   ├─> 对比不同初始构象的收敛行为
   └─> 若4mut导致界面loop更flexible → 结合可能不稳定

6. 结合能计算（MM-GBSA）
   └─> 至少3个独立重复，统计显著性

7. 对照实验
   ├─> 阴性对照：已知无相互作用的蛋白对
   └─> 阳性对照：文献中有Kd数据的TRIM41-底物
```

---

## 四、务实替代方案

若 AF3/Boltz-2 预测的复合物置信度低（ipTM < 0.6）：

1. **单体 MD + 对接 + 平衡 MD**（本项目的实际路径）：
   - 使用单体预测结构（高 pLDDT 区域可信）
   - 蛋白对接获得初始复合物构象
   - 延长 MD（200-500 ns）观察自然收敛态
   
2. **模型间一致性检验**：
   - 对比 AF3 和 Boltz-2 的界面几何
   - 若两者单体一致但界面不同（如 cGAS-TRIM41 案例），说明界面不可信
   - 物理合理性检查： clashes、buried surface area、氢键网络
   
3. **实验验证优先**：
   - 建议合作方尝试冷冻电镜
   - 或 Co-IP / 荧光共定位实验确认界面存在性

---

## 五、关键参数速查

| 项目 | 推荐 | 备注 |
|------|------|------|
| 结构预测 | AlphaFold3 + Boltz-2 交叉验证 | 低 ipTM 时界面不可靠，需交叉验证 |
| MD 力场 | Amber ff19SB + OPC / CHARMM36m | |
| MD 时长 | ≥ 500 ns（结合能计算） |
| 结合能方法 | MM-GBSA（`igb=5`, GB-OBC II） |
| 快速突变扫描 | FoldX / Rosetta ddG |
| 高精度 ΔΔG | FEP+（仅结构可靠时） |

---

## 六、项目当前进度（2026-05-01）

### 已完成

| 阶段 | 任务 | 结果 |
|------|------|------|
| 结构预测 | AlphaFold3 (4 systems) | ✅ 完成，详见 `docs/10-reports/af3_report.md` |
| 分子对接 | ClusPro (4 systems) | ✅ 完成，详见 `docs/10-reports/docking_report.md` |
| MD 生产 | OpenMM Hsap_WT 200ns × 3 reps | ✅ 完成 |
| MD 生产 | OpenMM Hgal_WT 200ns × 3 reps | ✅ 完成 |
| GROMACS 验证 | 旧转换（CMAP bug）Hsap_WT/Hsap_4mut | ✅ 完成，数据不可靠，不纳入分析 |
| 磷酸化 | S305-phos 3× replica MD | ✅ 200ns/200ns，完成 |
| 磷酸化 | S305E 体系构建 + 3× replica MD | 🔄 0ns/200ns，运行中 |
| GROMACS 验证 | GROMACS 2026 native amber19sb.ff | 🔄 ~83ns/200ns，运行中 |
| Boltz-2 验证 | 全长 + 截断 + 与 AF3 结构对比 | ✅ 完成 |

### 关键发现

1. **GROMACS CMAP 修复成功**: GROMACS 2026 原生 `amber19sb.ff` 的 COM/Rg 与 OpenMM 几乎完全相同。RMSD 差异从 4× 降至 1.3×。重要教训：分析 GROMACS 轨迹必须先修复 PBC 包裹。
2. **S305-phos 导致解离**: 3 个 replica 的 COM 距离（67-90 Å）均显著大于 WT（45 Å），rep2 完全解离（COM~110 Å）。与文献"促进结合"论断相反。可能解释：溶液环境差异、力场对 -2 电荷排斥估计过度、或解离为中间态。
3. **Boltz-2 vs AF3 界面分歧**: 对截断 cGAS-TRIM41，AF3 和 Boltz-2 的 TRIM41 相对取向完全不同（RMSD=21.34 Å, Jaccard=0.00）。cGAS 单体一致（RMSD=1.06 Å），但界面预测不可信。进一步验证了对接+MD 策略的必要性。

### 下一步

| 优先级 | 任务 | ETA |
|--------|------|-----|
| 🔴 高 | 等待 S305-phos 3× replica 完成 | ~6.6h |
| 🔴 高 | 等待 GROMACS 2026 验证完成 | ~20h |
| 🟡 中 | S305E 电荷对照体系构建 | 视 S305-phos 最终结果决定 |
| 🟡 中 | 200ns 后 MM-GBSA 能量分解 | 等 MD 完成 |
| 🟡 中 | 分析 Hsap_WT vs Hgal_WT vs 4mut 已完成数据 | 可随时开始 |
| 🟢 低 | 全长体系（S120-phos）是否需要 | 待决策 |
| 🟢 低 | 论文图表制作 | 等分析完成 |

---

*文档创建时间：2026-04-22*
*更新时间：2026-04-23*
