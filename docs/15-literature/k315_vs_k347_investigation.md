# cGAS 泛素化位点调研：K315 vs K347 结论与后续验证建议

> 调研日期：2026-05-11  
> 调研人：Kimi Code CLI  
> 背景：项目计算预测 K315 为 TRIM41 介导泛素化的唯一几何可达位点，与文献中提及的 K347 存在张力。本次调研旨在厘清文献事实、评估结论可靠性并提出验证路径。

---

## 一、核心结论

### 1.1 文献事实：不存在"TRIM41 泛素化 K347"的实验证据

| 论文 | E3 连接酶 | 泛素化位点 | 证据强度 | 备注 |
|------|----------|-----------|---------|------|
| **Seo et al. 2018** *Nat. Commun.* 9:613 [1] | **TRIM56** | 小鼠 K335 = **人源 K347** | ✅ 质谱（Lys278/335/350）+ K→R 突变 + 功能实验 | **位点已确证** |
| **Liu et al. 2018** *Cell Biosci.* 8:35 [2] | **TRIM41 (RINCK)** | **未指定位点** | ⚠️ 仅证明 TRIM41 能泛素化 cGAS 并促进其活性 | **无位点鉴定** |

- **K347 是 TRIM56 的实验确证位点**，功能为促进 cGAS 二聚化和 DNA 结合。
- **TRIM41 的泛素化位点从未被任何实验鉴定**。文献仅笼统表述为"mediates monoubiquitination of cGAS"，无质谱位点 mapping，无 K→R 扫描。
- Chen et al. (2025) *Science* 也未提及任何具体泛素化位点，仅描述为"weakening TRIM41-mediated ubiquitination"。

### 1.2 项目内部文献综述中的错误

| 文档 | 错误内容 | 正确事实 |
|------|---------|---------|
| `docs/15-literature/人cGAS泛素化位点.md` §4.1 及表格 | 将 K347 同时归给 TRIM56 和 TRIM41 | K347 仅由 **TRIM56** 介导；**TRIM41 位点未鉴定** |
| `paper/paper_final.md` §3.3 | "The TRIM41 literature identifies K347—not K315—as the primary monoubiquitination site [10,11]" | 文献 [10] 说的是 **TRIM56**→K347；文献 [11] 说的是 TRIM41→cGAS 但**未指定位点**。不存在"TRIM41 literature identifies K347"。 |
| `paper/paper_final.md` References | [10] 期刊信息错误（写成 *Cell Rep.* 23:1111）；[11] 期刊信息错误（写成 *Nat. Immunol.* 2021） | [10] 应为 *Nat. Commun.* 2018, 9, 613；[11] 应为 *Cell Biosci.* 2018, 8, 35 |

### 1.3 K315 结论的可靠性评估

**靠谱的部分：**
1. 在 **apo cGAS-TRIM41 四元复合物模型**中，K315 是唯一几何可达的赖氨酸（全 36 Lys 扫描，K315 < 25 Å，其余 > 30–70 Å）。
2. 4mut 通过长程变构效应将 K315 的 PMF 最小值偏移 −2.5 Å、催化准备度提升 4.3×——这一**机制不依赖于具体位点是 K315 还是 K347**。
3. TRIM41 的具体泛素化位点确实是文献空白，K315 作为计算预测的新假说具有科学价值。

**需要谨慎的部分：**
1. **模型基于 apo 单体**，而 Chen 2025 的实验 readout（抗 K48-Ub）发生在 **DNA-bound、可能二聚化的核内 context**。DNA-bound 时 SPRY 结合 C 端面（Boltz-2 预测），K315 反而不可达（36–57 Å）。
2. **K347 位于 cGAS-DNA 二聚化界面附近**，在 DNA-bound/二聚体状态下可能通过 **trans-泛素化**（一个 cGAS 分子被另一个分子上的 TRIM41 修饰）变得可达。
3. K315 作为泛素化位点**无任何实验文献支持**，在论文中必须明确标注为"计算预测"。

---

## 二、建议的进一步确认（按优先级）

### 🔴 高优先级：文献与文本修正

1. **修正 `paper/paper_final.md` §3.3 的表述**
   - 明确区分：TRIM56→K347（实验确证，二聚化相关）vs TRIM41→位点未鉴定。
   - 将 K315 定位为"我们的计算预测"，而非对文献的挑战。
   - 修正参考文献 [10][11] 的期刊信息。

2. **修正 `docs/15-literature/人cGAS泛素化位点.md`**
   - 表格中 K347 的 E3 列只保留 TRIM56，移除 TRIM41。
   - §4.1 中明确说明 TRIM41 的位点尚未被实验鉴定。

### 🟡 中优先级：计算层面验证

| 验证实验 | 目的 | 方法/资源 | 预期结果 |
|---------|------|----------|---------|
| **cGAS-DNA 二聚体 + TRIM41 建模** | 测试 K347 在二聚体/核内状态下是否可达 | Boltz-2 或 AF3（支持 DNA 输入），构建 cGAS-DNA 二聚体 + TRIM41 SPRY 复合物 | 如果 DNA-bound 界面在 C 端面，K315 不可达；K347 可能因二聚体界面暴露或 trans-泛素化而可达 |
| **TRIM56 vs TRIM41 RING 域结构比对** | 两者 E2~Ub 催化模块的取向是否不同，导致靶点偏好差异 | AlphaFold / PDB 结构比对 | 如果 RING 域结构/表面电荷差异显著，支持不同 E3 偏好不同位点的假说 |
| **K315R/K347R 双突变的 Rosetta ddG 扫描** | 测试两个位点各自对 TRIM41 泛素化能量学的贡献 | Rosetta FastRelax + ddG | K315R 若在能量上显著不利于催化几何，支持其为主要位点；K347R 若无显著变化，支持其不是 TRIM41 靶点 |
| **扩大 CC linker 长度的 US/PMF** | 当前 CC40 限制在 40–100 Å，可能人为排除了 K347 可达的构象 | 测试 CC60/CC80，观察 K347 距离分布是否出现新的低能量极小值 | 如果即使 CC80 K347 仍 > 40 Å，则进一步支持 K347 在 apo 状态下不可达 |
| **2D PMF（距离 × 二面角）** | 当前 1D PMF 沿 K315→Ub G76 距离，可能遗漏角度自由度 | 计算 K315 NZ→Ub G76 C 距离 + K315 NZ→E2 K85 NZ→Ub G76 C 二面角 | 验证 4mut 是否在更完整的反应坐标上仍偏向催化有利区域 |

### 🟢 低优先级：湿实验方向（供合作者参考）

| 实验 | 设计 | 预期结果 |
|------|------|---------|
| **体外泛素化（TRIM41 + cGAS WT/K315R/K347R）** | 纯化 TRIM41 RING+SPRY + E2(UbcH5b) + Ub + cGAS(WT/K315R/K347R) | 若 K315R 显著降低 TRIM41 介导的泛素化，而 K347R 无影响 → 支持 K315 为 TRIM41 位点 |
| **质谱位点 mapping** | TRIM41 体外泛素化后胰酶消化，LC-MS/MS 检测 K-ɛ-GG 残基 | 直接鉴定 TRIM41 的泛素化位点，终结假说争议 |
| **DNA 依赖性界面转换验证** | XL-MS 检测 cGAS-TRIM41 ± DNA | 验证 apo 时 N 端界面 vs DNA-bound 时 C 端界面的切换 |

---

## 三、对论文叙事的影响

当前论文的核心机制——"4mut 通过长程变构调控催化几何"——**不依赖于 K315 是否是真正的生理位点**。即使最终实验发现 TRIM41 的靶点是另一个赖氨酸，以下结论仍然成立：

1. 4mut 位点不在 TRIM41 物理界面上（30–39 Å）。
2. 4mut 诱导 N 端长程变构位移（12.3 Å）。
3. 4mut 不改变结合亲和力（MM-GBSA p = 0.50）。
4. 4mut 重塑构象动态（DCCM/PCA）。
5. 在 apo 模型中，K315 是唯一几何可达靶点，且 4mut 将其 PMF 偏移 −2.5 Å。

**建议叙事调整：**

> "TRIM56 has experimentally established K347 (murine K335) as the monoubiquitination site promoting cGAS dimerization on DNA [10]. However, **no study has mapped the TRIM41 ubiquitination site on cGAS** [11]. In our apo-state quaternary model, K315 is the only lysine geometrically accessible to the TRIM41 RING-E2~Ub catalytic center, while K347 remains >50 Å away. The catalytic geometry of K315 is significantly biased toward the catalytic center by 4mut (PMF shift −2.5 Å; readiness 4.3×). Whether K315 represents the physiological TRIM41 target in the apo state, or whether DNA binding exposes alternative lysines (e.g., K347 via trans-ubiquitination), remains to be tested experimentally."

---

## 四、References

[1] Seo, G. J.; Kim, C.; Shin, W. J.; Sklan, E. H.; Eoh, H.; Jung, J. U. TRIM56-Mediated Monoubiquitination of cGAS for Cytosolic DNA Sensing. *Nat. Commun.* **2018**, *9*, 613. DOI: 10.1038/s41467-018-02936-3

[2] Liu, Z. S.; Zhang, Z. Y.; Cai, H.; Zhao, M.; Mao, J.; Dai, J.; et al. RINCK-Mediated Monoubiquitination of cGAS Promotes Antiviral Innate Immune Responses. *Cell Biosci.* **2018**, *8*, 35. DOI: 10.1186/s13578-018-0231-4

[3] Chen, Y. et al. A cGAS-Mediated Mechanism in Naked Mole-Rats Potentiates DNA Repair and Delays Aging. *Science* **2025**, *390*, eadp5056.
