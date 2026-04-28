# 评审意见验证报告

> 日期：2026-04-27
> 验证人：Kimi Code CLI（当前会话）
> 验证方式：静态代码审查 + 实际系统数据验证

---

## 验证结论总览

| 评审来源 | 问题 | 评级 | 验证结果 | 影响范围 |
|---------|------|------|---------|---------|
| Gemini #1 | Docking Pose Bias | 🔴 致命 | ✅ 属实 | Hgal WT vs 4mut_rev 比较无效 |
| Gemini #2 | US 硬编码 RING 1-43 | 🔴 致命 | ⚠️ 部分属实 | `run_umbrella_sampling.py` 有；当前用 `run_us_simple.py` 也有 |
| Gemini #3 | 变构效应证据链断裂 | 🟡 高 | ✅ 属实 | 缺少 DCCM/PCA |
| Gemini #4 | 盐浓度注释不符 | 🟡 中 | ✅ 属实 | 所有系统仅为 neutralize |
| Gemini #5 | PDB split 逻辑脆弱 | 🟢 低 | ⚠️ 属实 | 对当前数据有效，但非通用 |
| Gemini #6 | Hsap 采样量不足 | 🟡 高 | ✅ 属实 | 单 replica 无统计力 |
| Kimi #1 | t-test 误用 | 🔴 致命 | ✅ 属实 | 所有跨系统 p 值不可信 |
| Kimi #2 | 突变位点未清理 | 🔴 致命 | ✅ 属实 | 旧 fasta + 旧脚本残留 |
| Kimi #3 | 文档 0.78Å 矛盾 | 🟡 中 | ✅ 属实 | docking_report.md 误导 |
| Kimi #4 | 方法学混杂 | 🟡 高 | 同 Gemini #1 | — |
| Kimi #5 | 盐浓度注释 | 🟡 中 | 同 Gemini #4 | — |
| Kimi #6 | 硬编码路径 | 🟢 低 | ✅ 属实 | rosetta_mutational_scan.py |
| Kimi #7 | 文档重复 | 🟢 低 | ✅ 属实 | project_log + interface_analysis |
| Kimi #8 | MMPBSA 参数 | 🟡 中 | ⚠️ 待确认 | `-yr/-yl` vs `-mr/-ml` |
| Kimi #9 | US 力常数过大 | 🟡 中 | ⚠️ 可能属实 | k=1000 偏大 |
| Kimi #10 | 截断局限 | 🟢 低 | ✅ 属实 | 需讨论 |

---

## 新发现的致命问题（评审未提及）

### N1. "Lys-334" 是根本不存在的位点 — 所有 US 标注完全错误

**验证过程**：
1. `run_us_simple.py` 第23行：`res.index + 1 == 334 and atom.name == "NZ"`
2. 用 OpenMM 加载 `Hsap_WT.prmtop`，遍历 residues：
   - 索引 333（1-based=334）→ `res.name = LYS`，有 NZ 原子
3. 但 PDB 中链 B resid 334 → `THR`（苏氨酸）
4. 序列验证：`Hsap_cGAS_WT.fasta` 第334位 = `T`（不是 K）

**根因**：Amber prmtop 把蛋白质重新编号为连续残基 1-N。拓扑 resid 334 ≠ 序列残基 334。

**映射关系**（Hsap 系统）：
- TRIM41 SPRY (413-630) = 218 aa → 拓扑 resid 1–218
- cGAS CTD (200-522) = 323 aa → 拓扑 resid 219–541
- **拓扑 resid 334 = cGAS-315 (LYS)** ✅
- **序列残基 334 = cGAS-334 (THR)** ❌ 没有 NZ 原子

**影响**：
- `run_us_simple.py` 约束的是 **TRIM41-413~455 ↔ cGAS-315**
- 但所有代码注释、文档、图表、PMF 标题都写的是 **"Lys-334"**
- `analyze_lys_ubiquitination.py` 报告 resid 334 是 top candidate → 实际是 cGAS-315
- `analyze_close_state.py` 分析的是 resid 334 → 实际是 cGAS-315
- `project_log.md` 所有 "Lys-334" 结论都基于错误编号

**数据有效性**：US 数据本身科学有效（约束的是真实的 top LYS），但标注错误会导致论文被拒。

**修复策略**：
- 方案 A（推荐）：统一将所有 "Lys-334" 改为 "Lys-315"，保留现有 US 数据
- 方案 B（不可行）：cGAS-334 是 THR，不可能作为泛素化位点

---

## 其他关键验证细节

### Kimi #1: t-test 误用 — 100% 属实

`compare_systems.py` 第33、65行和 `compare_four_systems.py` 第33、76、109行直接对整个 MD 时间序列（每 0.5 ps 一帧）做 Welch t-test。

MD 帧的自相关时间约 1-10 ps，200ns 轨迹有 ~40,000 帧，有效独立样本数可能只有 ~200-2000。直接 t-test 把标准误低估了 ~10 倍，p 值被严重压缩。

**修复**：计算自相关时间 → block averaging（block size > 自相关时间）→ 重新统计检验。

### Kimi #2: 突变位点未清理 — 100% 属实

- `sequences/Hsap_cGAS_4mut.fasta` 仍是旧错误 `C463S,K479E,L495Y,K498T`
- `sequences/Hsap_cGAS_4mut_corrected.fasta` 才是正确的 `D431S,K479E,L495Y,K498T`
- `calc_residue_distances.py` 第74行：`hsap_residues = {463: 'S', 479: 'E', 495: 'Y', 498: 'T'}` — 未更新
- `rosetta_mutational_scan.py` 第220行：`hgal_mutations = [(463, 'D'), ...]` — Hgal 也用了旧映射

### Gemini/Kimi #4: Docking Pose Bias — 100% 属实

- Hgal_WT：LightDock best_pose
- Hgal_4mut_rev：Rosetta docking_protocol
- Hsap_WT：Rosetta
- Hsap_4mut：Rosetta

Hgal WT 与 4mut_rev 的界面差异（N-端切换）可能主要由 docking 方法差异驱动，而非突变效应。

### Kimi #3: 0.78Å 矛盾 — 100% 属实

`docking_report.md` 第274行："RMSD = 0.78 Å（3592 atoms），整体折叠高度相似"

但 3592 atoms ≈ 450 CA atoms（两个蛋白 CTD），0.78Å 只可能是某个核心局部 domain 的对齐结果。全长/CTD 全局 CA RMSD 实际约 20Å（`interface_analysis_report.md`）。

---

## 待确认项

### MMPBSA.py 参数

`run_mmpbsa.py` 使用 `-yr` 和 `-yl`：
```python
"-yr", receptor_mask,
"-yl", ligand_mask,
```

AmberTools MMPBSA.py 的标准参数：
- `-mr` / `-ml`：receptor/ligand mask（用于从复合物拓扑中分割）
- `-yr` / `-yl`：receptor/ligand trajectory（用于单独轨迹）

当受体和配体是同一轨迹时，应使用 `-mr` 和 `-ml`。需要实际运行验证。

### US 力常数

当前 k=1000 kJ/mol/nm²。蛋白质 US 常用 100-500。k=1000 可能导致：
- 窗口采样过度约束
- CV 分布过窄
- PMF 峰谷被抹平

建议用 k=200-500 测试 1-2 个窗口对比。
