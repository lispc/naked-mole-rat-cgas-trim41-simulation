# cGAS-TRIM41 项目审查报告

> 审查人：Kimi Code CLI  
> 审查日期：2026-04-27  
> 审查范围：全部核心文档（README.md、project_log.md、af3_mutation_analysis.md、interface_analysis_report.md 等）及主要代码（build_system.py/v2、run_md.py、analyze_system.py、compare_systems.py、compare_four_systems.py、run_umbrella_sampling.py、rosetta_mutational_scan.py 等）  
> 审查方式：静态代码与文档审阅，未重新运行程序

---

## 🔴 严重问题（直接影响科学结论可信度）

### 1. 统计方法根本错误：对自相关时间序列误用独立样本 t-test

**位置**：`scripts/compare_systems.py`（第33、65行）、`scripts/compare_four_systems.py`（第33、76、109行）

**问题**：对 MD 轨迹的 RMSD 和 COM 距离时间序列直接使用 Welch's t-test。MD 轨迹是强自相关的时间序列（相邻帧高度相关），完全不满足独立同分布（i.i.d.）假设。对自相关数据使用 t-test 会严重低估标准误，导致虚假的极端 p 值（如文档中大量出现的 p < 1e-300、p ≈ 0）。

**后果**：README 和 project_log 中所有的"统计显著"结论（如"Hsap 4mut 使复合体更不稳定，p < 1e-30"）都是**不可信的**。如果直接写入论文，会被审稿人直接质疑甚至拒稿。

**修复建议**：
- 计算自相关时间（autocorrelation time）
- 对时间序列进行 block bootstrap 或子采样（subsampling，间隔大于自相关时间）后再做统计检验
- 重新计算所有跨系统比较的 p 值

---

### 2. 突变位点编号严重错误，且清理不彻底

**位置**：`sequences/Hsap_cGAS_4mut.fasta`、`scripts/calc_residue_distances.py`、`docs/af3_mutation_analysis.md`

**问题**：项目最初将人类 4mut 错误定义为 C463S/K479E/L495Y/K498T，但论文原文实为 **D431S/K479E/L495Y/K498T**。虽然后来发现了错误并生成了 `*_corrected.fasta`，但：
- 旧错误文件 `Hsap_cGAS_4mut.fasta` 仍然保留在项目中，未删除或重命名
- `calc_residue_distances.py` 第74行仍然使用旧 residues（Hsap: 463/479/495/498），没有更新为 corrected 版本
- `docs/af3_mutation_analysis.md` 中的大量数据（如 Hgal_rev 的 RMSD = 13.043Å）是基于**错误的反向突变**（S463C 而非正确的 S463D）

**后果**：如果误用旧文件或旧脚本，会重新产生错误数据；文档中存在基于错误突变的大量历史数据，容易造成混淆。

**修复建议**：
- 删除或归档旧错误 fasta，将 corrected 文件设为标准名称
- 更新所有脚本中的 residue 编号
- 在文档中明确标注哪些分析基于错误位点、哪些基于正确位点

---

### 3. 文档数据严重自相矛盾

**位置**：`docs/docking_report.md` 第274行 vs. 其他所有文档

**矛盾**：
- docking_report.md 图2描述：**"RMSD = 0.78 Å（3592 atoms），整体折叠高度相似"**
- interface_analysis_report.md / project_log / af3_mutation_analysis.md：**"Hgal vs Hsap 全局 RMSD = 20.96 Å，物种间差异极大"**

0.78Å 几乎肯定是 PyMOL `align` 在某个核心局部区域（如 β-sheet）的对齐结果，而 20.96Å 是全长/CTD 全局 CA RMSD。在同一项目中将 0.78Å 描述为"整体折叠高度相似"是**严重误导**，与项目自己的其他结论直接冲突。

**修复建议**：
- 在 docking_report.md 中明确说明 0.78Å 是**局部 domain 对齐 RMSD**，而非全局
- 统一所有文档中的 RMSD 数值描述

---

## 🟡 中等问题（方法学缺陷或代码隐患）

### 4. 方法学混杂因素未充分控制

**问题**：
- Hgal WT 使用 **LightDock** pose，Hgal 4mut_rev 使用 **Rosetta** pose
- 两者直接比较时，界面差异（如 Hgal_4mut_rev 界面切换到 N-末端）**混杂了 docking 方法差异和突变效应**
- Hsap WT 和 Hsap 4mut 各只有 **1 个 replica**，统计力极弱

**后果**：README 中列出的"Finding 3/4/5"等关键发现，其因果归因（WT vs 4mut 的差异）可能部分甚至主要来源于 pose 差异，而非突变本身。

**修复建议**：
- 对 Hgal 4mut_rev 重新用 **LightDock** 对接（与 Hgal WT 相同方法），或对 Hgal WT 重新用 **Rosetta** 对接
- 若无法重做，需在论文/文档中将"不同 docking 方法"作为 **major limitation** 明确声明
- Hsap 系统至少各补跑 2 个 replica

---

### 5. tleap 脚本注释与代码不一致

**位置**：`build_system.py`（第108行注释）、`build_system_v2.py`（第70行注释）

**问题**：注释写"Neutralize and add 150 mM NaCl"，但实际代码是 `addIonsRand complex Na+ 0` 和 `addIonsRand complex Cl- 0`，这表示**只 neutralize，不加额外盐**。若要 150 mM NaCl，应使用 `addIonsRand complex Cl- 0.15`。

**修复建议**：统一注释与代码，或明确决定是否加盐。

---

### 6. 硬编码旧项目路径

**位置**：`scripts/rosetta_mutational_scan.py` 第218行

```python
base = "/Users/zhangzhuo/repos/personal/naked-mole-rat-cgas-trim41-simulation"
```

当前项目路径是 `/Users/zhangzhuo/repos/personal/cgas`，此脚本若直接运行会失败。

**修复建议**：改为相对路径或命令行参数传入。

---

## 🟢 轻微问题（可优化项）

### 7. 文档重复与混乱

- `project_log.md` 第99–105行，"待办"列表有**重复条目**（轨迹分析、MM-GBSA、Rosetta 扫描等出现两次）
- `interface_analysis_report.md` 有两个**"七、图表索引"**章节（第254行和第266行），明显是复制粘贴错误

---

### 8. MMPBSA.py 参数需确认

**位置**：`run_mmpbsa.py` 第87–89行

使用了 `-yr` 和 `-yl` 参数，但 AmberTools MMPBSA.py 的标准参数通常是 `-mr`（receptor mask）和 `-ml`（ligand mask）。需确认所用版本的语法。

---

### 9. Umbrella Sampling 力常数可能过大

US 脚本默认 k = 1000 kJ/mol/nm²。对于蛋白质体系的伞形采样，这个力常数**偏大**，可能导致每个窗口的采样被过度约束、收敛变慢，PMF 峰谷被抹平。通常建议 100–500 kJ/mol/nm²。

**建议**：用 k = 200–500 重新测试 1–2 个窗口，对比 CV 分布的宽度。

---

### 10. 蛋白质截断的生物学局限

所有 MD 系统仅包含 cGAS CTD (200-554/522) + TRIM41 SPRY (413-630)。这排除了：
- cGAS N-terminal DNA-binding domain（已知参与 TRIM41 调控）
- TRIM41 RING/B-box/coiled-coil domains

引用28（Zhen et al.）明确指出 TRIM41 的 coiled-coil domain 对相互作用很重要。截断可能丢失关键的 domain-domain 相互作用，影响结论的外推性。

**建议**：至少在 Discussion 中讨论截断的局限性。

---

## 📋 后续计划建议

### 立即修复（1–2 天）

| 任务 | 优先级 | 说明 |
|------|--------|------|
| 统一并清理突变位点文件 | P0 | 删除旧错误 fasta，更新脚本 residue 编号 |
| 修复统计方法 | P0 | 加入自相关时间计算 + block averaging，重算所有 p 值 |
| 修正文档矛盾 | P0 | 统一 RMSD 描述，明确 0.78Å 是局部对齐值 |
| 修复硬编码路径 | P1 | rosetta_mutational_scan.py 等 |
| 清理文档重复 | P1 | project_log.md、interface_analysis_report.md |

### 短期补充（1–2 周）

| 任务 | 优先级 | 说明 |
|------|--------|------|
| 补充 Replica | P0 | Hsap WT 和 Hsap 4mut 各至少再跑 2 个 replica |
| 控制 pose 混杂因素 | P0 | Hgal 4mut_rev 用 LightDock 重新对接，或明确声明 limitation |
| 完成 US 并优化 | P1 | 跑完剩余窗口，用 WHAM/MBAR 构建 PMF；测试更低 k 值 |
| 修正 tleap 盐浓度 | P1 | 统一注释与代码 |

### 中期深化（2–4 周）

| 任务 | 优先级 | 说明 |
|------|--------|------|
| 变构通路的因果链验证 | P0 | 计算 DCCM + PCA，从轨迹中提取 C 端 → N 端变构通路 |
| MM-GBSA 能量分解 | P1 | WT vs 4mut 界面结合自由能差异，per-residue 分解 |
| 全长结构考量 | P2 | AF3 multimer 全长预测（定性参考），或讨论截断局限 |

### 论文写作策略调整

- 早期"紧凑补丁（18Å vs 28Å）驱动泛素化"假说已被自己的 AF3 突变验证实验**推翻**
- 应将论文重心转向：**"物种特异性整体折叠差异 → 4mut 位点通过变构效应改变 N 端界面动力学 → 影响 TRIM41 识别/泛素化效率"**
- 强调这是"in silico 探索可能性空间"（因为原论文无结构数据），而非"验证原论文机制"
- 最有发表价值的发现是：**Hsap 4mut 引起 N 端 12Å 位移的变构效应**，值得作为核心 Figure

---

## 总体评价

项目的核心科学推理（从功能筛选 → 结构预测 → 发现变构效应）整体方向是正确的，但存在**统计方法错误、数据矛盾、突变位点清理不彻底**等严重问题。修复这些问题后，项目的变构效应假说（Hsap 4mut 引起 N 端 12Å 位移）是最有发表价值的发现，值得继续深入。
