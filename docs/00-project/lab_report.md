# cGAS-TRIM41 4mut 计算研究：实验报告

> 项目周期：2026-04-22 — 2026-05-08  
> 目标：通过计算方法解释裸鼹鼠 cGAS 四个氨基酸变异如何影响 TRIM41 介导的泛素化

---

## 一、项目概览

### 1.1 核心科学问题

Chen et al. (Science, 2025) 发现裸鼹鼠 cGAS 的四个氨基酸变异（人源：D431S/K479E/L495Y/K498T）改变 TRIM41 介导的泛素化效率，进而影响 DNA 修复和衰老。**这些变异如何影响 cGAS-TRIM41 相互作用的分子机制完全未知。**

### 1.2 实验演进路线

```
原始假设：4mut 在 TRIM41 界面上，直接改变结合
    ↓
发现：界面在 N 端，4mut 在 C 端，距离 30-39 Å（推翻原始假设）
    ↓
新假设：4mut 通过长程变构效应影响 N 端界面
    ↓
四系统 MD：确认结合亲和力不变，但动态改变
    ↓
四元复合物：SPRY 锚定 cGAS，K315 是唯一可达靶点
    ↓
US/PMF：定量 4mut 的自由能偏移
```

### 1.3 最终产出的计算实验清单

| 实验 | 规模 | 结论 |
|------|------|------|
| AF3 结构预测 | 5 jobs | cGAS 单体高置信度，复合物 ipTM<0.25 |
| 蛋白对接 | LightDock + Rosetta + ClusPro | 界面在 N 端，4mut 不在界面上 |
| 四系统 MD | 4×3×200ns = 2.4μs | Hsap 4mut 不改变结合几何；Hgal 4mut_rev 破坏界面 |
| 磷酸化 MD | S305-phos + S305E 共 1.2μs | 完全解离（side story，未纳入主线） |
| GROMACS 验证 | 200ns ×1 | CMAP bug 发现与修复 |
| 四元 FULL MD | 50ns ×2 | SPRY 锚定 cGAS，K315 零漂移 |
| 全 Lys 扫描 | 38 Lys ×2 系统 | K315 唯一可达靶点 |
| Umbrella Sampling | 21 窗口 ×10ns = 210ns | PMF 最小值偏移 −2.5 Å |
| Boltz-2 cGAS+DNA | 5 models | K347 在所有模型中均背对 SPRY |
| CC40 MD | 50ns ×4 | K347 不可达，接受 K315 为唯一靶点 |

---

## 二、关键决策与转折点

### 2.1 界面不在 C 端（04-23）

**背景**：原始假设认为 4mut 位点在 TRIM41 结合面上，通过直接改变结合强度影响泛素化。

**发现**：LightDock + Rosetta 界面分析显示 80% 接触在 N 端（200-299），C 端（451-554）**零接触**。4mut 位点距界面 30-39 Å。

**影响**：整个项目方向从"界面结合"转向"变构效应"。

### 2.2 突变编号纠正（04-24）

**错误**：项目初期将人源 4mut 标记为 C463S/K479E/L495Y/K498T。

**纠正**：通过全局序列比对确认正确编号为 **D431S**/K479E/L495Y/K498T。

**影响**：早期基于 C463S 的 fasta、脚本、部分分析需要修正。Hgal 4mut_rev RMSD 从 13.04 Å（错误）修正为 0.98 Å（正确）。

### 2.3 S305 磷酸化 side story 的扩展与收敛（04-30 ~ 05-02）

**启动原因**：Zhen et al. 文献指出 CHK2 磷酸化 S305 促进 cGAS-TRIM41 相互作用。

**发现**：S305-phos 导致**完全解离**（与文献相反）。S305E 显示异质性行为。

**决策**：跑完 200ns×3 确认结论后**立即停止扩展**。共投入 ~1.2μs，在项目中占比过高。最终在论文中降为 1 小节。

### 2.4 GROMACS CMAP bug（05-01）

**发现**：parmed 将 ff19SB 的 14 种残基特异性 CMAP 压缩为 1 种，导致 GROMACS RMSD 被高估 4 倍。

**修复**：改用 GROMACS 2026 原生 amber19sb.ff，RMSD ratio 降至 1.2×。

**教训**：跨引擎验证时，力场转换工具可能引入系统性错误。始终用原生实现作为 gold standard。

### 2.5 Python 对齐 bug（05-02）

**发现**：`rotation_matrix()` 在未中心化坐标上计算，GROMACS RMSD 被高估 82%。

**修复**：改用 MDAnalysis `alignto()` 和 `gmx rms` 原生工具。

**教训**：分析脚本在用于正式结论前必须用独立方法验证。

### 2.6 四元复合物迭代（05-03 ~ 05-06）

**v1**（E2~Ub + cGAS，05-03）：
- E2~Ub 在 MD 中完全散开（tleap 未识别异肽键）
- K315→Ub 保持在 27-29 Å
- **失败原因**：缺少 RING 域稳定 E2~Ub closed 构象；缺少异肽键约束

**v2**（+RING + 异肽 bond restraint，05-06）：
- E2~Ub 稳定在 ~2.9 Å
- K315→Ub 降至 11-12 Å（初期）→漂移到 20-28 Å
- **失败原因**：缺少 SPRY 锚定 cGAS

**FULL**（+SPRY + CC linker，05-06）：
- K315 稳定在 ~20 Å，零漂移
- CC linker ~79 Å
- **成功**作为 US 的 baseline 模型

### 2.7 COM restraint 水原子 bug（05-06）

**发现**：CENTROID 计算包含了 12 万水原子，导致速度从 122 降到 29 ns/day。

**修复**：改用 protein-only atom ranges，速度恢复到 121 ns/day。

**教训**：性能异常时首先检查 force group 的原子选择。

### 2.8 K347 靶点探索（05-07）

**触发**：文献报告 TRIM41 泛素化 K347。我们一直聚焦的 K315 可能不是真正靶点。

**探索**：
1. 全 Lys 扫描：K315 唯一 <25 Å，K347 在 59 Å
2. cGAS-DNA 二聚体（4LEZ）：K347 仍在对立面
3. Boltz-2 cGAS+DNA+TRIM41：5 个 model 全部显示 SPRY 近 K315 远 K347
4. CC40 MD：缩短 CC 到 40-100 Å，K347 仍不可达（min 50 Å）

**结论**：K347 在 cGAS 拓扑上不可达。接受 K315 为唯一靶点。

---

## 三、技术细节与教训

### 3.1 Bug 清单

| Bug | 发现日期 | 影响 | 修复 |
|-----|---------|------|------|
| C463S→D431S 编号错误 | 04-24 | 早期 fasta/脚本/部分分析错误 | 全局序列比对纠正 |
| GROMACS CMAP 压缩 | 05-01 | RMSD 高估 4× | 原生 amber19sb.ff |
| Python 对齐 bug | 05-02 | RMSD 高估 82% | MDAnalysis alignto() |
| t-test 对自相关数据误用 | 04-27 | p 值不可信 | correlated t-test |
| PBC 包裹误判 | 05-01 | COM 虚报 | gmx trjconv -pbc mol |
| 4mut_rev NaN 崩溃 | 05-02 | MD 无法启动 | 深度重新最小化 |
| auto-launcher race condition | 05-03 | 重复启动进程 | 停用，改手动 |
| COM restraint 水原子 | 05-06 | 速度降 4× | protein-only selection |
| 异肽键 PDB 距离 2.41 Å（非共价） | 05-06 | E2~Ub 在 MD 中散开 | harmonic restraint |
| 磁盘瞬时满（DCD 写入失败） | 05-07 | 4×CC40 崩溃 | 重启恢复 |

### 3.2 方法学经验

1. **多方法交叉验证是不可替代的。** AF3、Boltz-2、Chai-1 对 cGAS-TRIM41 界面的预测定性不同（Jaccard=0.00），单一模型无法给出可靠结论。

2. **跨引擎验证救了项目。** GROMACS vs OpenMM 的不一致引导我们发现了 CMAP bug 和对齐 bug。没有这一步，大量 GROMACS 数据会被错误使用。

3. **Replica 设计至关重要。** S305E 的 3 个 replica 显示了从完全解离到稳定结合的极端异质性——单 replica 会给出完全误导的结论。

4. **Steered MD 的阻力本身就是 PMF 信号。** K315 被拉到 ~13 Å 后无法再前进——这个"阻力墙"的方向和高度在 WHAM 分析中得到了定量确认。

5. **PDB 格式是计算结构生物学的永恒痛点。** Bio.PDB 列对齐问题在项目中出现了至少三次。应该有一个经过列验证的 PDB writer。

### 3.3 资源分配反思

| 实验 | GPU 时间 | 论文贡献 | 性价比 |
|------|---------|---------|--------|
| 四系统 MD | ~2.4μs | 核心 Table 3 | 高 |
| S305 磷酸化 | ~1.2μs | 1 小节 side note | **低** |
| 四元 FULL + US | ~0.6μs | 核心 §2.6-2.7 | 高 |
| K347 探索 | ~0.3μs | negative result | 中 |
| CC40 | ~0.2μs | negative result | 低 |

**最大浪费**：S305 磷酸化。结论在 200ns 时已明确（完全解离），但跑了 600ns+600ns。"知道什么时候停"是计算项目最重要的技能之一。

---

## 四、文件索引

- 项目主日志：`docs/00-project/project_log.md`（§1-65）
- 历史日志：`docs/00-project/project_log_2026_04.md`
- 项目回顾：`docs/00-project/retrospective.md`
- 论文终稿：`paper/paper_final.md`
- 文档索引：`docs/README.md`

---

*最后更新：2026-05-08*
