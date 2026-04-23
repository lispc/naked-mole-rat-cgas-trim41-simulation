# Nice-to-Have Items（未来扩展计划）

> 本文档记录当前阶段未执行但值得在未来补充的实验/分析项目。
> 
> 原则：**Hgal 数据足以支撑论文发表**，以下项目为增强论证或拓展深度。

---

## 目录

| Item | 优先级 | 估计工作量 | 论文价值 |
|------|--------|-----------|---------|
| [1. Hsap 4mut 结构获取](#1-hsap-4mut-结构获取) | P1 | 2-3 天 | 直接对照实验（WT vs 突变） |
| [2. Hgal_4mut_rev 结构验证](#2-hgal_4mut_rev-结构验证) | P2 | 1-2 天 | 反向验证突变效应可逆 |
| [3. 全长蛋白 MD](#3-全长蛋白-md) | P3 | 1 周 | 验证结构域截断未引入偏差 |
| [4. 更多 MD 重复](#4-更多-md-重复) | P2 | 线性 | 提高统计显著性 |
| [5. 不同力场/水模型验证](#5-不同力场水模型验证) | P3 | 1 周 | 方法学稳健性 |
| [6. 全长蛋白对接](#6-全长蛋白对接) | P4 | 3-5 天 | 完整复合物结构 |
| [7. 温度/变性模拟](#7-温度变性模拟) | P4 | 2-3 天 | 热稳定性差异 |

---

## 1. Hsap 4mut 结构获取

### 现状
- AF3 Job 2 实际提交的是 Hsap_WT 序列（序列 bug）
- 因此目前没有 Hsap_cGAS_4mut 的预测结构
- LightDock 对 Hsap_WT 的对接结果：0/25 成功（几何约束，28.6Å 间距）

### 为什么重要
- **最直接的突变效应验证**：将 Hsap_WT 突变为 4mut 后，如果 495/498 位置的 loop 发生构象变化，可能改变空间几何
- 如果 Hsap_4mut 也形成紧凑补丁 → 支持"突变驱动几何改变"的假说
- 如果 Hsap_4mut 仍保持分散 → 说明还有其他因素

### 方法

#### 方案 A：PyMOL in-silico 突变（最快，已具备条件）

```bash
conda activate py311
pymol

# PyMOL 命令
load structures/af3_raw/job2_Hsap_4mut/cgas_fixed.pdb  # 实际这是 WT
# 或使用 domain truncation:
load structures/af3_raw/job2_Hsap_4mut/cgas_CT_200-554.pdb

# 4 个突变
wizard mutagenesis
# C463S, K479E, L495Y, K498T
# 选择最优 rotamer（避免 clash）
# save as hsap_4mut_model.pdb
```

**局限性**：
- PyMOL 的 rotamer 选择基于统计库，不能预测 loop 重排
- 突变后区域（特别是 495/498 附近）可能发生显著构象变化
- 适合作为 MD 起始点（MD 会自然弛豫），但不适合作为"最终结构"

**适用场景**：快速获得起始结构用于 MD

#### 方案 B：AF3 重新提交（最可靠，需等待）

需要重新提交 2 个 job：
- Hsap_cGAS_4mut + TRIM41_WT
- Hgal_cGAS_4mut_rev + TRIM41_WT

**序列确认**：
```
> Hsap_cGAS_4mut
[完整序列，确认包含 C463S, K479E, L495Y, K498T]

> Hgal_cGAS_4mut_rev
[完整序列，确认包含 S463C, E511K, Y527L, T530K]
```

**时间成本**：
- AF3 Server 排队：未知（通常 1-3 天）
- 加上结果处理：+0.5 天

**适用场景**：如果 reviewer 要求更严格的结构基础

#### 方案 C：ColabFold/AlphaFold2 本地（备选）

如果 AF3 Server 排队太长，可考虑 ColabFold（免费 GPU）：
- https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
- 支持 multimer，但精度略低于 AF3
- 适合快速验证

### 当前决策
**暂缓**。Hgal WT（即实际上的"突变体"）已经提供了充分的结构证据：
1. 18.4Å 紧凑补丁与 28.6Å 分散位点的对比
2. Docking 成功 vs 失败的几何解释
3. MD 轨迹将展示动态行为

如果后续 reviewer 要求对照实验，再执行方案 A 或 B。

---

## 2. Hgal_4mut_rev 结构验证

### 现状
- AF3 Job 4 也提交了 WT 序列（同 Job 2 的 bug）
- 没有真正的 Hgal_cGAS_WT + 反向突变结构

### 为什么重要
- **反向验证**：将裸鼹鼠突变回人类序列，观察空间几何是否恢复
- 如果反向突变后恢复分散 → 强支持突变是充分必要条件
- 如果反向突变后仍保持紧凑 → 提示其他物种特异性因素

### 方法
同 Hsap 4mut（PyMOL 或 AF3 重新提交）

### 当前决策
**暂缓**。与 Hsap 4mut 类似，属于增强论证而非必需。

---

## 3. 全长蛋白 MD

### 现状
- 当前 MD 使用结构域截断：cGAS CTD (200-554) + TRIM41 SPRY (413-630)
- 全长 cGAS = NTase (1-200) + CTD (200-554)，共 522/554 aa
- 全长 TRIM41 = RING + B-box + coiled-coil + SPRY，共 630 aa

### 为什么重要
- 结构域截断可能：
  a) 改变了蛋白间的相对取向
  b) 丢失了远距离的别构效应
  c) 影响了柔性 loop 的行为

### 方法

#### 方案 A：使用 AF3 全长结构（低置信度但可用）
- AF3 预测的复合物 ipTM <0.25，但单体 pLDDT ~60
- 可以提取单体后重新对接（ClusPro 全长已失败，需其他工具）
- 或直接用 AF3 的低分复合物作为 MD 起始，让 MD 弛豫

#### 方案 B：Rosetta docking（替代蛋白对接）
```bash
# Rosetta docking protocol
rosetta_scripts.default.macosclangrelease \
  -in:file:s full_length_complex.pdb \
  -parser:protocol docking.xml
```

#### 方案 C：从截断结构外推
- 假设 NTase 域和 RING/B-box 域不影响 CTD-SPRY 界面
- 在论文中明确说明此假设

### 当前决策
**暂不执行**。在论文 Methods 中明确声明：
> "由于 AF3 全长复合物预测置信度低，我们使用截断结构域（cGAS CTD 200-554 + TRIM41 SPRY 413-630）进行 MD 模拟。SPR Y域是 TRIM41 的底物识别域，CTD 包含全部 4 个突变位点，截断保留了功能核心。"

如果 reviewer 质疑，再执行全长 MD。

---

## 4. 更多 MD 重复

### 现状
- 计划：200ns × 3 重复
- 当前：0 个完成

### 为什么重要
- MD 是随机过程，需要重复来评估收敛性和统计显著性
- 3 重复是发表论文的最低要求
- 更多重复（5-10）可提高 p 值可信度

### 方法
- 从同一 minimized 结构出发，使用不同随机种子
- 或从 heating 的不同时间点生成不同起始构象

### 执行状态
**待执行**。Hgal_domain rep1 即将启动。

---

## 5. 不同力场/水模型验证

### 现状
- 当前：ff19SB + OPC

### 为什么重要
- 力场选择可能影响定量结果（结合能、RMSF 等）
- 不同水模型（TIP3P vs OPC）对界面水分子行为描述不同

### 可选组合

| 力场 | 水模型 | 特点 |
|------|--------|------|
| ff19SB | OPC | 当前选择，最新蛋白力场 |
| ff14SB | TIP3P | 经典组合，文献最常用 |
| ff19SB | TIP4P-Ew | 4 点水模型，更精确 |

### 当前决策
**暂不执行**。ff19SB + OPC 是 2020s 的标准选择，足以发表论文。如果 reviewer 要求，再跑一个 ff14SB+TIP3P 的对照。

---

## 6. 全长蛋白对接

### 现状
- ClusPro 全长对接失败（蛋白相距 80+Å）
- 未尝试其他全长对接工具

### 为什么重要
- 验证截断结构的对接 pose 是否与全长一致
- 完整复合物结构更有说服力

### 方法
- HDOCK (http://hdock.phys.hust.edu.cn/)
- RosettaDock
- ZDOCK

### 当前决策
**暂不执行**。结构域对接已经通过了活性残基验证（20/20 poses），截断是合理的。

---

## 7. 温度/变性模拟

### 为什么重要
- 验证复合物热稳定性差异
- 观察高温下界面是否解离

### 方法
- Replica Exchange MD (REMD) — 计算成本高
- 渐进升温 MD (steered MD) — 较简单
- 不同温度下的常规 MD（300K, 350K, 400K）

### 当前决策
**暂不执行**。属于深度分析，可在论文 revision 阶段补充。

---

## 总结：优先级与决策矩阵

| 优先级 | Item | 触发条件 | 预计投入 |
|--------|------|---------|---------|
| **P0** | Hgal MD 生产运行 | 立即执行 | 4-5 天/重复 |
| **P1** | Hsap 4mut 结构 | Reviewer 要求对照 | 2-3 天 |
| **P2** | 更多 MD 重复 | 按原计划执行 | 线性 |
| **P2** | Hgal_4mut_rev | Reviewer 要求反向验证 | 1-2 天 |
| **P3** | 全长 MD | Reviewer 质疑截断 | 1 周 |
| **P3** | 力场对照 | Reviewer 质疑方法 | 1 周 |
| **P4** | 全长对接 | Reviewer 要求 | 3-5 天 |
| **P4** | 温度模拟 | 拓展研究 | 2-3 天 |

**核心原则**：
> "以最少的工作量讲清楚科学问题。Hgal 的 18Å vs 28Å 几何差异已经足够有力。"
