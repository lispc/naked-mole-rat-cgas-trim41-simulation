# cGAS-TRIM41 MD 研究项目日志

> 本文档记录项目执行过程中的关键决策、数据、估算和推理，供后续论文撰写和复盘使用。
> 
> 最后更新：2026-04-22

---


## 一、AF3 结构预测结果

## 三、序列与突变映射

### 3.1 参考序列

| 蛋白 | UniProt | 长度 | 来源 |
|------|---------|------|------|
| 人 cGAS | Q8N884 | 522 aa | UniProt |
| 裸鼹鼠 cGAS | A0AAX6RS70 | 554 aa | UniProt (Heterocephalus glaber) |
| TRIM41 | Q8WV44 | 630 aa | UniProt |

**裸鼹鼠序列选择**：
- UniProt 上裸鼹鼠 cGAS 有多个注释版本（A0AAX6RS70 554 aa, A0AAX6RTF7 475 aa, A0AAX6RSF6 461 aa 等）
- 选择 **A0AAX6RS70**（554 aa）的原因：
  1. 长度最长，最可能包含完整功能域
  2. 论文中的突变位点编号（S463, E511, Y527, T530）在此序列上可直接对应
  3. 序列比对验证：保守区域 "NTGSYYEHVKI" 和 "GSPAVTLLI" 均存在

### 3.2 精确突变映射

通过全局序列比对（Biopython pairwise2，gap open=-10, gap extend=-0.5）确认：

| 论文标记 | 裸鼹鼠 aa | 人源对应位置 | 人源 aa | 突变方向 |
|---|---|---|---|---|
| S463 | S | **463** | C | C→S |
| E511 | E | **479** | K | K→E |
| Y527 | Y | **495** | L | L→Y |
| T530 | T | **498** | K | K→T |

**验证**：
- 人 cGAS 位置 463: C（在 `KDLGLCFDNCV` 上下文中）✅
- 人 cGAS 位置 479: K（在 `CLRTEKLENYF` 上下文中）✅
- 人 cGAS 位置 495: L（在 `LFSSNLIDKRS` 上下文中）✅
- 人 cGAS 位置 498: K（在 `SNLIDKRSKEF` 上下文中）✅

**重要说明**：论文中的编号是基于裸鼹鼠 cGAS 的坐标系统。由于裸鼹鼠序列比人源长 32 aa（主要在 N-端有延伸），导致位置 511→479、527→495、530→498 的偏移。这在论文 Methods 中需要明确说明。

---



## 二、AF3 复合物预测质量

## 六、风险预案

### 预案 A：AF3 预测置信度低（ipTM < 0.6）—— 已触发

**实际结果（全部 4 个 AF3 Job）**：

| Job | 体系 | ipTM (best) | pTM | 判定 |
|-----|------|------------|-----|------|
| 1 | Hsap_cGAS_WT + TRIM41 | 0.15 | 0.38 | ❌ 拒绝 |
| 2 | Hsap_cGAS_4mut + TRIM41 | 0.16 | 0.37 | ❌ 拒绝 |
| 3 | Hgal_cGAS_WT + TRIM41 | 0.23 | 0.37 | ❌ 拒绝 |
| 4 | Hgal_cGAS_4mut_rev + TRIM41 | 0.17 | 0.36 | ❌ 拒绝 |

所有复合物预测置信度均远低于 0.60 最低标准，这是 E3-底物瞬态相互作用的典型特征。

**单体结构质量**：
- 所有 job 的 Chain B (TRIM41): 630 residues, pLDDT ~62, ~54% 原子 ≥70
- Hsap cGAS (job1/2): 522 residues, pLDDT ~60
- Hgal cGAS (job3/4): 554 residues, pLDDT ~60
- **单体折叠可接受，但相对取向完全不可信**

**科学解释**：
1. TRIM41 无任何实验结构
2. E3 连接酶-底物互作通常是 transient/weak（AF 的已知弱点）
3. 互作可能依赖翻译后修饰或 DNA 损伤信号
4. 论文本身未直接测量物理结合强度

**已启动应对方案 A1：蛋白-蛋白对接 + 对接后 MD**

| 决策 | 内容 |
|------|------|
| 对接工具 | **ClusPro** (https://cluspro.bu.edu/) |
| 选择原因 | 无需注册，纯网页，用户可手动提交；HADDOCK/LightDock 无法本地自动化 |
| 受体 | TRIM41 WT 全长 或 SPRY 域 (413-630) |
| 配体 | cGAS 全长 或 C-端域 (~200-554) |
| 活性残基限制 | cGAS 的 4 个突变位点必须位于界面 |
| 打分模式 | Attraction（酶-底物互作优化） |

**TRIM41 结构域关键发现**（UniProt Q8WV44）：
- **B30.2/SPRY domain: residues 413-630** — 这是底物识别结构域！
- RING finger (N-端): E2 结合
- 多个无序区域: 51-86, 143-176, 503-538
- **这与我们的截取策略完全吻合**：SPRY 域负责识别 cGAS

**已准备文件**（全部 4 个 job）：

| Job | 目录 | cGAS 全长 | cGAS CTD | TRIM41 全长 | TRIM41 SPRY |


## 三、AF3 序列问题

通过检查 AF3 的 `job_request.json`：

| Job | 预期突变 | 实际序列 | 状态 |
|-----|---------|---------|------|
| Job 1 (Hsap_WT) | C, K, L, K | C, K, L, K | ✅ 正确 |
| **Job 2 (Hsap_4mut)** | **S, E, Y, T** | **C, K, L, K** | ❌ **WT 序列** |
| Job 3 (Hgal_WT) | S, E, Y, T | S, E, Y, T | ✅ 正确 |
| **Job 4 (Hgal_4mut_rev)** | **C, K, L, K** | **S, E, Y, T** | ❌ **WT 序列** |

**影响**：Job 2 和 Job 4 的 AF3 结构实际上是 WT，后续 MD 需要在 WT 基础上做 in-silico 突变。

echo "✅ af3_report.md created ($(wc -l < af3_report.md) lines)"
