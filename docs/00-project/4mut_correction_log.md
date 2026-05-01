# 4mut 位点修正记录

> 日期：2026-04-23
> 修正人：Kimi Code CLI（基于论文原文重新确认）

---

## 问题发现

通过 `pdftotext main.pdf` 重新提取论文原文，发现项目内部 4mut 定义存在严重错误。

### 论文原文（page 2）

> "Introduction of the four amino acid substitutions into **naked mole-rat cGAS (S463D+E511K+Y527L+T530K)** significantly diminished its stimulatory effect on HR repair, whereas the corresponding 4-aa mutations in **human cGAS (D431S+K479E+L495Y+K498T)**"

### 修正前的错误定义

| 来源 | 错误定义 | 问题 |
|------|---------|------|
| `sequences/Hsap_cGAS_4mut.fasta` | C463S, K479E, L495Y, K498T | ❌ 463 位点错误（应为 431） |
| `structures/af3_raw/job4_Hgal_4mut_rev/` | 与 WT 完全相同 | ❌ 未做任何突变 |

### 修正后的正确定义

| 物种 | 突变 | 方向 | 功能 |
|------|------|------|------|
| **人 cGAS** | **D431S + K479E + L495Y + K498T** | WT → 裸鼹鼠样 | 获得 HR 刺激 |
| **裸鼹鼠 cGAS** | **S463D + E511K + Y527L + T530K** | WT → 人样 | 丧失 HR 刺激 |

### 位点对应关系

由于 Hgal N 端插入导致编号偏移：
- Hgal **463** = Hsap **431**
- Hgal **511** = Hsap **479**
- Hgal **527** = Hsap **495**
- Hgal **530** = Hsap **498**

---

## 修正文件

| 文件 | 说明 |
|------|------|
| `sequences/Hsap_cGAS_4mut_corrected.fasta` | 正确的人 4mut：D431S, K479E, L495Y, K498T |
| `sequences/Hgal_cGAS_4mut_rev_corrected.fasta` | 正确的裸鼹鼠逆转：S463D, E511K, Y527L, T530K |

---

## 对之前分析的影响

### 需要重新做的分析

1. **Hsap 4mut 变构效应验证**
   - 之前基于错误结构（C463S+K479E+L495Y+K498T）
   - 正确的 4mut 位点（431/479/495/498）在错误结构中位移几乎为 0
   - **必须用修正后的序列重新做 AF3 预测**

2. **4mut 到界面距离**
   - 之前基于 463/511/527/530 计算的距离无效
   - 需要重新基于 431/479/495/498 计算

3. **Hgal MD 系统**
   - 现有 Hgal MD 模拟的是 WT（S463/E511/Y527/T530）
   - 这是正确的裸鼹鼠天然状态
   - 但需要额外构建"人样"突变体（S463D+E511K+Y527L+T530K）的 MD

### 仍然有效的分析

1. **界面位置（N 端）**
   - Hgal 和 Hsap 的界面都在 N 端/中间区域
   - 4mut 位点远离界面的结论仍然成立（只是具体距离数值需要更新）

2. **Hgal vs Hsap 结构差异**
   - 全局 RMSD 20.96Å 仍然有效
   - 物种特异性界面的结论仍然成立

---

## 下一步

1. [ ] 用修正后的 fasta 重新提交 AF3 单体预测（Hsap 4mut + Hgal 4mut_rev）
2. [ ] 基于新的 AF3 结构重新评估变构效应
3. [ ] 更新所有文档中的 4mut 位点引用
4. [ ] 重新计算 4mut 到界面的距离
