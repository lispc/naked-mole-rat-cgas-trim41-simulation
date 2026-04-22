# ClusPro 蛋白-蛋白对接提交指南

> 因 AF3 复合物预测置信度极低（ipTM ~0.17），启动风险预案 A1：蛋白-蛋白对接 + 对接后 MD。
> 
> 工具选择：**ClusPro** (https://cluspro.bu.edu/)
> - 不需要注册，直接提交
> - 支持活性残基限制（提高精度）
> - 全自动，通常 1-4 小时出结果

---

## 背景：为什么需要对接

AF3 对 cGAS-TRIM41 复合物的 ipTM 仅 **0.17**（接受标准为 >0.60）。这意味着 AF3 无法可靠预测两个蛋白的相对结合取向。但两个蛋白各自的单体折叠结构是可信的（pLDDT ~60），因此需要用蛋白-蛋白对接来探索可能的结合模式。

---

## TRIM41 结构域信息（关键）

从 UniProt (Q8WV44) 获取：

| 结构域 | 范围 | 功能 |
|--------|------|------|
| RING finger | ~1-50 | E2 泛素结合酶结合 |
| B-box | ~50-100 | 辅助寡聚化 |
| Coiled-coil | ~100-400 | 寡聚化/定位 |
| **B30.2/SPRY** | **413-630** | **底物识别（关键！）** |

**结论**：TRIM41 的底物识别由 C-terminal SPRY 域（413-630）负责。cGAS 的 4 个突变位点（C-端催化域）很可能与 SPRY 域相互作用。

---

## 提交方案

### 方案 A：全长对接（推荐先尝试）

| 参数 | 设置 |
|------|------|
| Receptor | TRIM41 WT (630 aa) |
| Ligand | cGAS (522 or 554 aa) |
| 优势 | 不丢失任何可能的长程相互作用 |
| 劣势 | 计算量大，ClusPro 可能需要更长时间 |

### 方案 B：结构域对接（更聚焦）

| 参数 | 设置 |
|------|------|
| Receptor | TRIM41 SPRY 域 (413-630, 218 aa) |
| Ligand | cGAS C-terminal 域 (200-554, 355 aa) |
| 优势 | 计算快，聚焦关键互作区域 |
| 劣势 | 可能丢失 linker 或远程残基的辅助作用 |

**建议**：两个方案都提交，比较结果一致性。

---

## 活性残基限制（Active Residues）

ClusPro 支持 "Attraction" 模式，可以指定哪些残基必须位于界面上。这基于生物学先验知识。

### cGAS 上的活性残基（4 个突变位点）

| 体系 | 活性残基位置 |
|------|-------------|
| Hsap cGAS WT / 4mut | **463, 479, 495, 498** |
| Hgal cGAS WT / 4mut_rev | **463, 511, 527, 530** |

说明：这 4 个位点在实验中被证明是 TRIM41 介导泛素化的功能关键位点，极可能位于或靠近结合界面。

### TRIM41 上的活性残基

SPRY 域的表面残基（基于 UniProt 和结构预测推断）：
- 选取 SPRY 域中溶剂可及的、带电荷的残基作为潜在活性位点
- 具体列表需提交前用 PDB 检查表面暴露度

**保守策略**：先不指定 TRIM41 的活性残基（让 ClusPro 全局搜索），仅限制 cGAS 的 4 个位点必须参与界面。

---

## 已准备好的提交文件

文件位置：`structures/af3_raw/job4_Hgal_4mut_rev/`

| 文件名 | 内容 | 用途 |
|--------|------|------|
| `cgas_fixed.pdb` | cGAS 单体（加氢、补全缺失原子）| 全长对接 |
| `trim41_fixed.pdb` | TRIM41 单体（加氢、补全缺失原子）| 全长对接 |
| `cgas_CT_200-554.pdb` | cGAS C-端域截取 | 结构域对接 |
| `trim41_SPRY_413-630.pdb` | TRIM41 SPRY 域截取 | 结构域对接 |

---

## 逐步提交步骤

### Step 1: 打开 ClusPro
访问 https://cluspro.bu.edu/

### Step 2: 上传受体（Receptor）
- 点击 "Receptor" 上传框
- 选择 `trim41_fixed.pdb`（全长）或 `trim41_SPRY_413-630.pdb`（结构域）
- 如果知道活性残基，在 "Active Residues" 框中输入（格式：每行一个残基编号）

### Step 3: 上传配体（Ligand）
- 点击 "Ligand" 上传框
- 选择 `cgas_fixed.pdb`（全长）或 `cgas_CT_200-554.pdb`（结构域）
- **关键**：在 "Active Residues" 框中输入 cGAS 的 4 个突变位点：
  ```
  463
  511
  527
  530
  ```
  （对于 Hgal cGAS；若用人源则输入 463, 479, 495, 498）

### Step 4: 选择打分模式
- **推荐**："Attraction"（对于酶-底物互作，惩罚排斥，奖励吸引）
- 备选："Balanced"（如果 Attraction 结果不合理）

### Step 5: 提交
- 点击 "Submit"
- 记录 Job ID
- 等待邮件通知（通常 1-4 小时）

### Step 6: 下载结果
- 下载 top 10 poses（PDB 格式）
- 放入 `structures/docking/cluspro/` 对应目录

---

## ⚠️ 关于截断后残基编号的关键说明

**截断后的 PDB 文件保留了原始序列的残基编号**。

例如：
- cGAS CTD (200-554) 的 PDB 中，第一个残基编号是 200，最后一个是 554
- 突变位点 **463, 511, 527, 530 在截断 PDB 中仍然是这些编号**
- TRIM41 SPRY (413-630) 同理

**结论**：无论提交全长还是截断结构，ClusPro 的 "Active Residues" 框中输入的编号**完全相同**，不需要转换。

---

## 需要提交的 Job 列表

### 现在就可以提交的（job4 AF3 结果已处理完毕）

| Job | Receptor | Ligand | cGAS 活性残基 | 状态 |
|-----|----------|--------|--------------|------|
| A | TRIM41_WT_full (630 aa) | Hgal_cGAS_4mut_rev_full (554 aa) | 463, 511, 527, 530 | **🟢 现在可提交** |
| B | TRIM41_SPRY_413-630 (218 aa) | Hgal_cGAS_CT_200-554 (355 aa) | 463, 511, 527, 530 | **🟢 现在可提交** |

### 等其余 3 个 AF3 job 完成后提交

| Job | Receptor | Ligand | cGAS 活性残基 | 状态 |
|-----|----------|--------|--------------|------|
| C | TRIM41_WT_full | Hsap_cGAS_WT_full (522 aa) | 463, 479, 495, 498 | 🟡 等待 job1 |
| D | TRIM41_WT_full | Hsap_cGAS_4mut_full (522 aa) | 463, 479, 495, 498 | 🟡 等待 job2 |
| E | TRIM41_WT_full | Hgal_cGAS_WT_full (554 aa) | 463, 511, 527, 530 | 🟡 等待 job3 |
| F | TRIM41_SPRY | Hsap_cGAS_CT | 463, 479, 495, 498 | 🟡 等待 job1/2 |

---

## 现在就可以提交的 2 个 Job（详细步骤）

### Job A: Hgal_cGAS_4mut_rev (全长) + TRIM41_WT (全长)

| 参数 | 值 |
|------|-----|
| Receptor PDB | `structures/af3_raw/job4_Hgal_4mut_rev/trim41_fixed.pdb` |
| Ligand PDB | `structures/af3_raw/job4_Hgal_4mut_rev/cgas_fixed.pdb` |
| Ligand Active Residues | 463, 511, 527, 530 |
| Scoring mode | Attraction |

### Job B: Hgal_cGAS_CT (截断) + TRIM41_SPRY (截断)

| 参数 | 值 |
|------|-----|
| Receptor PDB | `structures/af3_raw/job4_Hgal_4mut_rev/trim41_SPRY_413-630.pdb` |
| Ligand PDB | `structures/af3_raw/job4_Hgal_4mut_rev/cgas_CT_200-554.pdb` |
| Ligand Active Residues | 463, 511, 527, 530 |
| Scoring mode | Attraction |

**注意**：两个 job 的活性残基编号完全相同（463, 511, 527, 530），因为截断 PDB 保留了原始编号。

---

## 对接后处理流程

1. **提取 top 10 poses**
2. **RMSD 聚类**：用 `scripts/analyze_docking.py` 去除冗余 pose
3. **能量最小化**：每个 pose 用 OpenMM 做 1000 steps minimization
4. **短 MD 弛豫**：每个 pose 做 5-10 ns NPT 平衡
5. **稳定性筛选**：观察 RMSD，保留收敛的 pose（RMSD < 3Å）
6. **长 MD**：对最优的 1-3 个 pose 做 200 ns 生产模拟
7. **MM-GBSA**：计算每个 pose 的结合能

---

## 质量控制

- 若 ClusPro 的 top pose 中，cGAS 活性残基与 TRIM41 SPRY 域距离 >10Å → 对接失败，需调整参数或换工具
- 若 top 10 poses 的 RMSD 差异 < 2Å → 结合模式单一，可信度高
- 若 top 10 poses 分成 2-3 个明显簇 → 可能有多种结合模式，需分别测试

---

*文档创建：2026-04-22*
*对应风险预案 A1*
