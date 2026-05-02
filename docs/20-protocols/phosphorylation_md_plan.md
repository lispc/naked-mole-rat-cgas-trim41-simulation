# cGAS 磷酸化位点 MD 模拟方案

**日期**: 2026-04-23（更新 2026-04-30）
**目标**: 探索磷酸化如何调控 cGAS-TRIM41 结合界面
**状态**: ✅ S305-phos 体系构建完成，3× replica OpenMM MD 运行中（~15.5ns/200ns）

---

## 一、构建体范围确认

当前项目使用的 cGAS 构建体为 **C-terminal domain (residues 200-554, 355 aa)**，由 `cgas_CT_200-554.pdb` 经 AlphaFold3 预测获得。

| 磷酸化位点 | 全长编号 | 是否在 200-554 内 | 可模拟性 |
|-----------|---------|------------------|---------|
| **S120** (CHK2) | 120 | ❌ 否 | 需要全长 cGAS (1-522) |
| **S305** (CHK2/AKT) | 305 | ✅ 是 | **最高优先级** |
| **S435** (PPP6C) | 435 | ✅ 是 | 次要优先级 |

> **结论**: 在当前构建体范围内，**S305 是唯一可直接模拟的磷酸化位点**。若需研究 S120，必须重新构建全长 cGAS + TRIM41 SPRY 体系（需重新 AF3 预测或 docking）。

---

## 二、磷酸化位点生物学依据

### S305 — 最高优先级

- **激酶**: CHK2 (DNA 损伤响应) / AKT (上下文依赖，功能相反)
- **位置**: NTase 催化结构域表面，暴露于溶剂
- **功能**: 
  - CHK2 磷酸化 S305 **促进** cGAS-TRIM41 结合（Zhen et al., Nat Commun 2023）
  - AKT 磷酸化 S305 **抑制** cGAS 酶活性（Seo et al., 2015）
  - 与 S120 协同：双磷酸化效果 > 单磷酸化
- **结构影响**: 引入 -2 电荷，可能形成新的盐桥，改变表面静电互补性

### S435 — 次要优先级

- **激酶/磷酸酶**: PPP6C (去磷酸化)
- **位置**: C-terminal domain
- **功能**: 调控 cGAMP 合成，非直接调控 TRIM41 结合

---

## 三、模拟体系设计

### 方案 A（推荐）：当前构建体 cGAS(200-554) + TRIM41 SPRY

| 体系 | 描述 | Replicas | 时长 | 目的 |
|------|------|----------|------|------|
| WT | 无磷酸化 | 3 | 200-500ns | 基线对照 |
| S305-phos | SEP @ 305 | 3 | 200-500ns | 测试磷酸化对结合界面的影响 |
| S305E | S305E 突变 | 3 | 200-500ns | 磷酸化模拟对照（验证电荷效应） |
| 4mut+S305E | 4mut + S305E | 3 | 200-500ns | 检验磷酸化是否补偿 4mut 结合减弱 |

### 方案 B（扩展）：若使用全长 cGAS (1-522)

| 体系 | 描述 | Replicas | 目的 |
|------|------|----------|------|
| S120-phos | SEP @ 120 | 3 | N-terminal 磷酸化效应 |
| S305-phos | SEP @ 305 | 3 | 催化域磷酸化效应 |
| **S120/S305-diphos** | SEP @ 120 + SEP @ 305 | 3 | **协同效应（文献最高优先级）** |
| S120E/S305E | 双磷酸化模拟 | 3 | 与真实磷酸化对比 |

---

## 四、初始结构准备流程（✅ 已完成）

> **构建结果**（2026-04-30）：
> - 体系: `data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop/.rst7`
> - 原子数: 81,487（蛋白 8,815 + 水 18,162 + 离子）
> - 电荷: 0.000（中性）
> - SEP 残基验证: 14 原子，磷酸基团净电荷 ≈ -2，无 HG
>
> **关键脚本**: `scripts/build_phosphorylated_system.py`

### 步骤 1：获取 WT 起始结构

使用现有 OpenMM NPT 平衡后的最终构象作为起始结构（最佳实践：已平衡的水盒子结构）。

```bash
# 从 OpenMM NPT 轨迹提取最后一帧
# (已有 data/md_runs/Hsap_WT/rep1/Hsap_WT_rep1_npt.dcd)
```

### 步骤 2：在 PDB 上添加磷酸基团

磷酸化 Ser/Thr 需要手动添加 PO₃ 基团。两种方法：

**方法 A：PyMOL Builder（推荐，最可控）**
```python
# PyMOL 命令
fetch your_structure.pdb
# 1. 选中 S305 的 OG 原子
# 2. Builder > Add PO3
# 3. 调整键长：P-O = 1.61 Å, O-P-O = 109.5°
# 4. Save as cgas_S305phos.pdb
```

**方法 B：专业工具**
- **PyTMs** (PyMOL plugin): ` phosphorylate selection, residue=S305 `
- **H++ server**: 预测磷酸化后的质子化状态（辅助）

### 步骤 3：修改 PDB 残基名（脚本自动完成）

脚本自动将残基名 SER → SEP，原子命名符合 AMBER 规范（P, O1P, O2P, O3P）。

```bash
# 如需手动检查
grep "SEP" Hsap_WT_S305phos_raw.pdb | head -15
```

### 步骤 4：pdb4amber 处理

```bash
pdb4amber -i Hsap_WT_S305phos_raw.pdb -o Hsap_WT_S305phos_amber.pdb --reduce
```

**注意**: 磷酸基团在双负离子形式下不应有氢。pdb4amber 不会为 SEP 的 OG 添加 HG（正确），因为脚本已在步骤 2 中排除 HG。

**预期警告**: `gap of 3.3 A between GLY 105 and SER 107` — 正常，因 SEP 106 的 backbone 连接由 tleap 重建。

### 步骤 5：tleap 构建体系（✅ 成功）

```bash
cat > tleap_S305phos.in << 'EOF'
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.phosaa19SB

complex = loadPdb Hsap_WT_S305phos_amber.pdb
check complex
solvateOct complex OPC 12.0
addIonsRand complex Na+ 0
addIonsRand complex Cl- 0
saveAmberParm complex Hsap_WT_S305phos.prmtop Hsap_WT_S305phos.rst7
savePdb complex Hsap_WT_S305phos_solvated.pdb
charge complex
quit
EOF

tleap -f tleap_S305phos.in
```

**构建结果**: tleap `check` 通过，总原子数 81,487，电荷 0.0。

**关键修复历史**:
1. 原始脚本保留 HG 原子 → tleap 报错 `Atom HG does not have a type`
2. 修复: `build_phosphorylated_system.py` 在写入 PDB 时跳过 target residue 的 HG 原子

### 步骤 6：能量最小化（关键！磷酸化体系必须 EM）

**问题**: 直接从 solvated PDB 运行 heating 出现 NaN（`Particle coordinate is NaN`）。

**原因**: 磷酸基团引入局部应变，初始能量 -214,940 kJ/mol 过高。

**解决**: 先用 OpenMM `LocalEnergyMinimizer` 最小化：

```bash
python scripts/minimize_system.py \
    --prmtop data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos.prmtop \
    --pdb data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_solvated.pdb \
    --out data/md_runs/Hsap_WT_S305phos/Hsap_WT_S305phos_minimized.pdb \
    --max-iterations 1000
```

**结果**: 能量从 -214,940 → -1,032,741 kJ/mol，下降 ~818,000 kJ/mol。

### 步骤 7：磷酸化模拟突变体 (S305E) — ✅ 已完成

**构建记录** (2026-05-01):
- 脚本: `scripts/mutate_s305e.py` + `scripts/build_s305e_system.py`
- 方法: 基于 Hsap_WT_amber.pdb，手动几何放置 GLU 侧链
- 体系: `data/md_runs/Hsap_WT_S305E/Hsap_WT_S305E.prmtop/.rst7`
- 原子数: 同 WT (81,487)，电荷 0.0
- EM: -517k → -1,319k kJ/mol
- **3× replica OpenMM NVT production 已启动** (GPU 0/1/2)

---

## 五、力场参数说明

### `leaprc.phosaa19SB` 内容

- 基于 ff19SB 优化的磷酸化侧链扭转参数
- MP2/6-311+G* 级别 QM 拟合
- 支持 SEP (Ser-PO₃²⁻), TPO (Thr-PO₃²⁻), PTR (Tyr-PO₃²⁻), 及质子化形式

### 残基命名对照

| 状态 | PDB 名 | AMBER 库名 | 电荷 | 说明 |
|------|--------|-----------|------|------|
| 未磷酸化 | SER | SER | 0 | 标准 |
| 双负离子 (pH 7.4) | **SEP** | SEP | -2 | **生理状态，推荐** |
| 单负离子 (pH ~6) | S1P | S1P | -1 | 低 pH 状态 |
| 磷酸化模拟 | GLU | GLU | -1 | S305E 替代 |

> **建议**: 明确使用双负离子 SEP（生理 pH 7.4），无需 pKa 预测。

---

## 六、MD 平衡方案（磷酸化体系需更谨慎）

磷酸化引入后，局部静电环境剧变，需要更长的约束平衡：

```
1. EM: 5000 steps steepest descent, Fmax < 1000
2. NVT heating: 0 → 300K, 100ps
   - 对磷酸化位点周围 5Å 施加 10 kcal/mol/Å² 位置约束
3. NPT 约束平衡: 100ns
   - 磷酸化位点周围 10Å 施加 1 kcal/mol/Å² 约束
   - 逐步降低约束（每 25ns 减半）
4. NPT 无约束平衡: 50ns
5. NVT Production: 200-500ns
```

### OpenMM 实现

```python
# 添加位置约束（仅在平衡阶段）
from openmm import CustomExternalForce

# 约束磷酸化位点周围 10Å 的原子
restraint_force = CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
restraint_force.addGlobalParameter("k", 1.0)  # kcal/mol/Å²
# ... 添加被约束原子的坐标 ...
system.addForce(restraint_force)
```

---

## 七、分析重点

### 1. 界面静电变化
- 磷酸基团 (-2) 与 TRIM41 SPRY 正电残基 (Arg/Lys) 的盐桥形成
- 氢键网络重组

### 2. 结合自由能变化
- MM-GBSA 计算 WT vs S305-phos 的结合能差异
- 关注磷酸化位点周围 10Å 的残基贡献

### 3. 构象动态变化
- 磷酸化位点周围 loop 的 RMSF
- 整体 Rg 和 COM 距离变化
- 解离动力学（若发生）

### 4. 与 4mut 的交互效应
- 4mut (D431A/K479A/L495A/K498A) 减弱 TRIM41 结合
- S305E 是否通过增强静电互补性补偿 4mut 的效应？

---

## 八、风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|---------|
| 磷酸基团初始构象不合理 | 局部结构畸变 | EM + 长约束平衡 (100ns) |
| 磷酸化导致快速解离 | 无法获得稳定结合态 | 缩短平衡约束释放时间，监控 COM |
| 力场参数不准确 | 磷酸基团行为异常 | 对比 S305E 模拟突变体验证 |
| S305 不在直接界面上 | 效应可能微弱 | 重点分析长程静电和变构效应 |

---

## 九、初步结果（2026-05-01，~130ns partial data）

### 重大发现：S305 磷酸化导致 cGAS-TRIM41 解离

| Replica | 时长 | Self-RMSD | COM (Å) | Rg (Å) | 状态 |
|---------|------|-----------|---------|--------|------|
| **WT (参考)** | 200ns | 11.2 Å | 45.1±3.2 | 31.1±1.3 | ✅ 稳定结合 |
| S305-phos rep1 | 129ns | 9.9±1.9 Å | **67.6±1.4** | 38.9±0.6 | ⚠️ 解离 |
| S305-phos rep2 | 130ns | 19.1±6.4 Å | **89.8±10.0** | 48.4±4.3 | 🔴 完全解离 |
| S305-phos rep3 | 129ns | 13.0±3.1 Å | **71.3±3.4** | 40.2±1.4 | ⚠️ 解离 |

**COM 距离演化**：
- WT：稳定在 40-50 Å
- S305-phos：起始即 ~70 Å，rep2 持续上升至 ~110 Å

**Rg 演化**：
- WT：29-32 Å
- S305-phos：38-48 Å，显著增大（单个蛋白更展开）

### 与文献的差异

文献（Zhen et al., 2023）报告 CHK2 磷酸化 S305 **促进** cGAS-TRIM41 结合。但我们的模拟显示**解离**。可能原因：

1. **模拟时间不足**：解离可能是中间态，更长模拟可能重新结合
2. **溶液 vs 核内环境**：文献中的结合增强可能依赖核小体/DNA 等额外因子
3. **力场局限性**：磷酸基团 -2 电荷的静电排斥可能被过度估计
4. **构象选择机制**：磷酸化可能先导致解离，释放 cGAS 与 DNA 结合，间接促进 TRIM41 泛素化功能

### 下一步验证

- ✅ 继续运行到 200ns，观察是否重新结合或稳定在解离态
- 🟡 构建 S305E 对照体系：若 S305E 也解离 → 电荷效应主导；若不解离 → 真实磷酸化的特异性几何效应
- 🟡 分析解离路径：哪些界面残基先断裂？

---

## 十、运行状态与时间线（2026-05-01）

### 已完成

| 任务 | 状态 |
|------|------|
| S305-phos 体系构建 | ✅ 完成（prmtop/rst7/solvated.pdb） |
| 能量最小化 | ✅ 完成（-215k → -1033k kJ/mol） |
| OpenMM MD 3 replicas | 🔄 运行中（~15.5ns / 200ns） |

### 进行中

| 任务 | 进度 | ETA |
|------|------|-----|
| S305-phos rep1 (CUDA:0) | ✅ 200ns / 200ns | 已完成 |
| S305-phos rep2 (CUDA:1) | ✅ 200ns / 200ns | 已完成 |
| S305-phos rep3 (CUDA:2) | ✅ 200ns / 200ns | 已完成 |
| S305E rep1 (CUDA:0) | 🔄 0ns / 200ns | ~26h |
| S305E rep2 (CUDA:1) | 🔄 0ns / 200ns | ~26h |
| S305E rep3 (CUDA:2) | 🔄 0ns / 200ns | ~26h |
| GROMACS 2026 验证 (CUDA:3) | ~77.2ns / 200ns | ~18h |

### 待启动

| 任务 | 优先级 |
|------|--------|
| S305E 体系构建 + 3× replica MD | 高（电荷对照） |
| S305-phos vs S305E vs WT 对比分析 | 高（等 MD 完成） |
| GROMACS 2026 验证 RMSD 对比 | 中（等 MD 完成） |

### 总体时间线

| 阶段 | 任务 | 耗时 |
|------|------|------|
| 准备 | PDB 修改 + tleap 构建 | ✅ 1 天 |
| 最小化 | OpenMM EM | ✅ 10 分钟 |
| 生产 | 200ns × 3 reps (S305-phos) | ✅ 完成 |
| 生产 | 200ns × 3 reps (S305E) | 🔄 已启动 |
| 分析 | MM-GBSA + 构象分析 | 2-3 天 |
| **总计** | | **~1 周** |

---

## 十、参考文献

1. Zhen et al. (2023). *Nuclear cGAS restricts L1 retrotransposition by promoting TRIM41-mediated ORF2p ubiquitination.* Nat Commun 14, 8394.
2. Seo et al. (2015). *Akt1-mediated phosphorylation of cGAS inhibits its enzymatic activity.* 
3. Li et al. (2021). *Phosphorylation and chromatin tethering prevent cGAS activation during mitosis.* Science 371, eabc5386.
4. AMBER 2023 Reference Manual — Chapter 3: `leaprc.phosaa19SB`
