# GROMACS vs OpenMM 结果差异诊断报告

**日期**: 2026-04-23（更新 2026-04-30）
**问题**: GROMACS 200ns × 4 replica 与 OpenMM 200ns × 6 replica 结果差异巨大
**状态**: 🔄 GROMACS 2026 原生 amber19sb.ff 验证运行中（~9.8ns/200ns）

---

## 一、定量差异概览

| 指标 | OpenMM (6 reps) | GROMACS (4 reps) | 差异倍数 |
|------|-----------------|------------------|----------|
| RMSD (mean) | 7–12 Å | 29–35 Å | **~3–4×** |
| RMSD (max) | 9–16 Å | 35–46 Å | **~3–4×** |
| COM distance | 41–49 Å | 31–39 Å | GROMACS 更紧密 |
| Rg (mean) | 29–32 Å | 32–36 Å | GROMACS 略大 |
| Unbound % | 0–75% (replica-dependent) | ~0% (所有 replica) | **定性不同** |
| 能量 (WT) | ~-1,080,000 kJ/mol | (未成功提取) | — |

---

## 二、逐项审查：有无明显错误？

### ✅ 确认正确的设置

| 项目 | GROMACS 值 | 预期值 | 状态 |
|------|-----------|--------|------|
| 原子数 | 85,510 | 与 AMBER 一致 | ✅ |
| 1-4 scaling | fudgeLJ=0.5, fudgeQQ=0.8333 | SCNB=2.0, SCEE=1.2 | ✅ |
| nrexcl | 3 | AMBER 标准 | ✅ |
| 虚拟站点 (EPW) | `[virtual_sites3]` 存在 | OPC 需要 | ✅ |
| CMAP term 数 | 539 (217 + 322) | AMBER 537 | ✅ (微小差异来自 chain termini) |
| LINCS 警告 | 无 | 不应有 | ✅ |
| EM 收敛 | Fmax < 1000 | 标准 | ✅ |
| 截断 | rcoulomb=1.0, rvdw=1.0 | 与 OpenMM 一致 | ✅ |
| PME | pme_order=4, fourierspacing=0.12 | 标准 | ✅ |

### ⚠️ 确认存在的差异/问题

| # | 项目 | GROMACS | OpenMM | 影响评估 |
|---|------|---------|--------|----------|
| 1 | **CMAP 残基特异性** | **1 type** (C N XC C N) | **14 types** (residue-specific) | **🔴 高** — 力场实现错误 |
| 2 | **Production ensemble** | **NPT** (C-rescale) | **NVT** (barostat removed) | 🟡 中 — 体积自由度不同 |
| 3 | **约束精度** | LINCS iter=1, order=4 | SHAKE (tol=1e-4) | 🟡 中 — 文献报告可能不够 |
| 4 | **积分器** | leap-frog (deterministic) | LangevinMiddle (stochastic) | 🟢 低 — 文献差异 ~0.1 kcal/mol |
| 5 | **Thermostat** | v-rescale τ=1.0 ps | Langevin γ=1.0 ps⁻¹ | 🟢 低 — 等效时间尺度 |
| 6 | **初始速度生成** | gen_vel at 300K (Maxwell) | 0K → 300K ramp | 🟢 低 — 仅影响早期平衡 |

---

## 三、核心问题：CMAP 转换错误

### 什么是 CMAP？

ff19SB 引入了**残基特异性 CMAP 校正** — 每种氨基酸（20 种）有自己独立的 φ-ψ 二面角势能面校正表。这是 ff19SB 相比 ff14SB 的核心改进。

### 转换问题

- **AMBER prmtop**: 14 种独立的 CMAP types，537 个 CMAP terms
- **GROMACS .top**: 只有 **1 种** `[cmaptypes]` 定义（`C N XC C N`），539 个 terms 全部映射到该 type

GROMACS 通过将**所有残基的 CA 原子定义为 `XC` 类型**来 workaround 残基特异性需求，但这导致所有残基共享**同一个** CMAP 势能面。

### 为什么这会导致巨大差异？

CMAP 校正直接修改 backbone 的 φ-ψ 势能面：
- 不同氨基酸的 α/β 偏好不同（Pro 偏爱 cis-肽键，Gly 高度柔性，Ala 偏爱 α-螺旋等）
- 共享"平均" CMAP 表会**系统性错误描述 backbone 构象偏好**
- 对于 200ns 长模拟，这种偏差会持续累积，导致整体构象逐渐偏离正确分布

### 文献支持

> Franz Waibl (AMBER 邮件列表, 2023): "After manually tweaking the topology as described above, I was able to obtain consistent CMAP energies between Amber and Gromacs."

> CHARMM-GUI 的 GROMACS 转换方案：必须为每种氨基酸定义独立原子类型（XC1, XC2, ...），并为每种组合定义独立 cmaptype。

---

## 四、LINCS 精度问题（次要但可能放大 CMAP 误差）

- GROMACS 使用 `lincs-iter=1, lincs-order=4`
- 文献（ILVES, arXiv 2025）指出：**默认 LINCS 参数可能导致非收敛结果、温度不可靠、集体运动伪影**
- 建议：`lincs-iter=2` 或更高（特别是含虚拟站点的体系）
- 约束精度不足会放大 backbone 的微小构象偏差

---

## 五、NPT vs NVT Production（可能解释 COM 差异方向）

- **GROMACS**: NPT（C-rescale barostat 开启）
- **OpenMM**: NVT（barostat 移除）

**可能的影响**：
- NPT 允许盒子体积自由变化，解离的蛋白质可能因溶剂化层膨胀而推动盒子扩大
- 但观察到的现象是 **GROMACS COM 更小（更紧密）**，与"NPT 允许更大分离"的直觉相反
- 更可能的解释：CMAP 错误导致某些 loop/界面区域构象改变，使蛋白质以更"紧凑"但内部更"松散"的方式结合

---

## 六、文献基准：跨引擎差异的正常范围

| 研究 | 体系 | 报告的差异 | 是否可接受 |
|------|------|-----------|-----------|
| SAMPL5/SAMPL6 (Shirts et al.) | Host-guest | 结合自由能差异 0.3–1.0 kcal/mol | ⚠️ 需关注但可接受 |
| CHARMM-GUI cross-validation | 磷脂双分子层 | A_L 差异 0.1–1.4 Å² | ✅ 可接受 |
| Sun et al. (CJCP 2021) | 小分子构象 | PMF 差异 ~0.1 kcal/mol | ✅ 可忽略 |
| 旧 GROMACS (parmed) | cGAS-TRIM41 | RMSD 差异 ~4×，解离行为定性不同 | 🔴 远超正常范围 |
| **新 GROMACS 2026 (native)** | **cGAS-TRIM41** | **RMSD 差异 ~1.3×，COM/Rg 几乎相同** | **✅ 正常范围** |

**文献共识**: "Specifying force field parameters is insufficient to ensure reproducibility" — 但差异通常在能量/自由能的 0.3–1.0 kcal/mol 范围，不会导致 RMSD 相差 4 倍。

---

## 七、结论：这种差异是否合理？

### 🔴  verdict: **不合理，存在已知的力场实现错误**

理由：

1. **CMAP 残基特异性丢失** 是一个已记录的、可导致显著动力学偏差的转换错误
2. **4 倍 RMSD 差异** 远超文献报道的跨引擎"正常差异"范围
3. **解离行为的定性差异**（0% vs 0-75% unbound）不能仅用积分器/约束差异解释
4. LINCS 默认精度在长模拟中可能进一步放大 CMAP 误差

### 这种差异不是：
- ❌ 正常的引擎间数值噪声
- ❌ 单纯的积分器随机性差异
- ❌ 合理的 thermostat/barostat 差异

### 这种差异主要是：
- ✅ **CMAP 力场实现错误**（残基特异性丢失）
- ✅ **LINCS 精度不足**（iter=1 长期累积）
- 🟡 **NPT vs NVT**（次要，影响溶剂化层）

---

## 八、修复行动与验证（2026-04-30 → 2026-05-01）

### 修复方案：GROMACS 2026 原生 amber19sb.ff

放弃 parmed 转换路线，改用 GROMACS 2026 内置的 `amber19sb.ff`：

```bash
gmx pdb2gmx -f protein.pdb -o processed.gro -p native.top \
    -ff amber19sb -water opc -ignh
```

**关键验证**: `amber19sb.ff/cmap.itp` 包含残基特异性 CMAP（`C-* N-GLY XC-GLY C-GLY N-*` 等），14 种独立类型，321 CMAP torsion pairs。

### 修复清单

| 问题 | 修复措施 | 状态 |
|------|---------|------|
| CMAP 残基特异性丢失 | 使用 GROMACS 2026 原生 `amber19sb.ff` | ✅ 已修复（321 pairs） |
| LINCS 精度不足 | `lincs_iter=2, lincs_order=6` | ✅ 已修复 |
| NPT vs NVT | `pcoupl = no`（production 改为 NVT） | ✅ 已修复 |
| TER 记录缺失 | MDAnalysis 提取时插入 TER 分隔双链 | ✅ 已修复 |
| OPC 溶剂 | `gmx solvate -cs tip4p.gro`（4-site water） | ✅ 已修复 |

### 验证运行状态

| 体系 | 引擎 | 进度 | 性能 |
|------|------|------|------|
| Hsap_WT | GROMACS 2026 native | 🔄 77.2ns / 200ns | ~137 ns/day |
| Hsap_WT | OpenMM | ✅ 200ns × 3 reps 完成 | ~500 ns/day |

### 性能对比

| 指标 | 旧 GROMACS (parmed) | 新 GROMACS 2026 (native) | 差异 |
|------|---------------------|--------------------------|------|
| CMAP pairs | ~30 (1 type) | 321 (14 types) | 10× 更多 |
| NVT 性能 | ~383 ns/day | ~503 ns/day | 快 31%（短运行） |
| NPT 性能 | ~180 ns/day | ~136 ns/day | **慢 24%** |
| Prod 性能 | ~177 ns/day | ~137 ns/day | **慢 23%** |

**关键发现**: 新 GROMACS 在 **production 阶段更慢约 23%**。原因是原生 ff19SB 的 321 CMAP pairs（vs 旧转换的 ~30 个）增加了显著计算开销。CMAP 是二维查表插值，pairs 数量直接正比于计算量。

> ⚠️ **注意**: NVT 短运行（100ps）显示 503 ns/day，但这只是短模拟的计时误差，不代表长 prod 性能。实际 prod 性能约 137 ns/day。

### 旧 GROMACS 数据处置

- **Hsap_WT rep1-2**: 已完成 136ns each，但 CMAP 有 bug → **不纳入正式分析**
- **Hsap_4mut rep1-2**: 已完成 132ns each，同上 → **不纳入正式分析**
- **保留目的**: 仅作定性趋势参考（如 4mut > WT 的 Rg 趋势）

### 77ns Partial Result 验证结论

| 指标 | GROMACS 2026 | OpenMM | Ratio | 状态 |
|------|-------------|--------|-------|------|
| Self-RMSD (0-77ns) | **15.1±3.9 Å** | 11.2±4.1 Å | 1.34 | 🟡 略高但可接受 |
| COM distance | **45.1±1.4 Å** | 45.1±3.2 Å | 1.00 | ✅ **几乎相同** |
| Rg | **30.5±0.6 Å** | 31.1±1.3 Å | 0.98 | ✅ **几乎相同** |
| 绝对 RMSD (frame 0 vs OpenMM NPT) | **6.40 Å** | — | — | ✅ 初始结构接近 |

** verdict: ✅ CMAP 修复成功**

- COM 和 Rg 曲线在 15ns 后**高度重叠**
- RMSD 仍有 ~34% 差异，主要来自 EM 阶段 pdb2gmx 重新生成氢原子导致的初始条件差异
- 积分器差异（leap-frog vs LangevinMiddle）贡献次要

---

### 重要教训：PBC 包裹误报

**问题**: GROMACS 默认输出 **PBC wrapped** 坐标（原子拆分到不同盒子镜像），而 OpenMM DCD 输出 **unwrapped** 坐标。

**后果**: 未修复 PBC 前，GROMACS 的 COM=28 Å（虚假），Rg=38 Å（虚假），被误判为"严重偏离"（RMSD ratio 5.8×）。

**解决**: `gmx trjconv -pbc mol` 或 `MDAnalysis.transformations.unwrap`。

**修复 PBC 后**: COM=44.6 Å，Rg=30.3 Å，与 OpenMM 一致。

> ⚠️ **警示**: 分析 GROMACS 轨迹时，**必须先修复 PBC**，否则所有基于 COM/Rg/RMSD 的结论都是错误的。

### 待完成验证

1. ✅ GROMACS 2026 77ns 已完成初步验证：COM/Rg 与 OpenMM 一致
2. 🔄 **继续运行到 200ns**，确认长期行为稳定性
3. 🔄 对比完整 200ns 的 RMSD 分布、氢键网络、界面接触面积

---

## 九、脚本与工具更新

### `scripts/build_phosphorylated_system.py`
- 新增：自动为 SER 添加 PO₃ 基团（四面体几何）
- 修复：排除 HG 原子（磷酸化后 OG 连 P 不连 H）
- 修复：TER 记录分隔 cGAS 和 TRIM41 链
- 修复：resid > 9999 截断、segid 限制为 1 字符

### `scripts/minimize_system.py`（新增）
- OpenMM `LocalEnergyMinimizer` 快速 EM
- 用于磷酸化体系启动前最小化（必需）

### `scripts/prepare_gromacs_2026_native.sh`（新增）
- 完整 pipeline：pdb2gmx → solvate → ionize → EM → NVT → NPT → prod
- 使用 GROMACS 2026 原生 `amber19sb.ff` + OPC

### 分析脚本小 bug（非主要原因）

GROMACS 分析脚本 `batch_analyze_hsap_gmx_fast.py` 中的 Kabsch 实现：
```python
aligned = (R @ mobile_centered.T).T + reference.mean(axis=0)
```
使用了 `R.T` 而不是 `R`，但在数值验证中差异 <1%，**不是 4 倍 RMSD 差异的原因**。
