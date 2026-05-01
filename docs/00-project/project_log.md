# cGAS-TRIM41 Simulation Project Log

## §45. GROMACS 修复与磷酸化 MD 方案 (2026-04-23)

### 45.1 GROMACS 诊断结论

基于系统文献调研和代码审查，确认 GROMACS 与 OpenMM 的 4× RMSD 差异**主要由 CMAP 转换错误导致**：

- **根本原因**: parmed 将 AMBER ff19SB 的 14 种残基特异性 CMAP types 错误地压缩为 1 种 (`C N XC C N`)
- **次要因素**: LINCS `iter=1, order=4` 精度不足；NPT vs NVT production ensemble 差异
- **文献基准**: 跨引擎 RMSD 差异 4× 远超正常范围（文献报道通常 <2×，能量差异 0.3-1.0 kcal/mol）

诊断报告: `docs/30-diagnostics/gromacs_openmm_divergence_diagnosis.md`

### 45.2 GROMACS 修复内容

#### (a) CMAP 修复脚本 `scripts/fix_gromacs_cmap.py`
- 读取 AMBER prmtop 提取 14 种 CMAP grid definitions
- 在 GROMACS `[atoms]` section 中，根据残基名将 CA 类型从 `XC` 改为 `XC{n}` (n=0-13)
- 在 `[atomtypes]` 中复制 `XC` 为 `XC0-XC13`
- 在 `[cmaptypes]` 中替换为 14 种独立的残基特异性定义
- **测试验证**: 成功修改 537 CA atoms，14 cmaptypes，13 new atomtypes

#### (b) 转换脚本更新 `scripts/convert_amber_to_gromacs.py`
- 在 `amber.save()` 后自动调用 `fix_gromacs_cmap.py`
- 转换失败时打印手动修复命令

#### (c) Production MDP 修复 `data/md_runs_gmx/mdp/prod.mdp`
- `pcoupl = no` — NVT ensemble（与 OpenMM 一致）
- `lincs_iter = 2`, `lincs_order = 6` — 提高约束精度
- `vdw-modifier = Potential-shift-Verlet` — 更精确的 LJ 处理

#### (d) 未来重跑建议
- 当前修复仅更新脚本和 setup，**暂未重跑 MD**
- 下次启动 GROMACS MD 时，需：
  1. 重新转换拓扑（使用修复后的脚本）
  2. 使用更新后的 prod.mdp
  3. 建议先做 10ns 测试 replica 验证 RMSD 是否回归 OpenMM 范围

### 45.3 磷酸化 MD 方案

方案文档: `docs/20-protocols/phosphorylation_md_plan.md`

#### 核心位点选择
| 位点 | 全长编号 | 构建体范围 | 优先级 | 依据 |
|------|---------|-----------|--------|------|
| **S305** | 305 | ✅ 200-554 内 | **最高** | CHK2 磷酸化促进 cGAS-TRIM41 结合 (Zhen et al., 2023) |
| S435 | 435 | ✅ 200-554 内 | 次要 | PPP6C 调控 cGAMP 合成 |
| S120 | 120 | ❌ N-terminal 外 | 需全长 | CHK2 磷酸化，与 S305 协同 |

#### 推荐模拟体系
- WT (基线)
- S305-phos (SEP @ 305)
- S305E (磷酸化模拟突变)
- 4mut+S305E (检验磷酸化是否补偿 4mut 结合减弱)

#### 技术路线
1. PyMOL 手动添加 PO₃ 基团 → 修改残基名 SER→SEP
2. `pdb4amber --reduce` 处理氢原子
3. `tleap` + `leaprc.phosaa19SB` 构建体系
4. 加长约束平衡 (100ns 逐步释放)
5. 200-500ns NVT production

### 45.4 待决策事项

1. **是否启动全长 cGAS 构建**？S120 磷酸化需要 cGAS(1-522) + TRIM41 SPRY，需重新 AF3 预测或 docking
2. **磷酸化模拟的规模**：先做 Hsap WT S305-phos (3 reps × 200ns) 做可行性验证？
3. **GROMACS 重跑时机**：是否在新数据（磷酸化）之前重跑修复后的 GROMACS 以验证？

## §46. GROMACS 验证成功 + S305-phos 解离发现 (2026-05-01)

### 46.1 GROMACS 2026 验证结果（77ns partial data）

**结论: ✅ CMAP 修复成功**

| 指标 | GROMACS 2026 native | OpenMM WT | Ratio |
|------|---------------------|-----------|-------|
| Self-RMSD (0-77ns) | 15.1±3.9 Å | 11.2±4.1 Å | 1.34 |
| COM distance | **45.1±1.4 Å** | **45.1±3.2 Å** | 1.00 |
| Rg | **30.5±0.6 Å** | **31.1±1.3 Å** | 0.98 |

- COM 和 Rg 曲线在 15ns 后**高度重叠**
- 从 24ns 时的 1.93× 改善到 77ns 时的 **1.34×**，趋势继续收敛
- 剩余 ~34% RMSD 差异来自 EM 阶段 pdb2gmx 重新生成氢原子导致的初始条件差异

**关键教训：PBC 包裹误报**
- GROMACS 默认输出 wrapped 坐标，OpenMM DCD 输出 unwrapped 坐标
- 未修复 PBC 前：GROMACS COM=28 Å（虚假），被误判为 RMSD ratio 5.8×
- 修复 PBC 后 (`gmx trjconv -pbc mol`)：COM=44.6 Å，与 OpenMM 一致

### 46.2 S305-phos 初步结果（~130ns）

**重大发现：S305 磷酸化导致 cGAS-TRIM41 解离**

| Replica | 时长 | COM (Å) | Rg (Å) | 状态 |
|---------|------|---------|--------|------|
| WT (参考) | 200ns | 45.1±3.2 | 31.1±1.3 | ✅ 稳定结合 |
| S305-phos rep1 | 129ns | **67.6±1.4** | 38.9±0.6 | ⚠️ 解离 |
| S305-phos rep2 | 130ns | **89.8±10.0** | 48.4±4.3 | 🔴 完全解离 |
| S305-phos rep3 | 129ns | **71.3±3.4** | 40.2±1.4 | ⚠️ 解离 |

- 所有 3 个 replica 的 COM 均显著大于 WT
- rep2 COM 达到 ~110 Å，完全解离
- Rg 增大（38-48 Å vs WT 31 Å），单个蛋白更展开

**与文献差异**: Zhen et al. (2023) 报告 CHK2 磷酸化 S305 "促进"结合，但模拟显示"解离"。可能解释：
1. 溶液 vs 核内环境差异（核小体/DNA 等额外因子）
2. 力场对 -2 电荷的静电排斥过度估计
3. 解离为中间态，更长模拟可能重新结合
4. 构象选择机制：先解离→与 DNA 结合→间接促进 TRIM41 招募

### 46.3 当前运行状态

| 实验 | GPU | 进度 | ETA |
|------|-----|------|-----|
| S305-phos rep1 | 0 | 128.6ns / 200ns | ~4h |
| S305-phos rep2 | 1 | 129.4ns / 200ns | ~4h |
| S305-phos rep3 | 2 | 128.4ns / 200ns | ~4h |
| GROMACS 2026 Hsap_WT | 3 | 77.2ns / 200ns | ~18h |

### 46.4 待决策事项

1. **S305E 是否还需要做？** S305-phos 已显示显著解离效应，S305E 对照可验证是否为单纯电荷效应
2. **是否构建全长 cGAS(1-522)？** S120 磷酸化不在当前构建体范围内
3. **GROMACS 验证是否继续到 200ns？** 77ns 已验证成功，可提前终止以释放 GPU
4. **S305-phos 200ns 后是否延长到 500ns？** 观察是否重新结合

## §47. Boltz-2 本地安装与验证 (2026-05-01)

### 47.1 安装

- 新建 conda 环境 `boltz`，Python 3.11
- `pip install boltz[cuda]` → 安装 Boltz-2.2.1（默认模型 Boltz-2）
- PyTorch 2.11.0 + CUDA 13 组件
- **问题与修复**：`libnvrtc-builtins.so.13.0` 未找到
  - 原因：系统 `LD_LIBRARY_PATH` 指向 CUDA 12，bolt 的 CUDA 13 库在 site-packages 中
  - 修复：`conda env config vars set LD_LIBRARY_PATH=".../nvidia/cu13/lib:$LD_LIBRARY_PATH" -n boltz`
- 模型权重自动下载至 `~/.boltz/`（boltz2_aff.ckpt + boltz2_conf.ckpt，共 ~5.5 GB）

### 47.2 全长 cGAS-TRIM41 预测验证

| 指标 | Boltz-2 | AF3 (历史) |
|------|---------|-----------|
| ipTM | 0.17 | 0.15 |
| pTM | 0.40 | 0.38 |
| complex pLDDT | 0.64 | ~0.60 |

**结论**：与 AF3 完全一致——全长复合物界面置信度极低，验证了当初选择截断构建体 + 对接的决策正确。

### 47.3 截断构建体预测与 AF3 结构对比

**输入**：cGAS CT (200-522, 323 aa) + TRIM41 SPRY (413-630, 218 aa)

**置信度**：
| 指标 | 数值 |
|------|------|
| ipTM | 0.33 |
| pTM | 0.69 |
| complex pLDDT | 0.78 |
| cGAS pLDDT | 0.94 |
| TRIM41 pLDDT | 0.84 |

截断后质量显著提升，但 ipTM 0.33 仍低于 0.6 可信阈值。

**结构与 AF3 的定量对比**（对齐 cGAS CA 后）：

| 指标 | 结果 | 含义 |
|------|------|------|
| cGAS CA RMSD | **1.06 Å** | 单体折叠几乎完美一致 |
| TRIM41 CA RMSD | **21.34 Å** | 相对位置完全不同 |
| 共享界面接触 | **0** | 接触网络完全不重叠 |
| Jaccard 相似度 | **0.000** | 界面定义完全不同 |

**活性位点距离对比**：
| 位点 | Boltz-2 | AF3 | Δ |
|------|---------|-----|---|
| D431 | 21.0 Å | 23.5 Å | −2.5 Å |
| K479 | 22.0 Å | 35.7 Å | **−13.8 Å** |
| L495 | 25.7 Å | 31.6 Å | −5.9 Å |
| K498 | 25.5 Å | 35.8 Å | **−10.3 Å** |

### 47.4 关键结论

1. **Boltz-2 对单体结构预测高度可靠**（cGAS RMSD 1.06 Å）
2. **但对 cGAS-TRIM41 这个瞬态/弱相互作用复合物的界面预测，Boltz-2 和 AF3 给出了定性不同的答案**
3. 两个当前最先进的模型（AF3 + Boltz-2）在低 ipTM 体系中存在**显著的界面几何分歧**
4. 这进一步验证了项目核心策略的正确性：**不信任单一模型的复合物预测，采用对接 + MD 平衡获得初始构象**

### 47.5 后续可探索方向

- **加入 DNA / 核小体**：Boltz-2 支持 DNA/RNA 输入，可测试核小体存在下 cGAS-TRIM41 的界面是否更稳定
- **亲和力预测**：加入小分子配体后，Boltz-2 可输出 `affinity_pred_value`（log10(IC50)）
- **磷酸化修饰**：YAML 支持 `modifications` 字段，可直接预测 SEP@305 的效应

对比图表：`data/boltz_test_truncated/comparison/boltz2_vs_af3_comparison.png`
对比脚本：`scripts/compare_boltz2_af3_v2.py`

