# cGAS-TRIM41 Simulation Project Log

## §45. GROMACS 修复与磷酸化 MD 方案 (2026-04-23)

### 45.1 GROMACS 诊断结论

基于系统文献调研和代码审查，确认 GROMACS 与 OpenMM 的 4× RMSD 差异**主要由 CMAP 转换错误导致**：

- **根本原因**: parmed 将 AMBER ff19SB 的 14 种残基特异性 CMAP types 错误地压缩为 1 种 (`C N XC C N`)
- **次要因素**: LINCS `iter=1, order=4` 精度不足；NPT vs NVT production ensemble 差异
- **文献基准**: 跨引擎 RMSD 差异 4× 远超正常范围（文献报道通常 <2×，能量差异 0.3-1.0 kcal/mol）

诊断报告: `docs/gromacs_openmm_divergence_diagnosis.md`

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

方案文档: `docs/phosphorylation_md_plan.md`

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
