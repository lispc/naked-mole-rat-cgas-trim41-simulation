# cGAS-TRIM41 MD 研究项目日志

> 本文档为项目主日志索引，详细记录见各子文档。
> 
> 最后更新：2026-04-23

---

## 目录

| 文档 | 内容 | 状态 |
|------|------|------|
| [docking_report.md](docking_report.md) | 蛋白-蛋白对接完整实验记录 | ✅ 已完成 |
| [af3_report.md](af3_report.md) | AF3 结构预测结果与序列分析 | ✅ 已完成 |
| [hardware_benchmark.md](hardware_benchmark.md) | 硬件资源与性能基准 | ✅ 已完成 |
| [software_versions.md](software_versions.md) | 软件版本记录 | ✅ 已完成 |
| [execution_plan_v1.md](execution_plan_v1.md) | 执行方案 v1.0 | ✅ 已完成 |
| [computational_workflow.md](computational_workflow.md) | 原始方案设计 | ✅ 已完成 |
| [paper_notes_cgas_trim41.md](paper_notes_cgas_trim41.md) | 论文理解笔记 | ✅ 已完成 |
| [cluspro_submission_guide.md](cluspro_submission_guide.md) | ClusPro 提交指南 | ✅ 已完成 |

---

## 一、项目概况

基于 Chen et al., Science 2025，研究裸鼹鼠 cGAS 的 4 个氨基酸变异（C463S, K479E, L495Y, K498T，即人源→裸鼹鼠）如何通过改变 cGAS-TRIM41 相互作用影响泛素化效率。

**论文定位**：作为实验论文的 In silico validation（Cell Reports / Nature Communications 级别）

**计算目标**：通过分子动力学模拟和计算突变扫描，在原子水平解释 4 个变异如何影响 cGAS-TRIM41 相互作用。

---

## 二、关键决策摘要

| 决策 | 选择 | 理由 |
|------|------|------|
| 结构预测 | AF3 Server | 最准确，免费；本地运行对 36GB 内存压力大 |
| 结构域策略 | 截取结构域 | cGAS CTD (200-554) + TRIM41 SPRY (413-630)；原子数减半，可做 200ns×3 重复 |
| 突变扫描 | Rosetta | 完全免费，不需要 FoldX license |
| GPU 后端 | OpenCL | Metal 在 conda-forge 上有兼容性问题 |
| 蛋白对接 | LightDock | ClusPro/SDOCK2.0 失败后，LightDock 在 py311 conda env 中成功运行 |

---

## 三、突变映射

| 论文标记 | 裸鼹鼠 aa | 人源对应位置 | 人源 aa | 突变方向 |
|----------|-----------|-------------|---------|----------|
| S463 | S | **463** | C | C→S |
| E511 | E | **479** | K | K→E |
| Y527 | Y | **495** | L | L→Y |
| T530 | T | **498** | K | K→T |

论文中的编号基于裸鼹鼠 cGAS 坐标系统（554 aa），人源为 522 aa。

---

## 四、当前状态总览

### ✅ 已完成

1. **环境搭建**：`cgas-md` conda env (Python 3.13) + `py311` conda env (LightDock)
2. **AF3 结构预测**：4 个 job 全部完成（ipTM 均 <0.25，单体可接受）
3. **域截断**：cGAS CTD (200-554) + TRIM41 SPRY (413-630)，4 个 job 全部处理
4. **蛋白-蛋白对接**：LightDock 结构域对接全部完成
   - **Hgal**：20/20 poses 成功，最佳 pose 已保存 (`best_pose.pdb`)
   - **Hsap (无约束)**：0/25 成功，495/498 始终在界面外
   - **Hsap (restraints, 50sw/200step)**：0/25 成功，物理上不可能
5. **核心发现**：裸鼹鼠 4 个突变将活性位点从分散的 ~28.6Å 聚集成紧凑的 ~18.4Å 补丁
6. **发表级图表**：PyMOL + matplotlib 生成 4 张图，2400×1800 dpi=300

### ⏳ 运行中

无

### 📋 待办（按优先级）

- [ ] **MD 体系构建**（Hgal best_pose → Amber tleap → OpenMM）← **当前最高优先级**
- [ ] 生产 MD 模拟（200ns × 3 重复）
- [ ] Hsap 4mut 结构获取（AF3 重新提交 或 PyMOL in-silico 突变）
- [ ] MM-GBSA 结合能计算
- [ ] Rosetta 突变扫描
- [ ] 论文 Methods 撰写
- [ ] 全长对接（可选，评估是否需要）

---

## 五、文件速查

```
sequences/
  Hsap_cGAS_WT.fasta, Hsap_cGAS_4mut.fasta
  Hgal_cGAS_WT.fasta, Hgal_cGAS_4mut_rev.fasta
  TRIM41_WT.fasta

structures/af3_raw/
  job1_Hsap_WT/     # cgas_fixed.pdb, trim41_fixed.pdb, cgas_CT_200-554.pdb, trim41_SPRY_413-630.pdb
  job2_Hsap_4mut/   # (同上)
  job3_Hgal_WT/     # (同上)
  job4_Hgal_4mut_rev/  # (同上)

structures/docking/lightdock/
  Hgal_domain/best_pose.pdb     # ← MD 起始结构
  Hsap_domain/
  Hsap_restrained/

figures/                          # 发表级图表
  hgal_active_residues.png
  hsap_active_residues.png
  comparison_overlay.png
  distance_comparison.png
  distance_data.txt

scripts/
  process_af3_results.py
  analyze_lightdock.py
  calc_residue_distances.py       # PDB 精确距离计算
  quick_analyze_top_poses.py      # LightDock top poses 分析
  visualize_active_residues.pml   # PyMOL 自动化作图
  build_system.py, run_md.py, analyze_trajectory.py, run_mmpbsa.py

docs/
  project_log.md                # 本文档（主索引）
  docking_report.md             # 对接完整记录
  af3_report.md                 # AF3 结果
  hardware_benchmark.md         # 硬件基准
  software_versions.md          # 软件版本
  execution_plan_v1.md
  computational_workflow.md
  paper_notes_cgas_trim41.md
  cluspro_submission_guide.md
```

---

## 六、时间线

```
2026-04-21  环境搭建（Miniforge + conda env + OpenMM）
2026-04-22  AF3 提交（4 jobs），AF3 结果处理，ClusPro 提交
2026-04-23  ClusPro 结果分析（失败），SDOCK2.0 尝试（放弃）
            LightDock 安装（py311 env），Hsap/Hgal 结构域对接
            核心发现：空间几何差异（18Å vs 28Å）
2026-04-23  Hsap 增强对接启动（restraints, 50sw/200step）
            Hsap 增强对接完成（0/25 失败，确认几何约束）
            精确距离测量（PDB CA-CA: 18.43Å vs 28.63Å）
            PyMOL + matplotlib 图表生成（4 张发表级图）
            README.md + 文档全面更新
```

---

*文档创建：2026-04-22*
*最后更新：2026-04-23 ( docking 阶段完成，进入 MD 准备 )*
*维护者：Kimi Code CLI*
