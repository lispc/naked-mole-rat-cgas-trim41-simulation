# 文档索引

> 最后更新：2026-05-06

---

## 00-project — 项目元信息

| 文件 | 内容 |
|------|------|
| [`project_log.md`](00-project/project_log.md) | **主项目日志**（§45 起，持续更新） |
| [`project_log_2026_04.md`](00-project/project_log_2026_04.md) | 2026-04 历史日志（已归档） |
| [`summary_2026-04-23.md`](00-project/summary_2026-04-23.md) | 2026-04-23 阶段性总结（快照） |
| [`nice_to_have.md`](00-project/nice_to_have.md) | 未来计划 + 想法 backlog（合并自 good_ideas_backlog） |
| [`proposed_validation_experiments.md`](00-project/proposed_validation_experiments.md) | 变构假说的计算验证方案 |
| [`paper_notes_cgas_trim41.md`](00-project/paper_notes_cgas_trim41.md) | ⚠️ 早期论文笔记（含过时信息，见文件头警告） |
| [`4mut_correction_log.md`](00-project/4mut_correction_log.md) | 4mut 位点 C463S→D431S 修正记录 |

## 10-reports — 实验报告

| 文件 | 内容 |
|------|------|
| [`docking_report.md`](10-reports/docking_report.md) | 蛋白-蛋白对接完整记录（含 ClusPro/SDOCK2.0/LightDock/Rosetta） |
| [`af3_report.md`](10-reports/af3_report.md) | AlphaFold3 结构预测结果 |
| [`af3_mutation_analysis.md`](10-reports/af3_mutation_analysis.md) | AF3 突变体 + 嵌合体分析 |
| [`corrected_4mut_analysis.md`](10-reports/corrected_4mut_analysis.md) | 修正后 4mut 变构效应分析 |
| [`interface_analysis_report.md`](10-reports/interface_analysis_report.md) | 界面位置分析（N-端 vs C-端） |
| [`allosteric_mechanism_analysis.md`](10-reports/allosteric_mechanism_analysis.md) | 4mut 变构机制：三种功能假说 + 实验计划 |
| [`rosetta_docking_report.md`](10-reports/rosetta_docking_report.md) | Rosetta 对接报告 |
| [`rosetta_docking_4mut_report.md`](10-reports/rosetta_docking_4mut_report.md) | Rosetta 4mut 对接报告 |
| [`rosetta_mutational_scan_report.md`](10-reports/rosetta_mutational_scan_report.md) | Rosetta 突变扫描 |
| [`cite28_analysis.md`](10-reports/cite28_analysis.md) | 引用28（Zhen et al.）分析 |
| [`cluspro_submission_guide.md`](10-reports/cluspro_submission_guide.md) | ClusPro 提交指南 |

## 15-literature — 文献调研

| 文件 | 内容 |
|------|------|
| [`computational_ubiquitination_research_survey.md`](15-literature/computational_ubiquitination_research_survey.md) | E3 泛素化计算方法综述（~20,000 字） |

## 20-protocols — 方案与执行计划

| 文件 | 内容 |
|------|------|
| [`computational_workflow.md`](20-protocols/computational_workflow.md) | 计算实验总体设计方案 |
| [`execution_plan_v1.md`](20-protocols/execution_plan_v1.md) | 执行计划 v1.0 |
| [`phosphorylation_md_plan.md`](20-protocols/phosphorylation_md_plan.md) | 磷酸化 MD 方案与运行记录 |
| [`quaternary_complex_experimental_design.md`](20-protocols/quaternary_complex_experimental_design.md) | 四元复合物详细实验设计 |

## 30-diagnostics — 问题诊断

| 文件 | 内容 |
|------|------|
| [`gromacs_openmm_divergence_diagnosis.md`](30-diagnostics/gromacs_openmm_divergence_diagnosis.md) | GROMACS vs OpenMM 差异诊断（CMAP bug 发现与修复） |
| [`md_realism_gap.md`](30-diagnostics/md_realism_gap.md) | MD 模拟与真实生物过程的差距分析 |
| [`quaternary_mvp_structure_assessment.md`](30-diagnostics/quaternary_mvp_structure_assessment.md) | 四元 MVP 结构问题评估与重建计划 |

## 40-reviews — 外部评审

| 文件 | 内容 |
|------|------|
| [`kimi-0427.md`](40-reviews/kimi-0427.md) | Kimi 评审（2026-04-27） |
| [`gemini-0427.md`](40-reviews/gemini-0427.md) | Gemini 评审（2026-04-27） |
| [`gemini-0428.md`](40-reviews/gemini-0428.md) | Gemini 评审（2026-04-28） |
| [`gemini-0428-2.md`](40-reviews/gemini-0428-2.md) | Gemini 评审补充（2026-04-28） |
| [`0427-response.md`](40-reviews/0427-response.md) | 评审回复与逐条验证 |
| [`gemini-0428-response.md`](40-reviews/gemini-0428-response.md) | Gemini 评审回复 |

## 50-infra — 环境与基础设施

| 文件 | 内容 |
|------|------|
| [`hardware_benchmark.md`](50-infra/hardware_benchmark.md) | 硬件资源与性能基准 |
| [`software_versions.md`](50-infra/software_versions.md) | 软件版本记录 |
| [`cuda.md`](50-infra/cuda.md) | CUDA 环境相关 |

---

## scripts/ 目录速查

| 目录 | 用途 |
|------|------|
| `01_build/` | 体系构建、突变、最小化 |
| `02_md/` | MD 生产运行与重启 |
| `03_analysis/` | 轨迹分析（RMSD/RMSF/DCCM/PCA/cluster/MM-GBSA） |
| `04_validation/` | 跨引擎验证、AMBER→GROMACS 转换 |
| `05_docking/` | Rosetta 对接与突变扫描 |
| `06_structure/` | 结构预测分析（AF3/Chai-1） |
| `lib/` | **共享库**：`paths.py`（路径）、`stats.py`（统计）、`mda_tools.py`（MDAnalysis）、`plot_style.py`（绑图） |
| `archive/` | 旧版本、废弃实验、重复脚本 |

**新建脚本应优先使用 `lib/paths.py` 的路径常量，避免硬编码绝对路径。**
