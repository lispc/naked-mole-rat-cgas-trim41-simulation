# Agent Instructions

## 环境

| conda 环境 | 用途 |
|-----------|------|
| `cgas-md` | OpenMM、AmberTools、分析、PyMOL |
| `gmx` | GROMACS 2026.0 CUDA（独立，避免与 AmberTools 冲突） |
| `rosetta` | PyRosetta、Rosetta 对接 |

**默认用 conda 环境，不要直接用系统 Python。**

```bash
# 运行 Python 脚本
conda run -n cgas-md python3 scripts/<script>.py

# GROMACS
conda activate gmx && gmx mdrun ...

# PyMOL
conda run -n cgas-md pymol -cq script.pml
```

## 项目结构速查

| 目录 | 内容 |
|------|------|
| `scripts/` | 全部脚本（索引见 `scripts/README.md`） |
| `data/md_runs/` | OpenMM MD 运行目录 |
| `data/md_runs_gmx/` | GROMACS MD 运行目录 |
| `data/analysis/` | 分析输出 |
| `docs/` | 日志、报告、方案（见下方索引） |
| `sequences/` | FASTA 序列 |
| `structures/` | 结构文件、对接结果 |
| `figures/` | 分析图表 |

## Git 边界

**除非用户明确要求，否则不要执行 `git add`、`git commit`、`git push`、`git reset`、`git rebase`。**

- ✅ 可以生成/修改文件、查看 `git status`/`git diff`
- 🚫 不要修改 git 索引或历史

## 边界

- ✅ **Always**: 用 conda 环境运行脚本；修改代码后同步更新相关 `docs/` 文档
- ⚠️ **Ask first**: git commit、安装新包、删除运行数据、修改 `data/md_runs/*` 中的轨迹
- 🚫 **Never**: 用系统 Python、直接修改 `data/md_runs/` 中的二进制轨迹文件

## 参考索引

| 主题 | 文档 |
|------|------|
| 脚本说明 | `scripts/README.md` |
| 脚本规范 | `scripts/AGENTS.md` |
| MD 规范 & GROMACS 诊断 | `docs/30-diagnostics/gromacs_openmm_divergence_diagnosis.md` |
| 磷酸化方案 | `docs/20-protocols/phosphorylation_md_plan.md` |
| 项目日志 | `docs/00-project/project_log.md` |
| 项目文档总览 | `docs/README.md` |
