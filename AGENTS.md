# Agent Instructions

## 软件环境

### Conda 环境

本项目使用 Miniforge 作为 conda 发行版，安装位置：`~/miniforge3`

| 环境 | 路径 | 用途 |
|------|------|------|
| `cgas-md` | `~/miniforge3/envs/cgas-md` | MD 模拟、轨迹分析 |
| `rosetta` | `~/miniforge3/envs/rosetta` | PyRosetta、Rosetta 对接 |
| `base` | `~/miniforge3` | 基础环境 |

**默认使用 conda 环境中的 Python，不要直接使用系统 Python**。

正确示例：
```bash
~/miniforge3/bin/conda run -n cgas-md python3 script.py
# 或先激活环境
conda activate cgas-md
python3 script.py
```

### PyMOL

PyMOL 安装在 `cgas-md` 环境中：
```bash
~/miniforge3/bin/conda run -n cgas-md pymol -cq script.pml
```

---

## Git 工作流

**除非用户明确要求，否则不要执行 `git add` 或 `git commit`**。

- 可以生成文件、修改代码、编辑文档
- 可以查看 git status、git diff
- 不要执行任何会修改 git 索引或历史记录的操作（`git add`、`git commit`、`git push`、`git reset`、`git rebase` 等）
