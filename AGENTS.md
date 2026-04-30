# Agent Instructions

## 软件环境

### Conda 环境

本项目使用 Miniforge 作为 conda 发行版，安装位置：`~/miniforge3`

| 环境 | 路径 | 用途 |
|------|------|------|
| `cgas-md` | `~/miniforge3/envs/cgas-md` | MD 模拟、轨迹分析、OpenMM、AmberTools |
| `gmx` | `~/miniforge3/envs/gmx` | GROMACS 2026.0 CUDA（独立环境，避免与 AmberTools 冲突） |
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

### GROMACS

GROMACS 2026.0 (nompi_cuda) 安装在独立环境 `gmx` 中：
```bash
conda activate gmx
gmx --version
```

**关键限制**：
- `cgas-md` 环境同时装有 GROMACS 2025.4 CPU 版（`ambertools` 与 `gromacs=2026.0=nompi_cuda` 有依赖冲突，不能共存）
- GPU MD 必须在 `gmx` 环境中运行
- AMBER → GROMACS 拓扑转换（`convert_amber_to_gromacs.py`）使用 `cgas-md` 环境的 MDAnalysis + parmed

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

---

## 项目结构速查

```
├── docs/               # 项目日志、分析报告、计划文档
├── scripts/            # 全部 Python / Shell / PyMOL 脚本（见 scripts/README.md）
├── sequences/          # FASTA 序列、AF3 提交记录
├── structures/         # 结构文件、对接结果
├── figures/            # 分析图表
├── data/
│   ├── md_runs/        # OpenMM MD 运行目录 (Hsap_WT, Hsap_4mut, Hgal_WT, Hgal_4mut_rev)
│   ├── md_runs_gmx/    # GROMACS MD 运行目录
│   ├── analysis/       # 分析输出 (batch, clustering, c1_c4, etc.)
│   └── us/             # 伞形采样窗口数据
├── main.pdf            # 项目主文档
└── README.md           # 项目简介
```

---

## MD 运行规范

### OpenMM
- 力场：ff19SB + OPC
- 积分器：LangevinMiddle，摩擦 1.0/ps，步长 2fs
- 温度：300K，压力：1 bar (MonteCarloBarostat)
- 非键截断：1.0 nm，PME，ewaldErrorTolerance=0.0005
- 约束：HBonds

### GROMACS
- 力场：通过 parmed 转换保留 AMBER 参数（**转换后必须运行 `fix_gromacs_cmap.py` 修复 CMAP 残基特异性**）
- 积分器：md (leap-frog)，步长 2fs
- 温度：v-rescale，tau_t=1.0 ps
- 压力：Production 阶段 **NVT** (`pcoupl = no`)，与 OpenMM 一致
- 非键截断：PME，rcoulomb=1.0 nm，rvdw=1.0 nm，vdw-modifier = Potential-shift-Verlet
- 约束：LINCS，`lincs-iter=2, lincs-order=6`（提高精度）

### GROMACS CMAP 修复
- **问题**：parmed 将 ff19SB 的 14 种残基特异性 CMAP types 错误压缩为 1 种 (`C N XC C N`)
- **后果**：RMSD 系统性偏高 ~4×，解离行为定性不同
- **修复**：`python scripts/fix_gromacs_cmap.py --top input.top --prmtop ref.prmtop --out fixed.top`
- **原理**：将 CA 类型从 `XC` 改为 `XC{n}` (n=0-13)，在 `[cmaptypes]` 中定义 14 种独立的 grid

### 轨迹保存
- OpenMM：每 10 ps 一帧（`reportInterval=5000`），DCD 格式
- GROMACS：每 10 ps 一帧（`nstxout-compressed=5000`），XTC 格式

### GROMACS 分析注意事项
- GROMACS 的 `.top` 文件经 parmed 转换后缺少 `%VERSION` 头，**MDAnalysis 无法直接读取**。应使用 `gmx editconf -f prod.tpr -o topol.gro` 生成 gro 拓扑，再用 MDAnalysis 读取 `gro + xtc`
- GROMACS 中每个 moleculetype 的 resid 从 1 重新开始。因此 cGAS=resindex 0–217，TRIM41=resindex 218–540，不能用 `resid 219-541` 选择 TRIM41
- GROMACS 的 XTC 坐标精度为 0.001 nm（0.01Å），对 RMSD 分析影响可忽略
- `convert_amber_to_gromacs.py` 使用 MDAnalysis 直接从 DCD 读取坐标+box，避免 `cpptraj → rst7` 的 box 信息丢失问题

### GROMACS 分析性能陷阱
- **XTC + `in_memory=True` 是致命组合**：MDAnalysis 的 `AlignTraj(..., in_memory=True)` 对 XTC 压缩格式需要逐帧解压缩 85k 原子到内存，4 replica × 20,000 帧可耗时 **10+ 小时**
- **推荐方案**：使用 `batch_analyze_hsap_gmx_fast.py` 的流式读取 + numpy Kabsch 对齐，全程仅需 **~15 分钟**
- **DCD vs XTC**：OpenMM 的 DCD 是未压缩二进制，可以内存映射随机访问；GROMACS 的 XTC 需要解压缩，I/O 开销大

---

## 磷酸化 MD 规范

- 磷酸化残基命名（生理 pH 7.4，双负离子）：`SER → SEP`, `THR → TPO`, `TYR → PTR`
- 力场加载：`source leaprc.phosaa19SB`（ff19SB 优化参数）
- 结构准备：PyMOL Builder 添加 PO₃ → 修改残基名 → `pdb4amber --reduce` → `tleap`
- 关键位点：S305（cGAS NTase 域，CHK2 磷酸化促进 TRIM41 结合）
- 注意：`pdb4amber`/`reduce` **不会自动识别磷酸化**，必须手动添加磷酸基团并修改残基名
- 详细方案：`docs/phosphorylation_md_plan.md`

## 分析规范

- 批量分析脚本默认输出到 `data/analysis/` 下对应子目录
- 图表统一保存为 PNG (dpi=300)
- 数值结果优先保存为 `.npz`（NumPy）或 `.json`
- 关键发现需同步更新 `docs/project_log.md`

---

## 脚本索引

完整脚本说明见 `scripts/README.md`。
