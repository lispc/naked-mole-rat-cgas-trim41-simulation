# CUDA 机器迁移指南

> 从 Apple Silicon (macOS, OpenCL) 迁移到 NVIDIA GPU (Linux, CUDA) 的注意事项。
>
> **最后更新：2026-04-23** — 已成功完成迁移，3× RTX 3090 并行 MD 运行中。

---

## 一、环境清单

### 1.1 Conda 环境创建

```bash
# 创建环境（建议 Python 3.11-3.12，OpenMM 对 3.13 支持尚在完善）
conda create -n cgas-md python=3.11
conda activate cgas-md

# 核心包
conda install -c conda-forge openmm=8.1  # 或更新版本
conda install -c conda-forge ambertools=24
conda install -c conda-forge mdanalysis=2.10
conda install -c conda-forge mdtraj=1.11
conda install -c conda-forge matplotlib seaborn pandas
conda install -c conda-forge biopython=1.87
conda install -c conda-forge pymol-open-source  # 或 schrodinger/pymol

# Optional: MDTraj/MDAnalysis 加速
conda install -c conda-forge pyopencl  # 如需要 OpenCL fallback
```

### 1.2 CUDA 驱动验证

```bash
# 检查 NVIDIA 驱动
nvidia-smi

# 检查 CUDA 版本
nvcc --version

# 检查 PyTorch/TensorFlow 中的 CUDA（仅作参考）
python -c "import torch; print(torch.cuda.is_available())"
```

### 1.3 OpenMM CUDA 验证

```bash
conda activate cgas-md
python scripts/verify_openmm.py
```

预期输出：
```
✅ CUDA backend: 10k steps in X.XXs (XXXX steps/s)
   Potential energy: -XXXX kJ/mol
```

如果 CUDA 不出现：
1. 检查 `conda list | grep openmm` 确认是 `cuda` 构建版本
2. 可能需要 `conda install -c conda-forge openmm=*=*cuda*` 强制 CUDA 版本
3. 检查 `LD_LIBRARY_PATH` 是否包含 CUDA 库

---

## 二、代码适配

### 2.1 平台选择（已自动处理）

所有脚本已内置平台自动检测：

```python
# scripts/build_system.py 和 scripts/run_md.py
get_platform(platform_name='auto')
# 自动选择: CUDA > OpenCL > CPU
```

**无需修改代码**，直接运行即可。

### 2.2 显式指定 CUDA

如需强制使用 CUDA（例如多 GPU 环境）：

```bash
python scripts/run_md.py \
  --prmtop data/md_runs/Hgal_domain/Hgal_domain.prmtop \
  --pdb data/md_runs/Hgal_domain/Hgal_domain_minimized.pdb \
  --name Hgal_domain_rep1 \
  --outdir data/md_runs/Hgal_domain/rep1 \
  --prod-ns 200 \
  --platform CUDA
```

### 2.3 多 GPU 选择

如需指定特定 GPU（例如 RTX 3090 #1）：

```bash
export CUDA_VISIBLE_DEVICES=1  # 使用第2张卡
python scripts/run_md.py ... --platform CUDA
```

或在代码中修改 `CudaDeviceIndex`：

```python
# scripts/run_md.py, get_platform() 函数
props['CudaDeviceIndex'] = '1'  # 改为目标 GPU 编号
```

---

## 三、性能预期

### 3.1 RTX 3090 vs Apple M3 Pro

| 指标 | Apple M3 Pro (OpenCL) | RTX 3090 (CUDA) | 加速比 |
|------|----------------------|-----------------|--------|
| NVT Heating | ~46 ns/day | ~150-200 ns/day | **3-4x** |
| NPT Equil | ~21 ns/day | ~100-150 ns/day | **5-7x** |
| NVT Production | ~40-50 ns/day (est.) | ~200-300 ns/day | **5-6x** |

> 注：RTX 3090 有 24GB VRAM，本体系 116k atoms 约需 2-3GB，完全无压力。

### 3.2 200ns × 3 重复耗时估算

| 机器 | 单条 200ns | 3 重复总计 |
|------|-----------|-----------|
| M3 Pro | ~4-5 天 | ~12-15 天 |
| RTX 3090 | ~1 天 | ~3 天 |

---

## 四、文件同步检查

### 4.1 通过 git 同步的内容

git 已 add（未 commit）的内容包括：

```
.gitignore
README.md
docs/          # 所有文档（10+ .md 文件）
scripts/       # 所有 .py 和 .pml 脚本
sequences/     # FASTA 文件
figures/       # 发表级图表
structures/docking/lightdock/Hgal_domain/best_pose.pdb  # MD 起始结构
structures/docking/lightdock/*/setup.json               # LightDock 配置
```

### 4.2 需要重新生成的内容（不在 git 中）

```bash
# 1. MD 体系文件（在目标机器上重新构建）
data/md_runs/Hgal_domain/
  Hgal_domain.prmtop          # tleap 生成
  Hgal_domain.rst7            # tleap 生成
  Hgal_domain_minimized.pdb   # OpenMM EM 生成

# 构建命令
python scripts/build_system.py \
  --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
  --name Hgal_domain \
  --outdir data/md_runs/Hgal_domain \
  --platform CUDA

# 2. 分析输出（运行后生成）
data/md_runs/Hgal_domain/rep1/
  *.dcd, *.log, *.chk
```

### 4.3 AF3 原始结构（可选）

如果需要重新分析 AF3 结果或重新截断结构域：

```bash
# 下载 4 个 job 的 zip 文件到 structures/af3_raw/
# 然后运行
python scripts/process_af3_results.py --job-dir structures/af3_raw/job1_Hsap_WT
```

或直接复制现有机器上的 `structures/af3_raw/*/cgas_CT_200-554.pdb` 等截断后文件。

---

## 五、常见问题

### Q1: `CUDA platform not available`

```bash
# 检查 conda 安装的是否为 CUDA 版本
conda list openmm
# 应显示类似: openmm 8.1.1 py311h..._cuda12_0

# 如果不是 CUDA 版本，重装
conda install -c conda-forge openmm=*=*cuda* -y
```

### Q2: `mixed` precision 报错

RTX 3090 支持 `mixed` precision（默认已设置）。如果出现 NaN：

```python
# scripts/run_md.py, get_platform()
props['CudaPrecision'] = 'single'  # 改为 single 试试
```

### Q3: 显存不足 (OOM)

116k atoms 体系在 RTX 3090 上约需 2-3GB，24GB 绰绰有余。如遇 OOM：

1. 检查是否同时运行了其他 GPU 任务
2. 减小 PME grid 或增大 `ewaldErrorTolerance`
3. 尝试 `CpuPme`：在 platform properties 中加 `{'UseCpuPme': 'true'}`

### Q4: `run_md.py` 不支持 `--platform CUDA`

原脚本中 `--platform` 的 `choices` 只包含 `['OpenCL', 'CPU', 'Reference']`，缺少 `CUDA`。

**修复**（已提交）：
```python
parser.add_argument('--platform', default='auto', choices=['auto', 'CUDA', 'OpenCL', 'CPU', 'Reference'])
```

### Q5: 轨迹文件兼容性

DCD 轨迹是跨平台的。在 macOS 上生成的 DCD 可在 Linux 上直接用 MDAnalysis 读取：

```python
import MDAnalysis as mda
u = mda.Universe('system.prmtop', 'trajectory.dcd')
```

---

## 六、推荐运行顺序

```bash
# 1. 克隆/同步代码
git clone <repo> && cd <repo>
git checkout <branch>  # 或直接使用当前状态

# 2. 创建环境（见 1.1）
conda create -n cgas-md python=3.11 ...

# 3. 验证 CUDA
conda activate cgas-md
python scripts/verify_openmm.py

# 4. 构建体系
python scripts/build_system.py \
  --pdb structures/docking/lightdock/Hgal_domain/best_pose.pdb \
  --name Hgal_domain \
  --outdir data/md_runs/Hgal_domain \
  --platform CUDA

# 5. 生产运行（3 重复，后台）
mkdir -p data/md_runs/Hgal_domain/rep{1,2,3}

# Rep 1
python scripts/run_md.py \
  --prmtop data/md_runs/Hgal_domain/Hgal_domain.prmtop \
  --pdb data/md_runs/Hgal_domain/Hgal_domain_minimized.pdb \
  --name Hgal_domain_rep1 \
  --outdir data/md_runs/Hgal_domain/rep1 \
  --prod-ns 200 --platform CUDA > rep1.log 2>&1 &

# Rep 2 (用不同 seed)
python scripts/run_md.py ... --name Hgal_domain_rep2 ... > rep2.log 2>&1 &

# Rep 3
python scripts/run_md.py ... --name Hgal_domain_rep3 ... > rep3.log 2>&1 &

# 监控
tail -f data/md_runs/Hgal_domain/rep1/Hgal_domain_rep1_prod.log
```

---

## 七、注意事项

1. **不要 commit 轨迹文件**：DCD 文件通常数百 MB 到数 GB，已加入 `.gitignore`
2. **定期 checkpoint**：脚本已内置每 1ns 自动保存 checkpoint
3. **温度监控**：RTX 3090 长时间运行温度可能达到 70-80°C，确保机箱散热良好
4. **电源**：单卡满载功耗约 350W，4 卡建议 1600W+ 电源
5. **不同 seed**：3 重复使用不同随机种子，通过 `CUDA_VISIBLE_DEVICES` + 不同启动时间自动实现
6. **业务逻辑监控**：已部署 `scripts/monitor_md.sh`，每小时检查能量稳定性、温度、体积收敛、DCD 文件完整性等

## 八、实际迁移经验总结

### 迁移耗时
| 步骤 | 耗时 | 备注 |
|------|------|------|
| 安装 Miniforge | ~2 min | 自动下载安装 |
| 创建 conda env + 安装包 | ~10 min | 3GB 下载，mamba 加速 |
| OpenMM CUDA 验证 | <1 min | `verify_openmm.py` |
| 构建 MD 体系 (tleap + EM) | ~2 min | 含溶剂化、能量最小化 |
| 启动 3 重复 MD | <1 min | 并行后台运行 |
| **总计** | **~15 min** | 从零到生产运行 |

### 实测 vs 预估对比
| 指标 | 预估 | 实测 | 偏差 |
|------|------|------|------|
| RTX 3090 速度 | 200-300 ns/day | 152 ns/day | 偏低（Python I/O 开销大）|
| 单条 200ns 耗时 | ~1 天 | ~1.3 天 | 基本吻合 |
| 3 重复并行总耗时 | ~3 天 | ~1.3 天 | 优于预估（3 卡并行）|
