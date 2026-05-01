# cGAS-TRIM41 MD 研究项目日志

> 本文档记录项目执行过程中的关键决策、数据、估算和推理，供后续论文撰写和复盘使用。
> 
> 最后更新：2026-04-22

---


## 软件版本记录

### macOS / Apple M3 Pro (原始环境)

| 软件 | 版本 | 来源 | 用途 |
|------|------|------|------|
| Python | 3.13.13 | conda-forge (cgas-md) | MD 分析 |
| Python | 3.11.12 | conda-forge (py311) | LightDock |
| OpenMM | 8.5.1 | conda-forge (apple build) | MD 模拟 |
| AmberTools | 24.8 | conda-forge | 体系构建 |
| MDAnalysis | 2.10.0 | conda-forge | 轨迹分析 |
| MDTraj | 1.11.1 | conda-forge | 轨迹分析 |
| Biopython | 1.87 | conda-forge / pip | 结构处理 |
| LightDock | 0.9.4 | pip | 蛋白对接 |
| NumPy | 2.4.3 / 2.4.4 | conda-forge / pip | 数值计算 |
| SciPy | 1.17.1 | conda-forge / pip | 科学计算 |
| Matplotlib | 3.10.8 | conda-forge | 绘图 |
| FFTW | 3.3.11 | Homebrew | SDOCK2.0 (未使用) |
| GCC/Clang | Apple clang 21.0.0 | Xcode | 编译 |

### Linux / 4× RTX 3090 (迁移后环境，2026-04-23)

| 软件 | 版本 | 来源 | 用途 | 备注 |
|------|------|------|------|------|
| Python | 3.11.15 | conda-forge (cgas-md) | MD 模拟 | 迁移后统一用 3.11 |
| OpenMM | 8.5.1 | conda-forge (CUDA build) | MD 模拟 | CUDA 12.x |
| AmberTools | 24.8 | conda-forge | 体系构建 | tleap, cpptraj |
| MDAnalysis | 2.10.0 | conda-forge | 轨迹分析 | |
| MDTraj | 1.11.1 | conda-forge | 轨迹分析 | |
| Biopython | 1.87 | conda-forge | 结构处理 | |
| NumPy | 2.4.3 | conda-forge | 数值计算 | |
| SciPy | 1.17.1 | conda-forge | 科学计算 | |
| Matplotlib | 3.10.8 | conda-forge | 绘图 | |
| Seaborn | 0.13.2 | conda-forge | 绘图 | |
| Pandas | 2.3.3 | conda-forge | 数据处理 | |
| OpenMMTools | 0.26.0 | conda-forge | 测试体系 | verify_openmm.py 依赖 |
| CUDA Driver | 580.95.05 | NVIDIA | GPU 驱动 | |
| CUDA Runtime | 12.9 | conda-forge | CUDA 库 | |
| **Boltz-2** | **2.2.1** | **pip (MIT)** | **结构预测** | **本地 AF3-level，支持蛋白/DNA/小分子/亲和力预测** |
| PyTorch | 2.11.0 | pip (CUDA 13) | Boltz-2 后端 | |

**Boltz-2 环境配置：**
- 独立 conda 环境：`boltz`
- CUDA 13 库路径：`/home/scroll/miniforge3/envs/boltz/lib/python3.11/site-packages/nvidia/cu13/lib`
- 已设置 `LD_LIBRARY_PATH` 环境变量（通过 `conda env config vars set`）
- 模型权重缓存：`~/.boltz/`（~5.5 GB）

**Boltz-2 使用示例：**
```bash
conda activate boltz
boltz predict input.yaml --use_msa_server --out_dir output/
```

**输入格式：** YAML（支持蛋白、DNA、RNA、小分子、共价修饰、亲和力预测）

---

*文档创建：2026-04-22*
*最后更新：2026-04-23（添加 Linux/RTX 3090 环境）*
*维护者：Kimi Code CLI*
