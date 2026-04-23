# cGAS-TRIM41 MD 研究项目日志

> 本文档记录项目执行过程中的关键决策、数据、估算和推理，供后续论文撰写和复盘使用。
> 
> 最后更新：2026-04-22

---


## 十三、软件版本记录（更新）

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

---

*文档创建：2026-04-22*
*最后更新：2026-04-23*
*维护者：Kimi Code CLI*
