# Scripts Index

本目录包含项目所有分析、模拟、对接和监控脚本。按功能模块分组。

---

## 1. 体系构建 (System Building)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `build_system.py` | 从 LightDock pose 用 tleap 构建溶剂化体系 (ff19SB + OPC)，OpenMM 验证 | Active |
| `build_system_v2.py` | 改进版 builder，支持截角八面体、可配置缓冲层、保守 EM + 骨架位置限制 | **推荐** |
| `build_hsap_md_system.py` | 从 Rosetta docking pose 构建 Hsap 体系 (pdb4amber + tleap) | Active，基本被 v2 取代 |
| `convert_amber_to_gromacs.py` | AMBER prmtop + DCD → GROMACS `.top` + `.gro` (MDAnalysis + parmed，保留 EPW)。**转换后自动修复 CMAP 残基特异性** | Active |
| `fix_gromacs_cmap.py` | 修复 parmed 转换丢失的 ff19SB CMAP 残基特异性：14 types → 14 residue-specific `XC{n}` atom types + cmaptypes | Active |
| `prepare_rosetta_input.py` | 合并受体/配体 PDB，分配 chain ID，供 Rosetta 对接使用 | Active |
| `prepare_4mut_docking_input.py` | 将突变体 cGAS CTD 通过 Kabsch 对齐到 WT 坐标，再与 TRIM41 SPRY 合并 | Active |
| `verify_openmm.py` | 验证 OpenMM 安装并 benchmark 各平台 (CUDA/OpenCL/CPU) | Active |

---

## 2. MD 执行 (MD Execution)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `run_md.py` | OpenMM 生产流程：升温 (0→300K) → NPT 平衡 → NVT 生产 | Active |
| `run_replicas.py` | 多 replica 启动（不同随机种子、多 GPU） | **Broken**（传了 `--gpu` 但 `run_md.py` 不支持） |
| `run_gmx_replica.sh` | 单个 GROMACS replica 全流程：EM → NVT → NPT → Production | Active |
| `run_gmx_batch.sh` | 并行启动 4 个 GROMACS replica（WT/4mut × rep1/2，GPU 0–3） | Active |

---

## 3. 伞形采样 & WHAM (Umbrella Sampling)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `run_umbrella_sampling.py` | 单窗口 US：谐振子限制 RING–Lys315 距离，EM→NPT→NVT | Active |
| `run_us_simple.py` | 简化 US 窗口（EM→warm-up→production，无显式 NPT barostat） | Active |
| `launch_us_4mut.sh` | Hsap_4mut US 窗口批量启动器，每批 4 窗口 × 4 GPU | Active |
| `us_scheduler.sh` | 根据 GPU 利用率自动调度 US 窗口 | Active |
| `us_auto_launch.sh` | 增强版队列启动器：读取窗口列表、跳过繁忙 GPU、检测已完成窗口 | Active |
| `run_wham.py` | WHAM PMF 重建 + bootstrap 误差估计 | Active |
| `test_us_k_values.py` | 测试不同谐振子力常数，跑短模拟比较 CV 分布 | Active |
| `analyze_close_state.py` | 分析 US 近态构象 (CV < threshold)，提取 Lys-334 环境接触 | Active |

---

## 4. 对接 & Rosetta (Docking)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `rosetta_dock_4mut.py` | PyRosetta `DockingProtocol` 对接 Hgal_4mut_rev / Hsap_4mut | Active，硬编码 macOS 路径 |
| `rosetta_mutational_scan.py` | Pack+Min 突变扫描，计算 WT/单点/组合 4mut 界面 ΔΔG | Active |
| `rosetta_fastrelax_mutscan.py` | FastRelax 突变扫描，捕捉变构效应对界面能的影响 | Active，硬编码 macOS 路径 |
| `analyze_docking_results.py` | 解析 Rosetta scorefile，计算 CA-RMSD，分类 PASS/GRAY/FAIL | Active |
| `quick_analyze_top_poses.py` | 快速检查 LightDock top poses 活性残基是否在 TRIM41 界面 10Å 内 | Active |

---

## 5. AlphaFold3 分析

| 脚本 | 目的 | 状态 |
|------|------|------|
| `analyze_af3.py` | 提取 AF3 预测 ipTM、pTM、pLDDT、界面残基、突变位点置信度 | Active |
| `process_af3_results.py` | 自动解压 AF3 结果、解析置信度、提取单体、PDBFixer 修复、绘图 | Active |

---

## 6. 单体系轨迹分析 (Single-System Analysis)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `analyze_system.py` | 通用单体系分析：RMSD、RMSF、界面接触占有率、COM、活性位点距离 | Active |
| `analyze_trajectory.py` | 综合轨迹分析：RMSD、RMSF、H-bond、界面接触、PCA（部分功能 stub） | Active |
| `analyze_lys_ubiquitination.py` | 分析 cGAS 赖氨酸泛素化潜力：到 TRIM41 RING 距离、溶剂暴露、RMSF | Active |
| `calc_residue_distances.py` | 计算 Hgal/Hsap cGAS CTD 四个活性残基的 CA–CA 距离 | Active |
| `run_dccm.py` | 动态交叉相关图 (DCCM)，提取高相关残基对 | Active |
| `run_pca.py` | 主成分分析，生成特征值谱、PC 投影、极端构象 PDB | Active |
| `run_mmpbsa.py` | MM-GBSA/PBSA (AmberTools `MMPBSA.py`)，DCD→NetCDF 转换 | Active，但脚本需根据手动 workflow 重写 |

---

## 7. 批量 & 对比分析 (Batch & Comparative)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `batch_analyze_hsap.py` | 6 replica (WT+4mut × 3 rep) 完整批量分析，含 O(N²) 界面指标、GMM 聚类 | Active |
| `batch_analyze_hsap_v2.py` | 轻量优化版，仅用 CA 界面指标，避免重循环 | **推荐** |
| `batch_analyze_hsap_gmx.py` | GROMACS 版批量分析（适配 gro+xtc，`resindex` 替代 `resid`） | Active（慢，XTC + in_memory 瓶颈） |
| `batch_analyze_hsap_gmx_fast.py` | GROMACS 快速版：流式读取 + numpy Kabsch 对齐，避免 in_memory 陷阱 | **推荐** |
| `batch_analyze_hsap_gmx_fast_v2.py` | 待创建：修复 Kabsch 对齐 bug（原脚本使用 `R.T` 而非 `R`），支持修复后 CMAP 拓扑分析 | Planned |
| `quick_compare_rep1.py` | WT rep1 vs 4mut rep1 快速并排对比 (RMSD/COM/Rg/RMSF) | Active |
| `cluster_analysis.py` | 基于 GMM 的聚类分析（COM、min CA、Rg、RMSD、突变位点距离） | Active |
| `analyze_c1_c4.py` | C1 vs C4 亚稳态深度构象比较：界面残基对、突变位点环境、RMSF | Active |
| `compare_systems.py` | 双体系对比（WT vs mutant），有效样本量 t 检验 | Active |
| `compare_four_systems.py` | 四体系综合对比（Hgal_WT、Hgal_4mut_rev、Hsap_WT、Hsap_4mut） | Active |

---

## 8. 监控 & 工具 (Monitoring & Utilities)

| 脚本 | 目的 | 状态 |
|------|------|------|
| `monitor_md.sh` | MD 物理检查：温度曲线、NPT 体积收敛、能量漂移、NaN 检测、ETA | Active |
| `monitor_md_loop.sh` | 每小时循环运行 `monitor_md.sh` | Active |
| `watch_and_launch.py` | 监控 Hgal MD 完成，自动启动后续分析和下一批 MD | Active |
| `test_one_replica.py` | 单 replica 性能/加载测试：对齐、RMSD、COM 计时 | Active |

---

## 9. PyMOL 可视化脚本

| 脚本 | 目的 |
|------|------|
| `visualize_active_residues.pml` | 高亮 cGAS 活性残基与 TRIM41 界面 |
| `visualize_4mut_allostery.pml` | 4mut 变构效应可视化 |
| `visualize_4mut_allostery_v2.pml` | 改进版 4mut 变构可视化 |
| `visualize_corrected_4mut.pml` | 修正后的 4mut 构象可视化 |
| `visualize_hgal_interface.pml` | Hgal 界面可视化 |

---

## 备注

- **Active**：当前可用，近期维护
- **Broken**：已知问题，需修复后再用（如 `run_replicas.py`）
- **推荐**：同功能多版本中的首选（如 `build_system_v2.py`）
- 所有 Python 脚本应在 `cgas-md` conda 环境中运行
- GROMACS 相关脚本应在 `gmx` conda 环境中运行
