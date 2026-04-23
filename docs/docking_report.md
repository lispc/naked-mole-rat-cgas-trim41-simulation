# cGAS-TRIM41 MD 研究项目日志

> 本文档记录项目执行过程中的关键决策、数据、估算和推理，供后续论文撰写和复盘使用。
> 
> 最后更新：2026-04-22

---


## 十、蛋白-蛋白对接（完整实验记录）

> 因 AF3 复合物预测置信度极低（ipTM < 0.25），必须通过蛋白-蛋白对接获取 cGAS-TRIM41 复合物结构。

### 10.1 尝试的工具清单（完整版）

#### 10.1.1 ClusPro 网页版

| 项目 | 详情 |
|------|------|
| 工具 | ClusPro (https://cluspro.bu.edu/) |
| 提交时间 | 2026-04-22 |
| Job ID | 1423373 |
| 受体 | TRIM41 WT full (job2/trim41_fixed.pdb, 630 aa) |
| 配体 | Hsap_cGAS_4mut (job2/cgas_fixed.pdb, 522 aa) |
| Ligand Active Residues | 463, 479, 495, 498 |
| Scoring mode | Attraction |
| 返回 poses | 2 个 (model.000.00.pdb, model.000.01.pdb) |

**❌ 失败原因**：
- 两个蛋白质心距离：Pose 1 = 83.7 Å，Pose 2 = 100.5 Å
- 活性残基到最近 TRIM41 残基距离：46-104 Å（远超 10Å 界面标准）
- 0/4 活性残基在界面上
- **推测原因**：全长对接搜索空间过大（522+630 aa），ClusPro 无法找到合理界面

**ClusPro 输出特征**：将两个蛋白合并为单链 A，残基编号 1-630（重叠），需通过残基名称区分。

**AF3 序列问题**：Job 2 实际提交的是 WT 序列（C,K,L,K），不是突变序列（S,E,Y,T）。这意味着 ClusPro 对接的也是 WT 结构。

#### 10.1.2 SDOCK2.0（本地编译，已放弃）

| 项目 | 详情 |
|------|------|
| 来源 | GitHub: victorPKU/SDOCK2.0 |
| 编译时间 | 2026-04-23 |
| 依赖 | FFTW 3.3.11 (`brew install fftw`) |
| Makefile 修改 | `FFTW_DIR = /opt/homebrew` |
| sdock.c 修改 | `FFTW_MEASURE` → `FFTW_ESTIMATE`（计划阶段耗时太长） |
| 编译命令 | `make clean && make`（成功，仅 warnings） |

**❌ 放弃过程**：
1. 第一次运行：5 分钟超时（FFTW_MEASURE 模式）
2. 修改 FFTW_MEASURE → FFTW_ESTIMATE，重新编译
3. 第二次运行：30 分钟仍无输出，example (1AY7) 也无法完成
4. **结论**：FFTW 在 Apple Silicon macOS 上性能极差，SDOCK2.0 不适合本项目

**已安装但未使用的文件**：
```
tools/SDOCK2.0/
  bin/preprocess, watmap, sdock, build
  data/ATM, so3layer.qua
```

#### 10.1.3 LightDock（最终成功方案）

| 项目 | 详情 |
|------|------|
| 安装方式 | `conda create -n py311 python=3.11 && pip install lightdock` |
| 版本 | 0.9.4 |
| 算法 | GSO (Glowworm Swarm Optimization) |
| 依赖 | Biopython 1.87, ProDy 2.6.1, FreeSASA 2.2.1, Cython 3.2.4 |
| 成功安装时间 | 2026-04-23 |

**✅ 成功原因**：
- 纯 Python + Cython，跨平台兼容性好
- 不需要 FFTW
- GSO 算法不需要预计算 FFT 计划
- macOS arm64 上编译顺利

---

### 10.2 LightDock 运行详细记录

#### 10.2.1 Hsap (人源) 结构域对接 — 第一轮（无约束）

**Setup**：
```bash
lightdock3_setup.py trim41_SPRY_413-630.pdb cgas_CT_200-554.pdb \
  -s 20 -g 20 --noxt --noh
```

| 参数 | 值 |
|------|-----|
| Swarms | 20 |
| Glowworms | 20 |
| Surface density | TotalSASA/50.00 |
| Swarm radius | 10.00 Å |
| Ligand Max Diameter | 73.23 Å |
| Surface distance | 18.31 Å |
| Swarms after filter | 546 |

**Docking**：
```bash
lightdock3.py setup.json 100 -c 1 -s fastdfire
```

| 参数 | 值 |
|------|-----|
| Steps | 100 |
| Cores | 1 |
| Scoring | fastdfire |
| 运行时间 | ~5 分钟 |

**结果**（25 poses = 5 swarms × 5 glowworms，实际只生成了 5 个 swarms 的 poses）：

| Pose | 界面残基对 | 463(Å) | 479(Å) | 495(Å) | 498(Å) | MaxDist | 状态 |
|------|-----------|--------|--------|--------|--------|---------|------|
| swarm_0/lightdock_0 | 1169 | 8.7 | 1.3 | 19.4 | 11.1 | 19.4 | ❌ |
| swarm_0/lightdock_1 | 1338 | 8.7 | 1.3 | 19.4 | 11.1 | 19.4 | ❌ |
| swarm_0/lightdock_2 | 1251 | 8.7 | 1.3 | 19.4 | 11.1 | 19.4 | ❌ |
| swarm_0/lightdock_3 | 1661 | 8.7 | 1.3 | 19.4 | 11.1 | 19.4 | ❌ |
| swarm_0/lightdock_4 | 2168 | 8.7 | 1.3 | 19.4 | 11.1 | 19.4 | ❌ |

**关键发现**：
- 所有 25 个 poses 的活性残基距离**完全相同**（8.7, 1.3, 19.4, 11.1）
- GSO 收敛到了同一个局部最小值
- 495/498 始终远离界面（~19Å），无法与 TRIM41 接触

#### 10.2.2 Hgal (裸鼹鼠) 结构域对接 — 成功

**Setup**：
```bash
lightdock3_setup.py trim41_SPRY_413-630.pdb cgas_CT_200-554.pdb \
  -s 20 -g 20 --noxt --noh
```

| 参数 | 值 |
|------|-----|
| Swarms | 20 |
| Glowworms | 20 |
| Ligand Max Diameter | 69.44 Å |
| Surface distance | 17.36 Å |
| Swarms after filter | 593 |

**Docking**：
```bash
lightdock3.py setup.json 100 -c 1 -s fastdfire
```

**结果**（20 top poses）：

| Pose | 界面残基对 | 463(Å) | 511(Å) | 527(Å) | 530(Å) | MaxDist | 状态 |
|------|-----------|--------|--------|--------|--------|---------|------|
| top_1 | 2154 | 3.8 | 3.7 | 3.7 | 3.8 | 3.8 | ✅ |
| top_2 | 2151 | 3.8 | 3.7 | 3.7 | 3.8 | 3.8 | ✅ |
| ... | ~2150 | 3.8 | 3.7 | 3.7 | 3.8 | 3.8 | ✅ |
| top_20 | 2091 | 3.8 | 3.7 | 3.7 | 3.8 | 3.8 | ✅ |

**✅ 20/20 poses 全部成功！**

#### 10.2.3 Hsap (人源) 增强对接 — 带 restraints（运行中）

**Setup**：
```bash
cat > restraints.list << EOF
L A.CYS.463
L A.LYS.479
L A.LEU.495
L A.LYS.498
EOF

lightdock3_setup.py trim41_SPRY_413-630.pdb cgas_CT_200-554.pdb \
  -s 50 -g 30 --noxt --noh \
  -r restraints.list -spr 3
```

| 参数 | 值 |
|------|-----|
| Swarms | 50 |
| Glowworms | 30 |
| Ligand restraints | 4 个 |
| Swarms per restraint | 3 |

**Docking**：
```bash
lightdock3.py setup.json 200 -c 1 -s fastdfire
```

| 参数 | 值 |
|------|-----|
| Steps | 200 |
| 预计运行时间 | ~30-40 分钟 |
| 状态 | 运行中 |

---

### 10.3 空间几何分析 — 核心发现

#### 10.3.1 活性残基在 cGAS CTD 单体上的空间分布

**Hgal (裸鼹鼠) — 紧凑分布**：

| 位点对 | CA-CA 距离 | 说明 |
|--------|-----------|------|
| 463S-511E | 18.4 Å | 最大间距 |
| 463S-527Y | 11.3 Å | |
| 463S-530T | 13.7 Å | |
| 511E-527Y | 11.3 Å | |
| 511E-530T | 11.2 Å | |
| 527Y-530T | 9.9 Å | 最小间距 |

**空间特征**：4 个位点分布在 ~18Å 球体内，形成紧凑界面补丁。

**Hsap (人源) — 分散分布**：

| 位点对 | CA-CA 距离 | 说明 |
|--------|-----------|------|
| 463C-479K | 21.1 Å | |
| **463C-495L** | **28.6 Å** | **最大间距** |
| **463C-498K** | **27.2 Å** | **最大间距** |
| 479K-495L | 11.4 Å | |
| 479K-498K | 11.3 Å | |
| 495L-498K | 9.7 Å | 最小间距 |

**空间特征**：463 与 495/498 相距 ~28Å，495/498 位于蛋白另一面。

#### 10.3.2 生物学意义

**实验观察**：裸鼹鼠 cGAS 的 TRIM41 介导泛素化更弱。

**我们的解释**：

```
人源 cGAS (WT):
  463C/479K  ←── 界面 A (可被 TRIM41 接触)
  495L/498K  ←── 界面 B (也可被 TRIM41 接触，但与 A 相距 28Å)
  → TRIM41 可能有两种结合模式，分别接触 A 或 B
  → 但无法同时接触全部 4 个位点

裸鼹鼠 cGAS (4mut):
  463S/511E/527Y/530T  ←── 聚集成单一 ~18Å 补丁
  → 4 个位点同时被 TRIM41 接触
  → 结合几何改变，可能影响：
    a) E2-Ub 复合物的空间取向
    b) 底物赖氨酸的可及性
    c) 泛素化效率（而非结合亲和力）
```

---

### 10.4 失败尝试的完整记录

| 尝试 | 时间 | 失败原因 | 学到的教训 |
|------|------|---------|-----------|
| ClusPro 全长对接 | 2026-04-22 | 搜索空间过大，蛋白相距 80+Å | 全长对接对于大蛋白不现实，应优先尝试结构域对接 |
| SDOCK2.0 | 2026-04-23 | FFTW 在 macOS 上性能极差，example 5 分钟无法完成 | FFT-based docking 工具在 Apple Silicon 上可能有兼容性问题 |
| LightDock Hsap 第一轮 | 2026-04-23 | GSO 收敛到局部最小值，495/498 始终在界面外 | 这与空间几何一致（28Å 间距），不是 docking 失败，而是物理上不可能 |
| LightDock Hsap (20sw/100step) | 2026-04-23 | 同上 | 需要增加 swarms/steps 或使用 restraints |

---

### 10.5 环境配置详情

#### Conda 环境: `py311`

```bash
conda create -n py311 python=3.11
conda activate py311
pip install lightdock
```

安装的包：lightdock 0.9.4, biopython 1.87, prody 2.6.1, freesasa 2.2.1, cython 3.2.4, numpy 2.4.4, scipy 1.17.1

#### Homebrew 包

```bash
brew install fftw    # 3.3.11, libomp 22.1.4 (为 SDOCK2.0 安装，最终未使用)
```

---

### 10.6 文件清单（完整版）

#### ClusPro 结果
```
structures/docking/cluspro/
  job2_Hsap_4mut_full/
    cluspro.1423373.tar.bz2
    cluspro.1423373/
      model.000.00.pdb          # 11286 atoms, 单链 A
      model.000.01.pdb
      pose_1_split.pdb          # 手动拆分: TRIM41(603 res) + cGAS(522 res)
      pose_2_split.pdb
    docking_analysis_report.md
```

#### SDOCK2.0 安装
```
tools/SDOCK2.0/                  # GitHub: victorPKU/SDOCK2.0
  src/
    sdock.c                      # 修改: FFTW_MEASURE → FFTW_ESTIMATE
    Makefile                     # 修改: FFTW_DIR = /opt/homebrew
  bin/
    preprocess, watmap, sdock, build
  data/
    ATM, so3layer.qua
  example/
    1AY7/, 2REX/
```

#### LightDock 结果
```
structures/docking/lightdock/
  Hsap_domain/
    setup.json
    init/
    swarm_0/ 到 swarm_19/        # 20 swarms
      gso_0.out 到 gso_100.out
      lightdock_0.pdb 到 lightdock_4.pdb
    rank_by_scoring.list
    rank_by_luciferin.list
    rank_by_rmsd.list
    top_1.pdb 到 top_20.pdb

  Hgal_domain/
    setup.json
    init/
    swarm_0/ 到 swarm_19/
    rank_by_scoring.list
    top_1.pdb 到 top_20.pdb
    best_pose.pdb                # 复制自 top_1.pdb

  Hsap_restrained/               # 运行中
    setup.json
    restraints.list              # L A.CYS.463, L A.LYS.479, L A.LEU.495, L A.LYS.498
    swarm_0/ 到 swarm_49/        # 50 swarms

  analyze_lightdock.py           # 分析脚本
```

---

