# cGAS-TRIM41 MD 研究项目日志

> 本文档记录项目执行过程中的关键决策、数据、估算和推理，供后续论文撰写和复盘使用。
> 
> 最后更新：2026-04-22

---

## 一、项目背景与目标

### 1.1 科学问题
基于 Chen et al., Science 2025（*A cGAS-mediated mechanism in naked mole-rats potentiates DNA repair and delays aging*），裸鼹鼠（Naked Mole-Rat, NMR）cGAS 蛋白存在 **4 个氨基酸变异**（位于 C-端结构域），导致其功能发生逆转：
- 人类/小鼠 cGAS：**抑制** 同源重组（HR）修复
- 裸鼹鼠 cGAS：**促进** HR 修复

实验因果链：
```
4个氨基酸变异
    ↓
TRIM41 介导的泛素化减弱（裸鼹鼠）/ 增强（人类）
    ↓
K48-泛素化水平变化
    ↓
P97 提取效率变化
    ↓
染色质滞留时间变化
    ↓
FANCI-RAD50 相互作用
    ↓
HR 修复效率变化
```

### 1.2 计算目标
通过分子动力学模拟和计算突变扫描，在原子水平上解释/验证这 4 个变异如何影响 cGAS-TRIM41 相互作用，进而影响泛素化效率。

### 1.3 论文定位
**选项 A：作为实验论文的 In silico validation**（已确认）
- 配合实验数据发 Cell Reports / Nature Communications 等级别
- 计算量适中，5 周内可完成
- 不作为独立计算论文（后者需要更多方法创新和对照，时间紧张）

---

## 二、关键决策记录

### 决策 1：结构预测途径 → AlphaFold3 Server（网页版）

| 选项 | 优劣 | 选择 |
|------|------|------|
| AF3 Server | 最准确，免费，但需手动提交 | **✅ 选中** |
| 本地 ColabFold | M3 Pro 36GB 对双链预测紧张，慢且 OOM 风险 | ❌ 放弃 |
| 本地 Boltz-1 | 开源 AF3 级精度，可本地运行 | ❌ 备用方案 |

**Rationale**：
- AF3 Server 是目前蛋白-蛋白复合物预测的 SOTA
- 免费账户每日限 20 个预测，4 个 job 绰绰有余
- 本地运行大模型对 36GB 统一内存压力过大

### 决策 2：MD 体系策略 → 智能截取结构域

| 选项 | 优劣 | 选择 |
|------|------|------|
| 截取结构域 | 原子数减半，可做 200ns×3 重复，统计 robust | **✅ 选中** |
| 全长复合物 | 更完整，但只能 100ns×2-3 重复，统计弱 | ❌ 放弃 |

**Rationale**：
- 全长 cGAS (522 aa) + TRIM41 (630 aa) + 水 ≈ 25-35 万原子
- 在 M3 Pro OpenCL 上预估速度：30-60 ns/day
- 4 体系 × 3 重复 × 200 ns = 2400 ns 总轨迹
- 全长串行需 **50-80 天**，远超"几周"预算
- 截取后 (~500-600 aa 保留互作结构域) ≈ 8-12 万原子，速度可翻倍
- E3 连接酶的底物识别通常发生在 C-terminal 结构域，截断科学上合理

### 决策 3：突变扫描工具 → Rosetta（不申请 FoldX）

| 选项 | 优劣 | 选择 |
|------|------|------|
| Rosetta | 完全免费，功能强，安装包 ~2GB | **✅ 选中** |
| FoldX | 需申请学术 license（1-2 天），操作方便 | ❌ 不申请 |

**Rationale**：
- FoldX 学术 license 需机构邮箱和 PI 信息，流程繁琐
- Rosetta `cartesian_ddg` 精度相当，完全免费
- 减少对外部审批的依赖，保证项目进度可控

### 决策 4：OpenMM GPU 后端 → OpenCL（放弃 Metal）

**Rationale**：
- conda-forge 上 OpenMM 的 Apple Metal 专用构建（`*apple*`）在 arm64 Python 3.13 环境下无法通过 `mamba install` 正常解析依赖
- 尝试直接安装 `openmm=8.5.1=py313h5024a7d_0_apple --no-deps` 后，Metal 后端仍不可用（仅列出 Reference/CPU/OpenCL）
- **结论**：conda-forge 的 Apple Silicon Metal 支持目前存在兼容性问题
- **替代方案**：OpenCL 后端正常工作，GPU 加速有效
- 性能差异：Metal 理论上比 OpenCL 快 10-20%，但在本项目尺度下不影响总体可行性

---

## 三、序列与突变映射

### 3.1 参考序列

| 蛋白 | UniProt | 长度 | 来源 |
|------|---------|------|------|
| 人 cGAS | Q8N884 | 522 aa | UniProt |
| 裸鼹鼠 cGAS | A0AAX6RS70 | 554 aa | UniProt (Heterocephalus glaber) |
| TRIM41 | Q8WV44 | 630 aa | UniProt |

**裸鼹鼠序列选择**：
- UniProt 上裸鼹鼠 cGAS 有多个注释版本（A0AAX6RS70 554 aa, A0AAX6RTF7 475 aa, A0AAX6RSF6 461 aa 等）
- 选择 **A0AAX6RS70**（554 aa）的原因：
  1. 长度最长，最可能包含完整功能域
  2. 论文中的突变位点编号（S463, E511, Y527, T530）在此序列上可直接对应
  3. 序列比对验证：保守区域 "NTGSYYEHVKI" 和 "GSPAVTLLI" 均存在

### 3.2 精确突变映射

通过全局序列比对（Biopython pairwise2，gap open=-10, gap extend=-0.5）确认：

| 论文标记 | 裸鼹鼠 aa | 人源对应位置 | 人源 aa | 突变方向 |
|---|---|---|---|---|
| S463 | S | **463** | C | C→S |
| E511 | E | **479** | K | K→E |
| Y527 | Y | **495** | L | L→Y |
| T530 | T | **498** | K | K→T |

**验证**：
- 人 cGAS 位置 463: C（在 `KDLGLCFDNCV` 上下文中）✅
- 人 cGAS 位置 479: K（在 `CLRTEKLENYF` 上下文中）✅
- 人 cGAS 位置 495: L（在 `LFSSNLIDKRS` 上下文中）✅
- 人 cGAS 位置 498: K（在 `SNLIDKRSKEF` 上下文中）✅

**重要说明**：论文中的编号是基于裸鼹鼠 cGAS 的坐标系统。由于裸鼹鼠序列比人源长 32 aa（主要在 N-端有延伸），导致位置 511→479、527→495、530→498 的偏移。这在论文 Methods 中需要明确说明。

---

## 四、硬件资源评估与性能基准

### 4.1 硬件规格

| 项目 | 规格 |
|------|------|
| CPU | Apple M3 Pro (12 核) |
| GPU | 18 核 Apple Silicon GPU (通过 OpenCL 访问) |
| 内存 | 36 GB 统一内存（无独立显存） |
| 操作系统 | macOS 15.3 (Darwin 25.3.0) |
| Python | 3.13.13 (conda-forge) |

### 4.2 OpenMM GPU 性能基准测试

#### 测试 1：小分子真空体系（~20k atoms, NoCutoff）
- 100k steps @ 2fs = 200 ps
- **OpenCL**: 0.2 ns in 0.7s → **23,667 ns/day**
- 该数字为理想上限，不代表真实 MD 性能

#### 测试 2：小分子显式溶剂 + PME（~60k atoms）
- 50k steps @ 4fs = 200 ps
- **OpenCL**: 200 ps in 127.6s → **135 ns/day**
- 包含 PME、溶剂、离子等全部真实 MD 开销

#### 测试 3：稳态 4ns benchmark（已终止）
- 1M steps @ 4fs = 4 ns，60k atoms，PME，HBonds 约束
- **超时终止**（900s limit），4ns 实际需 >15 分钟
- 推算速度：~< 384 ns/day（但此估算不精确，因未排除 Python overhead）
- **决定**：不再跑更长的 benchmark，之前 200ps 测试已足够可靠

**采用保守性能估计**（基于 Test 2）：
- 60k atoms, PME, HBonds, 4fs: **135 ns/day**
- 80k atoms (domain-truncated): **~95 ns/day** (scale factor 0.7)
- 200ns 轨迹: **~2.1 天**
- 12 trajectories serial: **~25 天**
- 12 trajectories (2 parallel): **~13 天**

### 4.3 体系规模与耗时估算

#### 情景 A：智能截取结构域（推荐方案）
- 假设保留 cGAS C-端域 (~250 aa) + TRIM41 C-端域 (~300 aa) = ~550 aa
- 加水后预计 **~80k atoms**
- 按 60k→135 ns/day 线性缩放（保守估计 0.7×）：
  - **预估速度：~95 ns/day**
  - 单轨迹 200 ns ≈ **2.1 天**
  - 4 体系 × 3 重复 = 12 轨迹
  - 串行：~25 天
  - **并行 2 个作业（内存允许）：~13 天**

#### 情景 B：全长复合物
- ~1150 aa，加水后 **~250k atoms**
- 按 60k→135 ns/day 线性缩放（250k/60k ≈ 4.2× 更慢）：
  - **预估速度：~32 ns/day**
  - 单轨迹 200 ns ≈ **6.3 天**
  - 12 轨迹串行 ≈ **75 天**，超出预算

**结论**：截取结构域是唯一能在"几周"内完成的方案。

### 4.4 内存限制
- 36 GB 统一内存需同时容纳体系数据 + 轨迹缓存 + 操作系统
- 单轨迹 200ns，60k atoms，每 10ps 一帧 = 20,000 帧
- DCD 格式：~3 bytes/atom/frame × 60k × 20k ≈ **3.6 GB/轨迹**
- 12 轨迹总存储：~43 GB（需及时归档到外部存储或压缩）
- 并行 2 个作业时内存占用：~2 × (体系 + 轨迹缓冲) ≈ 15-20 GB，安全

---

## 五、工作流程与时间线

### 5.1 总体时间线（~5 周）

```
Week 1: 环境搭建 ✅ | 结构预测（用户提交 AF3）⏳
Week 2: 体系构建 + 平衡化 + 第一批生产模拟启动
Week 3: 生产模拟（持续运行）+ 初步分析
Week 4: 生产模拟收尾 + 系统分析 + 突变扫描
Week 5: 对照验证 + 论文图表 + 文档撰写
```

### 5.2 详细阶段

#### 阶段 0：环境搭建（Day 1-2）✅
- Miniforge (arm64) 安装
- conda 环境 `cgas-md` 创建
- OpenMM 8.5.1 + AmberTools 24 + 分析库安装
- OpenCL GPU 后端验证

#### 阶段 1：结构预测（Day 3-5）⏳
- 工具：AlphaFold3 Server
- 4 组复合物：
  1. 人 cGAS_WT + TRIM41_WT
  2. 人 cGAS_4mut + TRIM41_WT
  3. 裸鼹鼠 cGAS_WT + TRIM41_WT
  4. 裸鼹鼠 cGAS_4mut_rev + TRIM41_WT
- 质量控制：ipTM > 0.70 接受，0.60-0.70 边缘，< 0.60 拒绝

#### 阶段 2：预分析与体系决策（Day 6-8）
- AF3 结构分析（界面大小、突变位点位置、pLDDT）
- 决定是否截取结构域
- Amber tleap 体系构建（ff19SB + OPC 水模型）
- 转换为 OpenMM XML

#### 阶段 3：MD 模拟（Day 9-30）
- 每个体系：能量最小化 → NVT 升温 → NPT 平衡（5ns）→ 生产（200ns）
- 3 个独立重复/体系，不同随机种子
- 参数：4fs HMR, LangevinMiddle, PME, MonteCarloBarostat

#### 阶段 4：结合能计算（Day 28-35）
- MM-GBSA（igb=5, saltcon=0.15M）
- 能量分解（per-residue）
- ΔΔG_bind = G(突变体) - G(WT)

#### 阶段 5：突变扫描（Day 33-38）
- Rosetta cartesian_ddg
- 4 个位点各自双向突变 + 组合突变
- Ala-scan 对照

#### 阶段 6：论文产出（Day 38-45）
- 图表：结构 overview、RMSD/RMSF、结合能箱线图、界面细节、突变扫描
- Methods 段落
- Results 初稿

---

## 六、风险预案

### 预案 A：AF3 预测置信度低（ipTM < 0.6）—— 已触发

**实际结果（job4_Hgal_4mut_rev + TRIM41）**：
- **ipTM = 0.17**（Model 0, best），其余 models 0.13-0.18
- **pTM = 0.36-0.37**
- **判定：❌ 拒绝**（远低于 0.60 最低标准）
- **fraction_disordered = 0.34-0.36**：约 1/3 的序列被预测为无序

**单体结构质量**：
- Chain A (cGAS Hgal 4mut_rev): mean pLDDT = 60.6, 51% 原子 ≥70
- Chain B (TRIM41): mean pLDDT = 62.2, 54% 原子 ≥70
- **单体折叠可接受，但相对取向完全不可信**

**科学解释**：
1. TRIM41 无任何实验结构
2. E3 连接酶-底物互作通常是 transient/weak（AF 的已知弱点）
3. 互作可能依赖翻译后修饰或 DNA 损伤信号
4. 论文本身未直接测量物理结合强度

**已启动应对方案 A1：蛋白-蛋白对接 + 对接后 MD**

| 决策 | 内容 |
|------|------|
| 对接工具 | **ClusPro** (https://cluspro.bu.edu/) |
| 选择原因 | 无需注册，纯网页，用户可手动提交；HADDOCK/LightDock 无法本地自动化 |
| 受体 | TRIM41 WT 全长 或 SPRY 域 (413-630) |
| 配体 | cGAS 全长 或 C-端域 (~200-554) |
| 活性残基限制 | cGAS 的 4 个突变位点必须位于界面 |
| 打分模式 | Attraction（酶-底物互作优化） |

**TRIM41 结构域关键发现**（UniProt Q8WV44）：
- **B30.2/SPRY domain: residues 413-630** — 这是底物识别结构域！
- RING finger (N-端): E2 结合
- 多个无序区域: 51-86, 143-176, 503-538
- **这与我们的截取策略完全吻合**：SPRY 域负责识别 cGAS

**已准备文件**：
- `cgas_fixed.pdb` / `trim41_fixed.pdb`（全长，加氢补全）
- `cgas_CT_200-554.pdb`（355 residues）/ `trim41_SPRY_413-630.pdb`（218 residues）

**截断后残基编号的关键验证**：
- 截断 PDB 保留原始编号（cGAS CTD 从 200 开始，TRIM41 SPRY 从 413 开始）
- 突变位点 463, 511, 527, 530 在截断 PDB 中编号**不变**
- ClusPro Active Residues 可直接使用相同数字，无需转换

**对接后流程**：
1. ClusPro 生成 top 10-50 poses
2. RMSD 聚类去冗余
3. 每个 pose：OpenMM minimization → 5-10 ns 平衡 MD
4. 稳定性筛选（RMSD 收敛）
5. 最优 1-3 poses：200 ns 生产 MD + MM-GBSA

**备选**：若 ClusPro 也失败 → 方案 A2（单体 MD + 表面分析）

**单体 MD 同步准备**：
- 即使对接成功，仍需单体 MD 作为对照
- 可检测 allosteric 效应：4 个突变是否改变了远端区域动态？
- 比较人 WT vs 4mut 的表面静电势差异

### 预案 B：MD 轨迹不收敛（RMSD 持续漂移）
- **应对**：
  1. 延长截取的柔性 linker
  2. aMD/GaMD 增强采样
  3. 缩短分析窗口（只用后 100ns）

### 预案 C：结合能趋势与实验相反
- **即**：MM-GBSA 预测裸鼹鼠结合更强，但实验显示泛素化更弱
- **科学解释**（可写入论文）：
  - 物理结合 ≠ 泛素化效率：TRIM41 可能仍能结合，但 E2-Ub 传递几何不利
  - 结合姿态改变：突变未削弱结合，但改变了 Lys 位点相对于 RING 的空间取向
  - MD 可检测 "same binding, different geometry"

### 预案 D：OpenCL 性能低于预期
- **应对**：
  - 缩短轨迹至 150ns
  - 减少重复至 2 个/体系
  - 进一步截取结构域

---

## 七、文件清单

### 已创建文件

```
sequences/
  cgas_trim41_sequences.fasta      # 合并 FASTA
  Hsap_cGAS_WT.fasta
  Hsap_cGAS_4mut.fasta
  Hgal_cGAS_WT.fasta
  Hgal_cGAS_4mut_rev.fasta
  TRIM41_WT.fasta
  AF3_server_submissions.md         # AF3 提交指南

scripts/
  analyze_af3.py                    # AF3 结构质量分析
  build_system.py                   # MD 体系构建（Amber/OpenMM）
  run_md.py                         # 生产 MD 模拟
  analyze_trajectory.py             # 轨迹分析（RMSD/RMSF/氢键/接触）
  run_mmpbsa.py                     # MM-GBSA/PBSA 结合能

structures/af3_raw/                 # AF3 结果存放位置
  job1_Hsap_WT/
  job2_Hsap_4mut/
  job3_Hgal_WT/
  job4_Hgal_4mut_rev/

data/
  md_runs/                          # MD 系统文件和轨迹
  analysis/                         # 分析输出
  rosetta/                          # Rosetta 突变扫描

docs/
  computational_workflow.md          # 原始方案设计
  paper_notes_cgas_trim41.md         # 论文理解笔记
  execution_plan_v1.md               # 执行方案 v1.0
  project_log.md                     # 本文档
  cluspro_submission_guide.md        # ClusPro 对接提交指南

structures/af3_raw/job4_Hgal_4mut_rev/
  ranked_0_chain_A.pdb               # cGAS 单体 (AF3 Model 0)
  ranked_0_chain_B.pdb               # TRIM41 单体 (AF3 Model 0)
  cgas_fixed.pdb                     # cGAS 单体 (PDBFixer 修复)
  trim41_fixed.pdb                   # TRIM41 单体 (PDBFixer 修复)
  cgas_CT_200-554.pdb                # cGAS C-端域截取 (355 residues)
  trim41_SPRY_413-630.pdb            # TRIM41 SPRY 域截取 (218 residues)
  plddt_profile.png                  # pLDDT 分布图
  confidence_scores.png              # 5 个 model 的 confidence 柱状图

README.md                            # 项目概览与快速开始
```

---

## 八、待补充数据

- [x] AF3 预测结果（job4 已完成，ipTM=0.17 极低，单体可接受）
- [ ] AF3 预测结果（job1-3 等待中）
- [ ] 稳态 OpenCL benchmark 结果（4ns 测试，后台运行中）
- [ ] 结构域截取决策（基于 AF3 单体分析：TRIM41 SPRY 413-630, cGAS CTD ~200-554）
- [ ] ClusPro 对接结果（top poses, 聚类分析）
- [ ] 对接 pose 稳定性筛选（短 MD 结果）
- [ ] 实际 MD 速度（首个体系首次运行后更新）
- [ ] MM-GBSA ΔG 和 ΔΔG 数值
- [ ] Rosetta ddG 突变扫描结果

---

## 九、软件版本记录

| 软件 | 版本 | 来源 |
|------|------|------|
| Python | 3.13.13 | conda-forge |
| OpenMM | 8.5.1 | conda-forge (apple build) |
| AmberTools | 24.8 | conda-forge |
| MDAnalysis | 2.10.0 | conda-forge |
| MDTraj | 1.11.1 | conda-forge |
| Biopython | 1.87 | conda-forge |
| NumPy | 2.4.3 | conda-forge |
| SciPy | 1.17.1 | conda-forge |
| Matplotlib | 3.10.8 | conda-forge |
| Seaborn | 0.13.2 | conda-forge |
| Pandas | 2.3.3 | conda-forge |

---

*文档创建：2026-04-22*
*维护者：Kimi Code CLI*
