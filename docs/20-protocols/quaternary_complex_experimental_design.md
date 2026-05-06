# E2~Ub-TRIM41-cGAS 四元复合物：详细实验设计方案

> 目标：构建一个包含 E2~Ub（UBE2D1-Ub 硫酯模拟物）、TRIM41 RING 二聚体、TRIM41 SPRY 域和 cGAS 的完整四元催化复合物模型，通过 MD 模拟评估 WT vs 4mut 在催化几何上的差异。

---

## 1. 设计原则与可行性分析

### 1.1 为什么这个模型是可实现的

| 组件 | 是否有实验结构 | 置信度 | 获取方式 |
|:---|:---|:---|:---|
| E2~Ub (UbcH5b~Ub) | ✅ 有多个晶体结构 | 高 | PDB 直接下载 |
| TRIM RING 二聚体 | ✅ TRIM25/32/2/5α 均有 | 高 | 同源建模 + 二聚体移植 |
| TRIM41 SPRY | ⚠️ 无实验结构 | 中 | AF3 预测 + 我们已有 docking pose |
| cGAS | ✅ AF3 预测 (pTM 0.87) | 高 | AF3 单体 |
| cGAS-SPRY 复合物 | ⚠️ 无实验结构 | 中 | 我们的 Rosetta/LightDock docking |

**关键突破点**：Zhen et al. (2014, *PLoS One*) 已经证明，可以通过**结构叠合**的方法，将 E3-E2-Ub 三元结构与底物-E2-E3 结构拼接，构建出 E3-E2-Ub-底物 四元复合物。他们的模型经过了 35 ns MD 平衡和 QM/MM 验证，是可靠的参考。

### 1.2 核心建模策略

采用 **"模块化拼装 + MD 平衡"** 策略：

```
Step 1: 获取/构建各模块
  ├── E2~Ub 模块：从 PDB 5FER (TRIM25 RING + UBE2D1-Ub) 提取
  ├── TRIM41 RING 二聚体：同源建模（以 TRIM25/2 RING 为模板）
  ├── cGAS-TRIM41(SPRY) 模块：使用我们现有的 docking pose
  └── 连接区：coiled-coil + B-box（简化处理，见下文）

Step 2: 拼装四元复合物
  ├── 将 TRIM41 RING 叠合到 TRIM25 RING 上
  ├── 保留 E2~Ub 相对于 RING 的位置
  ├── 将 cGAS-SPRY 复合物整体定位到合适位置
  └── 手动调整 / MD 平衡使各组件达到合理几何

Step 3: MD 模拟
  ├── 多步能量最小化
  ├── 渐进式加热和平衡
  └── 生产运行（3 reps × 200-500 ns）

Step 4: 后处理分析
  ├── 催化几何指标
  ├── 动态/变构指标
  └── WT vs 4mut 统计比较
```

---

## 2. 阶段一：结构获取与准备（预计 3-5 天）

### 2.1 E2~Ub 模块

**选择 E2 亚型**：

TRIM 家族通常与 **UBE2D (UbcH5)** 家族互作。文献证据：
- TRIM25 晶体结构 (PDB: 5FER) 使用 **UBE2D1**
- TRIM2 晶体结构 (PDB: 7ZJ3) 使用 **UBE2D1~Ub 硫酯模拟物**
- TRIM32 使用 UbcH5 (Albor et al., 2006)
- TRIM58 使用 UbcH5B (Xu et al., 2003)

**推荐选择：UBE2D1 (UbcH5a) 或 UBE2D2 (UbcH5b)**，两者序列高度相似 (>90%)，且均有丰富的晶体结构数据。

**硫酯键建模**：

真实 E2~Ub 硫酯键在 MD 力场中无法直接模拟（共价键断裂/形成需要 QM/MM）。有两种策略：

| 策略 | 方法 | 优点 | 缺点 |
|:---|:---|:---|:---|
| **A. Isopeptide mimic** | 用 K85-Ub 异肽键代替 C85-Ub 硫酯键 | 力场支持，可直接模拟 | 键的化学性质略有不同（酰胺 vs 硫酯） |
| **B. Cysteine thioester** | K85→C85 突变，手动构建硫酯键拓扑 | 更接近真实化学 | 需要自定义力场参数，较复杂 |
| **C. 共价约束** | 用 harmonic restraint 约束 C85-S 与 Ub-G76-C 的距离 | 灵活性高 | 人为约束可能影响动力学 |

**推荐策略 A（isopeptide mimic）**，因为：
1. TRIM25 晶体结构 (5FER) 本身就使用了 isopeptide-mimic（K85 与 Ub C-末端通过异肽键连接）
2. 文献中广泛使用这种策略（Pruneda 2012; Plechanovova 2011）
3. 对于构象采样目的，isopeptide 足以捕捉 E2~Ub 的 closed/open 几何

**具体操作**：
```bash
# 从 PDB 5FER 下载 TRIM25 RING + UBE2D1-Ub 复合物
# 提取 chain B (UBE2D1) 和 chain C/D (Ub)
# 注意：5FER 是 TRIM25 RING 二聚体 + 两个 UBE2D1-Ub 分子
# 我们只需要一个 RING + 一个 E2~Ub
```

**备用结构**：
- PDB 7ZJ3 (TRIM2 RING + UBE2D1~Ub) — 分辨率 2.53 Å
- PDB 4AP4 (RNF4 RING + UbcH5A-Ub) — 经典 RING-E2~Ub 结构
- PDB 7R71 (Ark2C RING + UbcH5b~Ub) — 包含 regulatory Ub

### 2.2 TRIM41 RING 域

**序列信息**：
- 人类 TRIM41 (UniProt: Q8WV44)
- RING 域：约 residues 1-60（具体边界需用 CD-search / Pfam 确认）

**同源建模模板选择**：

| 模板 | PDB | 分辨率 | 与 TRIM41 RING 序列一致性 | 备注 |
|:---|:---|:---|:---|:---|
| TRIM25 RING | 5FER | 2.6 Å | ~40-50% | 有 E2~Ub 共晶 |
| TRIM2 RING | 7ZJ3 | 2.53 Å | ~40-50% | 有 E2~Ub 共晶 |
| TRIM32 RING | 5FEY | ~2.5 Å | ~40-50% | RING 二聚体 |
| TRIM5α RING | 4TKP | ~2.5 Å | ~35-45% | 高阶寡聚体 |

**建模流程**：

```
1. 用 ClustalOmega / MAFFT 做 TRIM41 RING 与模板的多序列比对
2. 用 Modeller 或 SWISS-MODEL 构建 TRIM41 RING 单体结构
3. 用 TRIM25/2/32 的 RING 二聚体结构指导 TRIM41 RING 二聚化：
   - 将两个 TRIM41 RING 单体叠合到模板的两个单体上
   - 保留模板中的 RING-RING 二聚化界面
4. 用 GROMACS / Amber 进行短 MD 平衡（10-20 ns），优化二聚体界面
```

**关键考量**：
- TRIM RING 二聚化由两个短螺旋介导（Yudina et al., 2015; Kiss et al., 2019）
- 不同 TRIM 的 RING 二聚体相对取向有差异（Fig 3B in TRIM25 论文）
- TRIM41 的 RING 可能类似于 TRIM25（低亲和力，需要 E2~Ub 稳定）或 TRIM32（组成型二聚）
- **保守策略**：用 TRIM25 的 RING 二聚体构象作为起点（因为有 E2~Ub 共晶），然后在 MD 中允许其 relax

### 2.3 cGAS-TRIM41(SPRY) 复合物

**直接使用我们现有的模型**：
- Hsap_WT + TRIM41 SPRY 的 docking pose
- Hsap_4mut + TRIM41 SPRY 的 docking pose

**需要做的调整**：
- 确认 docking pose 中 cGAS 的 **K315** 位置（这是预期的泛素化位点）
- 如果 K315 被埋在 interface 内部，需要检查是否合理（Chen et al. 的实验表明 K315 是泛素化位点，因此它应该在可接近的位置）

### 2.4 连接区（Coiled-Coil + B-Box）

**这是最大的不确定性来源**。

TRIM41 的完整 N-端区域为：RING → B-box1 → B-box2 → Coiled-Coil → SPRY

**三种处理策略**：

| 策略 | 描述 | 复杂度 | 推荐度 |
|:---|:---|:---|:---|
| **A. 完全省略** | 只模拟 RING + cGAS-SPRY，用柔性 linker 连接 | 低 | ⭐⭐⭐ |
| **B. 刚性连接** | 用短 linker（如 10-20 个 Gly/Ser）连接 RING 和 SPRY | 中 | ⭐⭐⭐⭐ |
| **C. 完整建模** | 用 AlphaFold3 预测 RING-BB-CC-SPRY 全长，然后嵌入 E2~Ub 和 cGAS | 高 | ⭐⭐⭐⭐⭐ |

**推荐策略 C 的简化版**：

1. 用 AlphaFold3 预测 **TRIM41(RING-BB-CC-SPRY)** 单体结构
2. 用 AF3 的 **multimer** 功能预测 **TRIM41(RING-BB-CC-SPRY) + cGAS** 复合物
3. 将 E2~Ub 从 PDB 5FER 移植到 AF3 预测的 RING 域上
4. 用 MD 平衡整个系统

**但 AF3 对 TRIM41 的预测置信度可能很低**（之前的多聚体预测 ipTM < 0.25）。

**折中方案（推荐）**：

```
1. 用 AlphaFold3 预测 TRIM41(RING-BB-CC) — 这是 N-端三聚体模序
2. 用 AlphaFold3 预测 TRIM41(SPRY) — 单独预测置信度更高
3. 用我们现有的 cGAS-SPRY docking pose
4. 用以下方式连接：
   - RING-BB-CC 的 C-端 ↔ SPRY 的 N-端
   - 中间用一段柔性的 Gly/Ser-rich linker (~20 residues)
   - 或者干脆将 RING-BB-CC 和 SPRY 作为两个独立的 rigid body，
     用弱约束（harmonic restraints）限制它们的相对位置范围
```

**更务实的方案（强烈推荐）**：

> 不要试图构建"完美"的全长 TRIM41。相反，构建一个 **"功能核心"模型**：
> - **核心模块**：RING 二聚体 + E2~Ub（从 PDB 5FER）
> - **底物模块**：cGAS-TRIM41(SPRY)（从我们现有的 docking）
> - **连接**：将两个模块放置在合理的相对位置上，使 cGAS K315 面向 E2~Ub 催化中心
> - **在 MD 中允许 RING 和 SPRY 之间的相对运动**

这种"功能核心"模型的优势：
- 避免了不确定的连接区建模
- 聚焦于催化几何（K315 到 Ub-G76 的距离）
- 计算成本低（原子数更少）
- 与文献中的建模策略一致（Zhen 2014 用的也是类似方法）

### 2.5 四元复合物拼装

**参照 Zhen et al. (2014) 的方法**：

```python
# 伪代码
# Step 1: 加载模板结构
trim25_ring_dimer = load("5FER")  # 提取 chain A,B (RING 二聚体)
e2_ub = load("5FER")              # 提取 chain E,F (UBE2D1-Ub)

# Step 2: 构建 TRIM41 RING 二聚体
trim41_ring = homology_model(trim41_ring_seq, template=trim25_ring_dimer)
trim41_ring_dimer = superpose(trim41_ring, trim25_ring_dimer)

# Step 3: 将 E2~Ub 移植到 TRIM41 RING 上
# E2~Ub 相对于 RING 的位置通过叠合保留
e2_ub_on_trim41 = superpose(e2_ub, trim25_ring_dimer, 
                             align_atoms="RING interface residues")

# Step 4: 加载 cGAS-SPRY 复合物
cgas_spri = load("our_docking_pose.pdb")

# Step 5: 定位 cGAS-SPRY
# 关键：使 cGAS K315 的侧链指向 E2~Ub 催化中心
# 可以通过手动平移/旋转，或基于以下约束自动对接：
#   - SPRY 的 N-端与 RING 的 C-端距离在合理范围内
#   - K315 NZ 到 Ub-G76 C 距离 < 15 Å（初始 guess）

# Step 6: 能量最小化
minimize(quaternary_complex)
```

**更自动化的方案：用 HADDOCK 或 LightDock 做数据驱动的对接**：

```
约束定义：
1. Active residues (必须接触)：
   - RING 上的 E2 结合位点残基（从 5FER 映射）
   - SPRY 上的 cGAS 结合位点残基（从我们的 docking 映射）

2. Passive residues (可接触)：
   - RING 和 SPRY 的表面残基

3. 模糊约束 (ambiguous restraints)：
   - RING C-端 ↔ SPRY N-端（距离 10-30 Å）
   - cGAS K315 ↔ Ub-G76（距离 5-15 Å）
```

**推荐工具链**：
- **PyMOL / ChimeraX**：手动可视化调整
- **HADDOCK2.4/2.5**：数据驱动的蛋白-蛋白对接（支持多种约束）
- **Rosetta docking_protocol**：如果 HADDOCK 效果不佳

---

## 3. 阶段二：MD 模拟（预计 2-3 周）

### 3.1 系统构建

**预估系统规模**：

| 组件 | 残基数 | 原子数（含水/离子）|
|:---|:---|:---|
| TRIM41 RING 二聚体 | ~120 | ~2,000 |
| UBE2D1~Ub | ~230 | ~3,500 |
| cGAS (522 aa) | 522 | ~8,000 |
| TRIM41 SPRY | ~200 | ~3,000 |
| **蛋白总计** | ~1,070 | **~16,500** |
| 溶剂 + 离子 (10 Å buffer) | — | **~45,000** |
| **系统总计** | — | **~62,000** |

这比我们的现有系统（85,510 原子）**更小**！因为 SPRY 域（~200 aa）比完整的 TRIM41 SPRY 域 + 大量溶剂化水要小。

实际上，如果我们用完整系统（包含连接区），原子数可能在 80,000-100,000 左右，仍然可管理。

**力场选择**：
- **蛋白**：Amber ff19SB（与现有模拟一致）
- **水**：OPC（与现有模拟一致）或 TIP3P
- **E2~Ub 异肽键**：ff19SB 天然支持（lysine side chain - C-terminal carboxylate isopeptide bond）

**溶剂化**：
- 溶剂盒：truncated octahedron
- Buffer：12 Å（比现有 10 Å 稍大，因为复合物更延展）
- 离子：Na⁺/Cl⁻，0.15 M，中性化

### 3.2 模拟流程

```
Step 1: 能量最小化
  ├── 5000 steps，蛋白重原子约束 (500 kJ/mol/nm²)
  ├── 5000 steps，蛋白 CA 约束 (100 kJ/mol/nm²)
  └── 10000 steps，无约束

Step 2: 加热 (NVT)
  ├── 0 → 100 K，50 ps，蛋白 CA 约束
  ├── 100 → 200 K，50 ps，蛋白 CA 约束
  └── 200 → 300 K，100 ps，蛋白 CA 约束

Step 3: 密度平衡 (NPT)
  ├── 200 ps，蛋白骨架约束
  └── 200 ps，蛋白 CA 约束

Step 4: 预平衡 (NPT)
  └── 1 ns，无约束

Step 5: 生产运行 (NPT)
  ├── 3 independent replicas
  ├── 每个 replica：200-500 ns
  ├── 保存间隔：每 10 ps 一帧（20,000-50,000 帧/rep）
  └── 总计：600 ns - 1.5 μs
```

**模拟参数**（与现有模拟一致）：
```python
# OpenMM 参数示例
integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtoseconds)
system.addForce(MonteCarloBarostat(1.0*bar, 300*kelvin))
nonbonded_method = PME
nonbonded_cutoff = 1.0*nanometer
constraints = HBonds
```

### 3.3 计算资源估算

| 系统 | 原子数 | 单 rep 速度 (4× RTX 3090) | 500 ns 用时 |
|:---|:---|:---|:---|
| 四元复合物 | ~80,000 | ~15-20 ns/day | ~25-35 天/rep |
| 3 reps | — | 可并行（3 GPUs） | **~25-35 天 total** |

**注意**：
- 如果系统更大（包含完整连接区，>100,000 原子），速度可能降至 ~10 ns/day
- 可以考虑 **protein-only 模拟**（隐式溶剂）做快速初筛，但 explicit solvent 更可靠
- 或者先用 **50 ns 短模拟** 测试模型稳定性，确认合理后再跑长模拟

---

## 4. 阶段三：后处理分析（预计 1 周）

### 4.1 催化几何指标（核心）

#### 指标 1：K315 NZ 到 Ub-G76 C 距离

```python
# MDAnalysis 代码示例
k315_nz = u.select_atoms("resid 315 and name NZ")  # cGAS K315
g76_c = u.select_atoms("resid XXX and name C")      # Ub G76 C-末端羧基碳
dist = np.linalg.norm(k315_nz.positions - g76_c.positions, axis=1)
```

- **< 5 Å**：可以直接发生亲核攻击（理想催化距离）
- **5-8 Å**：需要小幅构象调整
- **8-15 Å**：需要显著的构象变化
- **> 15 Å**：不太可能发生催化

#### 指标 2：K315 NZ 到 E2-C85 Sγ 距离

反映底物赖氨酸到硫酯键催化中心的距离。

- 与指标 1 联合使用，评估整体催化几何

#### 指标 3：E2~Ub "Closed" 构象分数

基于 Pruneda et al. (2012) 和 Dou et al. (2012) 的定义：

```python
# Closed 构象的关键特征：
# 1. Ub Ile44 CA 到 E2 α2 helix 的 COM 距离 < 15 Å
# 2. Ub Gly76 C 到 E2 Cys85 Sγ 距离 < 8 Å
# 3. Ub 的 C-末端 tail 指向 E2 活性位点

ile44 = u.select_atoms("resid XXX and resname UB and name CA")  # Ub I44
e2_a2 = u.select_atoms("resid XXX-XXX and name CA")             # E2 α2 helix
g76 = u.select_atoms("resid XXX and resname UB and name C")     # Ub G76
c85 = u.select_atoms("resid XXX and name SG")                   # E2 C85

dist_i44_a2 = distance_array(ile44.positions, e2_a2.center_of_mass())
dist_g76_c85 = distance_array(g76.positions, c85.positions)

is_closed = (dist_i44_a2 < 15) & (dist_g76_c85 < 8)
closed_fraction = np.mean(is_closed)
```

#### 指标 4：RING-RING 二聚化界面稳定性

```python
ring1 = u.select_atoms("chain A and resid 1-60")
ring2 = u.select_atoms("chain B and resid 1-60")
ring_com_dist = distance_array(ring1.center_of_mass(), ring2.center_of_mass())
```

- TRIM RING 二聚体界面如果被破坏，E2~Ub 的 closed 构象可能无法稳定

#### 指标 5：底物赖氨酸的溶剂可及性（SASA）

```python
# 计算 K315 侧链的 SASA
k315_sidechain = u.select_atoms("resid 315 and (name NZ or name CE or name CD)")
sasa = calc_sasa(k315_sidechain)
```

- 如果 K315 被埋在 cGAS-TRIM41 interface 内部，它无法被泛素化
- 4mut 可能通过改变 cGAS 构象来暴露 K315

### 4.2 动态/变构指标

| 指标 | 方法 | 用途 |
|:---|:---|:---|
| **RMSF** | MDAnalysis `RMSF` | 识别各组件的柔性区域 |
| **PCA / Essential Dynamics** | MDAnalysis `PCA` 或 ` encore` | 提取主导运动模式，比较 WT vs 4mut 的差异 |
| **DCCM** | 自写脚本或 `mdentropy` | 识别 allosteric 耦合路径 |
| **H-bond 占有率** | MDAnalysis `HydrogenBondAnalysis` | 界面关键相互作用的稳定性 |
| **Salt bridge 分析** | MDAnalysis | 静电相互作用的持久性 |
| **Contact map** | MDAnalysis `contacts.Contacts` | 残基接触频率的变化 |

### 4.3 构象集合分析

#### 自由能面（FEL）

```python
# 使用两个主成分作为 CV
from sklearn.decomposition import PCA

# 收集所有 reps 的 CA 坐标
coords = ...  # shape: (n_frames, n_ca_atoms * 3)
pca = PCA(n_components=2)
pc_proj = pca.fit_transform(coords)

# 计算自由能面
import numpy as np
H, xedges, yedges = np.histogram2d(pc_proj[:,0], pc_proj[:,1], bins=50)
F = -kB * T * np.log(H + 1e-10)

# 可视化 WT vs 4mut 的 FEL
```

#### 聚类分析

```python
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans

# 基于催化几何特征聚类
features = np.column_stack([
    k315_g76_dist,
    k315_c85_dist,
    closed_score,
    ring_com_dist,
])

gmm = GaussianMixture(n_components=4, covariance_type='full')
labels = gmm.fit_predict(features)
```

### 4.4 统计比较

对每个指标，比较 WT vs 4mut：

```python
from scipy import stats

# 合并 3 reps 的数据
wt_data = np.concatenate([wt_rep1, wt_rep2, wt_rep3])
mut_data = np.concatenate([mut_rep1, mut_rep2, mut_rep3])

# Welch t-test
t_stat, p_value = stats.ttest_ind(wt_data, mut_data, equal_var=False)

# 效应量（Cohen's d）
pooled_std = np.sqrt((np.std(wt_data)**2 + np.std(mut_data)**2) / 2)
cohens_d = (np.mean(wt_data) - np.mean(mut_data)) / pooled_std
```

---

## 5. 预期结果与可检验假设

### 5.1 核心假设

> **H1**: 4mut 不改变 cGAS-TRIM41 的结合亲和力（MM-GBSA p=0.5，已验证），但改变 cGAS 的构象集合，使 K315 在 E2~Ub-TRIM41 存在时更频繁地出现在催化距离内。

### 5.2 预期观察

| 预测 | WT | 4mut | 生物学意义 |
|:---|:---|:---|:---|
| K315-UbG76 平均距离 | ~12 Å | ~7 Å | 4mut 使底物赖氨酸更接近催化中心 |
| Closed E2~Ub 构象分数 | ~30% | ~60% | 4mut 通过 allosteric 效应稳定 closed 构象 |
| K315 SASA | 低（埋藏） | 高（暴露） | 4mut 暴露泛素化位点 |
| RING 二聚体稳定性 | 相似 | 相似 | 4mut 不影响 RING 域本身 |
| cGAS N-端 RMSF | 较低 | 较高 | 4mut 增加 N-端柔性，有利于采样催化构象 |

### 5.3 阴性对照

如果结果不支持假设：
- WT 和 4mut 的 K315-UbG76 距离无显著差异 → 4mut 的作用机制可能不涉及催化几何
- 4mut 的 closed E2~Ub 分数反而更低 → 4mut 可能通过其他机制（如改变 E2 招募）发挥作用

---

## 6. 风险与应对

| 风险 | 可能性 | 影响 | 应对措施 |
|:---|:---|:---|:---|
| TRIM41 RING 同源建模不准确 | 中 | 高 | 用多个模板（TRIM25/2/32）分别建模，比较结果一致性 |
| 四元复合物在 MD 中解离 | 高 | 高 | 初始阶段用弱约束（flat-bottom restraints）维持关键距离，逐步释放 |
| 系统太大导致模拟太慢 | 中 | 中 | 先用 "功能核心" 模型（省略连接区），验证概念后再扩展 |
| E2~Ub 的 isopeptide mimic 不准确 | 低 | 中 | 与文献对比，确认 closed/open 构象的关键几何特征是否被保留 |
| K315 不在预期的催化位置 | 中 | 高 | 检查 Chen et al. 的实验数据，确认 K315 是主要泛素化位点；如果不是，换用实验验证的位点 |
| AF3 预测失败 | 中 | 中 | 有 PDB 模板可用，不依赖 AF3 |

---

## 7. 时间线估算

| 阶段 | 任务 | 预计时间 |
|:---|:---|:---|
| **Week 1** | 结构准备：TRIM41 RING 同源建模、E2~Ub 提取、cGAS-SPRY 复用 | 3-5 天 |
| **Week 1-2** | 四元复合物拼装：HADDOCK/Rosetta 对接 + 手动调整 + 最小化测试 | 5-7 天 |
| **Week 2-3** | 短 MD 测试：50 ns × 1 rep，检查稳定性 | 3-5 天 |
| **Week 3-6** | 生产 MD：3 reps × 200-500 ns | 3-4 周 |
| **Week 6-7** | 后处理分析：催化几何、动态分析、统计比较 | 1 周 |
| **总计** | | **~6-8 周** |

**与现有 Hgal MD 并行**：Hgal 的 3 reps × 200 ns 正在跑（预计 2-3 周完成），可以在这期间完成四元模型的构建和测试，待 Hgal 完成后切换 GPU 资源到四元复合物模拟。

---

## 8. 工具与脚本清单

### 8.1 结构准备

```bash
# 1. 下载 PDB
wget https://files.rcsb.org/download/5FER.pdb
wget https://files.rcsb.org/download/7ZJ3.pdb
wget https://files.rcsb.org/download/4AP4.pdb

# 2. 序列比对 (ClustalOmega)
clustalo -i trim41_ring.fa -o alignment.fa --templatefile=5FER_A.pdb

# 3. 同源建模 (Modeller)
python model_trim41_ring.py

# 4. 结构叠合与拼装 (PyMOL / ChimeraX)
# 手动或使用脚本

# 5. 加氢、能量最小化 (pdb4amber, tleap)
pdb4amber -i complex.pdb -o complex_H.pdb -y
```

### 8.2 MD 模拟

```bash
# OpenMM 脚本（基于现有 run_production.py 修改）
python run_quaternary_md.py \
  --system Hsap_WT_quaternary \
  --prmtop complex.prmtop \
  --inpcrd complex.inpcrd \
  --nsteps 250000000 \
  --output rep1/
```

### 8.3 后处理分析

```python
# 分析脚本（基于现有 analyze_s305e.py 扩展）
scripts/03_analysis/analyze_quaternary.py
# - 催化几何指标
# - PCA / FEL
# - 统计比较
```

---

## 9. 与现有工作的衔接

### 9.1 复用现有数据

| 现有资源 | 复用方式 |
|:---|:---|
| cGAS WT/4mut AF3 结构 | 直接作为四元复合物的 cGAS 组件 |
| cGAS-TRIM41(SPRY) docking pose | 作为四元复合物的底物模块定位参考 |
| MD 模拟流程 (OpenMM) | 直接使用，仅需修改拓扑文件 |
| 分析脚本 | 扩展 `analyze_s305e.py` 添加四元复合物专用指标 |

### 9.2 新增工作

| 新增内容 | 工作量 |
|:---|:---|
| TRIM41 RING 同源建模 | 中等 |
| E2~Ub 硫酯模拟物建模 | 小（已有模板） |
| 四元复合物拼装 | 较大（需要多次迭代优化） |
| 系统构建（拓扑/坐标） | 小（用 tleap/pdb4amber） |
| 四元复合物 MD（3 reps × 200-500 ns） | 大（3-4 周 GPU 时间） |
| 催化几何分析脚本 | 中等 |

---

## 10. 最终建议

### 推荐执行顺序

```
Phase 1 (现在就可以开始):
  ├── TRIM41 RING 序列分析 + 同源建模
  ├── 下载并处理 PDB 5FER/7ZJ3
  └── 用 PyMOL 手动尝试初步拼装

Phase 2 (Hgal MD 完成前):
  ├── 用 HADDOCK 或 Rosetta 优化四元复合物拼装
  ├── 50 ns 测试 MD（检查稳定性）
  └── 根据测试结果调整模型

Phase 3 (Hgal MD 完成后):
  ├── 3 reps × 200-500 ns 生产 MD
  └── 催化几何 + 动态分析

Phase 4 (论文修订):
  ├── 将四元模型结果融入 Discussion
  └── 作为 "未来方向" 或 "验证假设"
```

### 最小可行版本 (MVP)

如果完整四元复合物太复杂，可以先做一个 **"概念验证" 版本**：

1. 只做 **E2~Ub + cGAS**（省略 TRIM41 RING 和 SPRY）
2. 将 E2~Ub 放置在 cGAS K315 附近（基于文献中的典型催化距离）
3. 做短 MD（50-100 ns）
4. 比较 WT vs 4mut 的 K315 到 Ub-G76 距离分布

这个 MVP 虽然缺少 E3 的变构调控细节，但可以**直接验证 4mut 是否改变了 cGAS 的催化几何**——这是论文中最核心的主张。

---

*文档版本：v1.0*
*日期：2026-05-03*
