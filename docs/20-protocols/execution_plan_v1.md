# cGAS-TRIM41 计算研究执行方案 v1.0

> 目标：以可发表论文标准，通过分子模拟解释裸鼹鼠 cGAS 4 个氨基酸变异对 TRIM41 互作的影响。
> 
> 计算资源：Apple M3 Pro (36 GB 统一内存, 12 核), 可用时间 ~4-5 周。

---

## 一、资源现实评估

### 1.1 硬件边界
| 项目 | 现状 | 影响 |
|------|------|------|
| Apple M3 Pro GPU | Metal 后端，18 核 GPU | OpenMM 支持良好，但性能约为同代 RTX 的 30-50% |
| 36 GB 统一内存 | 无独立显存 | 体系原子数上限约 25-30 万（全长复合物接近上限） |
| 可用时间 | 3-5 周 | 无法完成 μs 级全长复合物 × 多重复的常规 MD |

### 1.2 体系规模估算
- 人 cGAS (522 aa) + TRIM41 (630 aa) ≈ 1150 aa
- 显式水 + 离子 ≈ **25-35 万原子**
- 在 M3 Pro Metal 上，25 万原子体系预估速度：**30-60 ns/day**（4 fs HMR）
- 即 **100 ns 需 2-3 天**，500 ns 需 10-15 天

**结论**：若做全长复合物的常规 MD，4 个体系 × 3 重复 × 200 ns = 2400 ns 总轨迹，需 **50-80 天**，超出预算。

### 1.3 可行路径
采取 **"结构预测 → 结构域聚焦 → 短而多重复 MD + 突变扫描"** 策略，这是目前资源下唯一能产出 robust 统计数据的方案。

---

## 二、总体时间线（~5 周）

```
Week 1: 环境搭建 + 结构预测（需您协助 AF3）
Week 2: 体系构建 + 平衡化 + 第一批生产模拟启动
Week 3: 生产模拟（持续运行）+ 初步分析
Week 4: 生产模拟收尾 + 系统分析 + 突变扫描
Week 5: 对照验证 + 论文图表 + 文档撰写
```

---

## 三、分阶段详细方案

### 阶段 0：环境搭建与序列准备（Day 1-2）

**我负责执行：**
1. 安装 Miniforge (arm64 conda-forge)
2. 创建 `cgas-md` conda 环境：
   - `openmm` (Metal 后端)
   - `ambertools` (tleap, cpptraj, MMPBSA.py)
   - `mdanalysis`, `mdtraj`
   - `pdbfixer`
   - `matplotlib`, `seaborn`, `pandas`
   - `biopython`
3. 精确确认序列与突变映射：
   - 人 cGAS: `Q8N884` (522 aa)
   - 裸鼹鼠 cGAS: 确定最佳参考序列（目前候选 `A0AAX6RS70`, 554 aa）
   - TRIM41: `Q8WV44` (630 aa)
   - **关键任务**：将论文中的 4 个位点（S463, E511, Y527, T530）精确映射到上述 UniProt 序列坐标，确认人源对应氨基酸
4. 准备 FASTA 文件（6 条序列）：
   - `Hsap_cGAS_WT`
   - `Hsap_cGAS_4mut`（人→裸鼹鼠）
   - `Hgal_cGAS_WT`
   - `Hgal_cGAS_4mut_rev`（裸鼹鼠→人）
   - `TRIM41_WT`

**风险与应对**：
- 若裸鼹鼠 cGAS 序列存在多个注释版本，可能需要从 NCBI 或论文 Supplementary 获取作者使用的精确序列。

---

### 阶段 1：结构预测（Day 3-5，需您协助）

**目标**：获得 cGAS-TRIM41 复合物的可靠三维结构。

**推荐工具：AlphaFold3 Server**（目前蛋白-蛋白复合物预测最准确）
- 网址：`https://alphafoldserver.com/`
- 需要 Google 账号登录
- 免费账户每日限 20 个预测
- 输入：两条 FASTA 序列（cGAS + TRIM41）

**需要您执行的操作：**
1. 注册/登录 Google 账号（若尚无）
2. 在 AF3 Server 上提交以下 4 组复合物（我可以准备好序列和说明，您复制粘贴即可）：
   1. 人 cGAS_WT + TRIM41_WT
   2. 人 cGAS_4mut + TRIM41_WT
   3. 裸鼹鼠 cGAS_WT + TRIM41_WT
   4. 裸鼹鼠 cGAS_4mut_rev + TRIM41_WT
3. 下载每个任务的 `ranked_0.pdb` 和 `confidence.json`

**备选工具（若 AF3 额度不足）**：
- 本地 ColabFold（AF-Multimer v3）：我可以安装，但 36GB 内存跑双链蛋白可能较慢且存在 OOM 风险。
- Boltz-1（开源 AF3 替代）：精度接近 AF3，可本地运行。

**质量控制**：
- 接受标准：`ipTM > 0.7` 且 `pTM > 0.8`
- 边缘标准：`ipTM 0.6-0.7`，仍可尝试但需在论文中明确说明预测置信度限制
- 拒绝标准：`ipTM < 0.6`，复合物预测不可靠，需切换为**替代方案 B**（见下方）

---

### 阶段 2：预分析与体系决策（Day 6-8）

**我负责执行：**

#### 2.1 界面初筛
- 用 `Arpeggio` 或自定义脚本计算：
  - 界面氢键、盐桥、疏水接触
  - 界面面积（ buried surface area ）
  - 4 个突变位点是否落在 5 Å 接触范围内
- 可视化检查：PyMOL 生成界面 overview 图

#### 2.2 关键决策：是否截取结构域？

基于 AF3 结构，判断界面是否集中在特定结构域：

| 情景 | 决策 | 依据 |
|------|------|------|
| 界面集中在 cGAS CTD + TRIM41 C-terminal 结构域 | **截取结构域做 MD** | 大幅减少原子数，可在几周内完成高质量重复 |
| 界面涉及 cGAS N-terminal 或 TRIM41 RING，且远程残基参与 | **保留全长** | 可能涉及变构，截断会丢失信息 |
| 界面非常弥散（多区域接触） | **保留全长，缩短轨迹** | 截取无法代表真实互作 |

**我的建议**：
- 对于 E3 连接酶（TRIM41），底物识别通常发生在 **C-terminal PRY-SPRY 或 coiled-coil 结构域**，而 RING 域负责 E2 结合。
- cGAS 的 4 个突变在 C-terminal 区域。
- **大概率可以安全截取 cGAS 的 C-terminal domain（约 250-300 aa）和 TRIM41 的底物识别结构域（约 200-300 aa）**，保留 interface ± 15-20 aa 的柔性 linker。
- 这样体系可从 ~1150 aa 降至 **~500-600 aa**，加水后约 **8-12 万原子**，速度可提升至 **~100 ns/day**。

#### 2.3 体系构建（Amber 流程）
1. `pdb4amber` 处理 PDB（加 missing atoms, 二硫键等）
2. `tleap` 加载 `ff19SB` 力场 + `OPC` 水模型
3. 溶剂化（12 Å 缓冲）+ 加离子（150 mM NaCl, 中性化）
4. 输出 `prmtop` / `inpcrd`
5. 转换为 OpenMM 格式（`AmberPrmtopFile`）

---

### 阶段 3：分子动力学模拟（Day 9-30，核心计算）

#### 3.1 模拟参数（OpenMM）
| 参数 | 设置 |
|------|------|
| 力场 | Amber ff19SB + OPC |
| 积分器 | LangevinMiddleIntegrator, 300K, 1/ps friction |
| 时间步长 | 4 fs（使用 HMR: hydrogen mass repartitioning） |
| 约束 | HBonds |
| 静电 | PME, 1.0 nm cutoff |
| 范德华 | 1.0 nm cutoff, switching at 0.9 nm |
| 压力耦合 | MonteCarloBarostat, 1 bar |
| 轨迹输出 | 每 10 ps 一帧，NetCDF 格式（via MDTraj/MDAnalysis） |

#### 3.2 模拟流程（每个体系）
1. **能量最小化**：5000 steps steepest descent
2. **升温**：NVT, 0→300K, 100 ps
3. **密度平衡**：NPT, 1 ns
4. **平衡**：NPT, 5 ns（位置约束主链，力常数 1.0 kcal/mol/Å²）
5. **生产运行**：NPT, 200 ns × 3 个独立重复（不同随机种子/初始速度）

#### 3.3 计算规模与耗时预估

**情景 A：截取结构域（~10 万原子，推荐）**
- 单轨迹 200 ns ≈ 2 天
- 4 体系 × 3 重复 = 12 轨迹
- 若串行：24 天；若并行 2-3 个作业（受内存限制）：**~10-14 天**

**情景 B：全长复合物（~25 万原子）**
- 单轨迹 200 ns ≈ 4-5 天
- 12 轨迹串行 ≈ 48-60 天，**超出预算**
- 若强行执行，需缩短至 100 ns × 2 重复，统计力度不足

**我的强烈建议：选择情景 A**。如果审稿人质疑截断，我们可以解释：
> "Based on AlphaFold3 prediction, the cGAS-TRIM41 interface is localized to the C-terminal domains. Domain-focused simulations (residues X-Y of cGAS and A-B of TRIM41) allow sufficient sampling for statistical robustness within computational constraints, while preserving the complete interface and adjacent flexible regions."

---

### 阶段 4：结合能计算与分析（Day 28-35）

#### 4.1 MM-GBSA
- 工具：`MMPBSA.py` (AmberTools) 或 `gmx_MMPBSA`（若用 GROMACS，但我们用 OpenMM，因此使用 AmberTools）
- 方法：
  - 从每轨迹抽取 5000-10000 帧（每隔 20-40 ps）
  - 计算 `ΔG_bind = G_complex - G_receptor - G_ligand`
  - 使用 `igb=5` (GB-OBC II) 或 `igb=8` (GB-Neck2)
- **关键输出**：
  - 4 个体系的绝对结合能趋势
  - **ΔΔG_bind = ΔG(突变体) - ΔG(WT)**
  - 若 |ΔΔG| > 1-2 kcal/mol，认为有显著差异

#### 4.2 能量分解（Per-residue decomposition）
- 分解到每个残基对结合能的贡献
- **重点**：4 个突变位点各自的 ΔΔG 贡献
- 识别界面 "hot spots"

#### 4.3 结构动力学分析
- **RMSD**：主链 RMSD 时间序列，评估收敛性
- **RMSF**：残基柔性，特别关注突变位点及其周围 loop
- **氢键/盐桥占据率**：界面关键相互作用的稳定性
- **界面面积变化**：SASA 时间序列
- **主成分分析（PCA）**：比较 WT 与突变体的构象空间差异

#### 4.4 静电与形状分析
- 用 `APBS` 或 PyMOL `APBS plugin` 计算分子表面静电势
- 比较人 WT vs 4mut 在 TRIM41 结合面上的静电变化

---

### 阶段 5：对照与突变扫描（Day 33-38）

#### 5.1 计算突变扫描（In silico mutagenesis）
工具选择（按推荐顺序）：
1. **Rosetta ddG**（`cartesian_ddg`）：免费学术版，安装包大但精度好
2. **FoldX**（`BuildModel` / `AnalyzeComplex`）：需要免费学术 license（在线申请，通常 1-2 天批准）
3. **PyRosetta**（若上述安装困难）：Python 接口，可脚本化

**执行内容**：
- 4 个位点各自的人↔裸鼹鼠双向突变
- 组合突变（4mut）
- 对照：界面上的非突变残基 Ala-scan（3-5 个位点）

**目的**：
- 与 MD 的 MM-GBSA 结果交叉验证
- 判断 4 个突变是协同作用还是独立作用

#### 5.2 阳性/阴性对照
- **阳性对照**：若文献中有 TRIM41 与其他底物的复合物结构（如 TRIM41-MAVS 等，如有），用同样方法预测并比较界面特征。
- **阴性对照**：将 cGAS 的一个非界面表面残基突变，预测应不影响结合（ΔΔG ≈ 0）。

---

### 阶段 6：论文产出（Day 38-45）

#### 6.1 图表清单（目标：4-5 张主图）
| Figure | 内容 |
|--------|------|
| Fig 1 | AlphaFold3 预测结构 + 界面 overview（4 个突变位点高亮） |
| Fig 2 | MD 收敛性（RMSD）+ RMSF（突变位点用箭头标出） |
| Fig 3 | 结合能比较：MM-GBSA ΔG 箱线图 + ΔΔG 热图 |
| Fig 4 | 代表帧界面细节：氢键网络、盐桥、突变位点侧链构象对比 |
| Fig 5 | Rosetta/FoldX 突变扫描结果与 MD 的相关性 |

#### 6.2 文稿
- 撰写 **Methods** 完整段落（软件版本、力场、模拟参数、分析流程）
- 撰写 **Results** 初稿
- 准备 **Supplementary**：各体系模拟参数表、RMSD 完整图、残基分解数据表

---

## 四、风险预案

### 预案 A：AF3 预测置信度低（ipTM < 0.6）
- **原因**：TRIM41 无实验结构，且 cGAS-TRIM41 互作可能是 transient/weak。
- **应对**：
  1. 不放弃，转用 **HADDOCK 2.4 / 2.5** 或 **ClusPro** 做蛋白-蛋白对接，以 AF3 的单个蛋白结构为输入，限制活性位点为 cGAS C-terminal。
  2. 同时执行 **单体 MD**：模拟人 WT vs 4mut cGAS 单体，看突变是否引起表面构象或动态变化（变构效应）。
  3. 论文表述改为："In the absence of a high-confidence complex structure, we employed an integrative approach combining docking, domain-focused MD, and monomer dynamics to probe the mechanism..."

### 预案 B：MD 轨迹不收敛（RMSD 持续漂移）
- **判断**：200 ns 后 RMSD 仍在上升，未出现平台期。
- **应对**：
  1. 检查是否因为截断结构域导致末端不稳定 → 延长截取的柔性 linker
  2. 使用 **aMD（accelerated MD）** 或 **GaMD** 增强采样，在更短时间内探索界面构象
  3. 若仅轻微漂移，可缩短分析窗口（如只用后 100 ns）

### 预案 C：结合能趋势与实验相反
- **即**：MM-GBSA 预测裸鼹鼠 cGAS 结合更强，但实验显示泛素化更弱。
- **科学解释**（可写入论文）：
  - 物理结合 ≠ 泛素化效率：TRIM41 可能仍能结合，但突变导致 E2-Ub 传递的几何不利（催化效率降低）。
  - 结合姿态（binding pose）改变：突变未削弱结合，但改变了 cGAS 上泛素化位点（Lys）相对于 TRIM41 RING 的空间取向。
  - MD 可以检测这种 "same binding, different geometry" 的现象。

---

## 五、需要您现在确认的关键决策

### 决策 1：结构预测途径（必选）
| 选项 | 说明 | 推荐度 |
|------|------|--------|
| **A. AlphaFold3 Server（网页版）** | 最准确，免费，但需您手动提交序列和下载结果 | **★★★ 强烈推荐** |
| B. 本地 ColabFold | 我可以全自动执行，但 M3 Pro 36GB 对 AF-Multimer 较紧张，速度慢 | ★★☆ |
| C. Boltz-1 本地 | 开源 AF3 级精度，可本地运行，我可以装 | ★★★ 若您不愿网页操作 |

### 决策 2：MD 体系策略（必选）
| 选项 | 说明 | 推荐度 |
|------|------|--------|
| **A. 智能截取结构域** | 若 AF3 显示界面集中，截取互作结构域做 200ns × 3 重复 | **★★★ 强烈推荐** |
| B. 全长复合物 | 更完整，但只能做 100ns × 2-3 重复，统计和收敛性较弱 | ★★☆ |

### 决策 3：目标论文定位
| 选项 | 说明 | 影响 |
|------|------|------|
| **A. 作为实验论文的 In silico 验证** | 如配合您的实验数据发 Cell Reports / Nature Communications 等 | **工作量适中，最可行** |
| B. 独立计算论文 | 发 PLOS Comp Bio / J. Chem. Inf. Model. 等 | 需要更多 method 创新和对照，5 周内紧张 |

### 决策 4：是否同时申请 FoldX 学术 License
- FoldX 是商用软件，学术免费但需在线申请（通常 1-2 工作日）。
- 若同意，我可以在 Day 1 就帮您准备申请材料（需您的机构邮箱和 PI 信息）。
- 若不愿申请，改用 **Rosetta**（完全免费，安装包 ~2GB）。

---

## 六、安装清单（执行首日）

若您确认方案，我将在 Day 1 完成：
- [ ] 安装 Miniforge (arm64)
- [ ] 安装 OpenMM + AmberTools + 分析库
- [ ] 安装 PyMOL (open-source) 或 ChimeraX
- [ ] 确认 M3 Pro Metal 后端正常工作（跑一个测试体系）
- [ ] 准备所有 FASTA 文件
- [ ] 精确映射 4 个突变位点

---

## 七、预期科学产出

1. **原子级解释**：4 个突变位点是否直接接触 TRIM41，或以变构方式影响界面。
2. **定量趋势**：WT vs 4mut 的 ΔΔG_bind 估算（半定量，~1-2 kcal/mol 精度）。
3. **机制假设**：
   - 若结合能显著降低 → "物理结合减弱导致泛素化减少"
   - 若结合能不变但界面几何改变 → "结合姿态改变导致泛素化催化效率下降"
   - 若柔性增加 → "动态构象变化影响 E2-Ub 传递"
4. **可发表的数据集**：12 条 200ns 轨迹 + 结合能统计 + 突变扫描对照。

---

---

## 附录：执行状态更新（2026-05-02）

| 原计划阶段 | 实际状态 | 备注 |
|-----------|---------|------|
| 阶段 0–2（环境 + 结构预测 + 预分析） | ✅ 完成 | AF3, LightDock, Rosetta docking 全部完成 |
| 阶段 3（MD 模拟） | 🔄 进行中 | Hsap_WT/4mut (6 reps) ✅；Hgal old (5 reps) ✅；Hgal NEW (0/6 reps) 刚启动 |
| 阶段 4（结合能 + 分析） | 🔄 部分完成 | MM-GBSA (Hsap 6 reps) ✅；深度分析 (WT vs S305-phos) ✅；Hgal 批量分析待 NEW 数据 |
| 阶段 5（对照与突变扫描） | ✅ 部分完成 | Rosetta FastRelax 突变扫描完成 |
| 阶段 6（论文产出） | ⏳ 待启动 | 需 Hgal NEW 6 reps 完成后开始 |

**关键调整**：
1. 硬件从 Apple M3 Pro 迁移至 Linux + 4× RTX 3090
2. 新增磷酸化研究（S305-phos/S305E）作为 side validation，不纳入主叙事
3. Hgal NEW 构建（Apr 28, Rosetta pose + solvateOct）取代 OLD 构建（LightDock + solvateBox）
4. GROMACS 交叉验证完成 1 rep，与 OpenMM 0–150ns 高度一致

*方案撰写时间：2026-04-22*
*最后更新：2026-05-02*
*基于资源评估和文献现状制定*
