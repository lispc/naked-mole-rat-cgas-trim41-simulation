# Good Ideas Backlog

> 高优先级任务完成后，回过头来做的实验/分析
> 创建：2026-04-24

---

## 计算类

### 1. Rosetta ΔΔG 精确计算
- **方法**: `cartesian_ddg` 或 `ddg_monomer`
- **目的**: 精确计算每个突变对稳定性的影响
- **优先级**: P1
- **依赖**: Rosetta 环境已就绪

### 2. Rosetta Per-Residue Energy Decomposition
- **方法**: 对 WT 和 4mut 复合物分别做能量分解
- **目的**: 看 4mut 位点周围的能量扰动如何传播到 N 端界面
- **优先级**: P1
- **依赖**: 现有 Rosetta 对接结构

### 3. Markov State Model (MSM)
- **方法**: pyEMMA 或 OpenPathSampling
- **目的**: 如果轨迹延长到 500ns+，建立 MSM 看亚稳态转换
- **优先级**: P2
- **依赖**: 轨迹长度 > 500ns，或增强采样数据

### 4. AlloSigMA / gRINN 变构网络分析
- **方法**: AlloSigMA（基于结构）或 gRINN（基于 MD 轨迹）
- **目的**: 专门的变构效应检测工具，计算残基相互作用网络
- **优先级**: P2
- **依赖**: 需要安装额外包

### 5. 增强采样（Metadynamics / Umbrella Sampling）
- **方法**: OpenMM + PLUMED，或 PySAGES
- **目的**: 沿变构反应坐标采样，计算自由能面
- **优先级**: P2
- **依赖**: 需要定义合适的集体变量（CV）

### 6. MM-GBSA / MM-PBSA 结合能计算
- **方法**: MMPBSA.py (AmberTools)
- **目的**: 计算 WT vs 4mut 的结合自由能差异
- **优先级**: P1
- **依赖**: Hsap MD 轨迹

---

## 湿实验建议

### 7. K414R 定点突变验证
- **方法**: 将 cGAS Lys-414 突变为 Arg
- **目的**: 验证 Lys-414 是否为 TRIM41 泛素化的真实位点
- **优先级**: P1（如果合作者有资源）
- **背景**: Lys-414 距离 N 端界面仅 ~6.7Å，是最可能的泛素化位点

### 8. 嵌合体实验（修正后的 4mut 位点）
- **方法**: 基于正确的 D431S 位点重新设计嵌合体
- **目的**: 测试长 loop（462-494）对紧凑几何的贡献
- **优先级**: P2
- **背景**: 之前的嵌合体设计基于错误的 C463S 位点

---

## 文献调研

### 9. TRIM41 coiled-coil domain 是否参与 cGAS 识别
- **来源**: 引用 28（Zhen et al.）提到 coiled-coil 对 TRIM41-ORF2p 必需
- **问题**: 我们的 MD 只包含 SPRY 域（413-630），是否遗漏了关键相互作用？
- **优先级**: P1
- **影响**: 如果 coiled-coil 参与，需要重新考虑结构域截断策略

### 10. cGAS 已知泛素化位点完整列表
- **来源**: UniProt Q8N884 + PhosphoSitePlus
- **目的**: 确认是否有其他 Lys 距离 N 端界面 < 10Å
- **优先级**: P2

---

## 可视化

### 11. 变构效应路径动画
- **方法**: PyMOL 或 VMD 制作 WT → 4mut 的形变动画
- **目的**: 论文 supplementary video
- **优先级**: P3

### 12. 动态交叉相关矩阵（DCCM）热图
- **方法**: MDAnalysis + matplotlib
- **目的**: 展示 C 端 ↔ N 端动态耦合
- **优先级**: P0（Hgal MD 完成后立即做）

---

*状态：待完成 | 维护者：Kimi Code CLI*
