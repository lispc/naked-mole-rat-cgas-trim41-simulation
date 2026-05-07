# 项目回顾：得失总结

> 项目周期：2026-04-22 — 2026-05-07（17 天）  
> 基于 project_log.md（§1-61）和 project_log_2026_04.md 的完整回顾

---

## 一、时间线概览

| 阶段 | 日期 | 关键事件 |
|------|------|---------|
| 启动 | 04-22 | 阅读 Chen et al. Science 2025，制定计算方案 |
| 结构预测+对接 | 04-23 | AF3 预测、LightDock/Rosetta/ClusPro 对接、发现界面在 N 端 |
| MD 构建与首批模拟 | 04-24~30 | tleap 构建、OpenMM MD 启动、GROMACS 设置 |
| 调试与修复 | 05-01~02 | CMAP bug 发现与修复、Python 对齐 bug、MM-GBSA、S305-phos |
| 四系统 MD | 05-02~05 | Hgal/Hsap WT/4mut 四系统 200ns×3 完成 |
| 四元复合物迭代 | 05-03~06 | v1 失败 → v2 (+RING+异肽键) → FULL (+SPRY+CC) |
| US/PMF + 论文 | 05-06~07 | Steered MD → umbrella sampling → PMF + 催化准备度分析 |

---

## 二、做得好的

### 1. 文档文化
3600+ 行的项目日志，按 § 编号记录每一次实验、每一个 bug 修复、每一个待决策事项。失败记录和成功记录同样详细——ClusPro 失败、SDOCK2.0 超时、auto-launcher 的 race condition、4mut 的 NaN 崩溃——全部可追溯。这种诚实度使得后来者（包括几天后的自己）可以精确理解每个决策的背景。

### 2. 多方法交叉验证
从不信任单一模型或单一方法。AF3 + Boltz-2 + Chai-1 三者交叉验证复合物界面（结果：界面预测互相矛盾，验证了对接+MD 策略的必要性）。LightDock + Rosetta + ClusPro 三种对接方法。OpenMM + GROMACS 两个引擎。MM-GBSA + Rosetta I_sc + MD 几何分析三个层面互相印证。这一条是项目结论可信度的基础。

### 3. 自我纠错文化
- **CMAP bug**（05-01）：发现 parmed 将 14 种残基特异性 CMAP 压缩为 1 种，GROMACS RMSD 被高估 4 倍。修复后 ratio 降至 1.2×。
- **Python 对齐 bug**（05-02）：rotation_matrix 在未中心化坐标上计算，RMSD 被高估 82%。修正后与 `gmx rms` gold standard 一致。
- **突变编号错误**：C463S → D431S。发现后全局修正了 fasta、脚本、文档。
- **t-test 误用**：对自相关 MD 轨迹直接用独立样本 t-test → 改用 correlated t-test。
- **COM restraint 水原子 bug**（05-06）：CENTROID 计算包含了 12 万水原子，导致速度从 122 降到 29 ns/day。修复后恢复正常。

每次发现错误 → 承认 → 修复 → 更新文档的循环执行得很彻底。

### 4. 渐进式模型构建
不是一步到位，而是逐层加复杂度：
1. 二元复合物（cGAS + TRIM41 SPRY）——验证界面稳定性和 4mut 效应
2. v1 最小四元（E2~Ub + cGAS）——发现 E2~Ub 散开，确认识别到需要 RING 和异肽键
3. v2（+RING + 异肽 bond restraint）——E2~Ub closed，K315 从 29→11 Å
4. FULL（+SPRY + CC linker）——K315 稳定在 ~20 Å，零漂移
5. US/PMF（沿 K315→Ub 距离）——定量自由能面 + 催化准备度

每一层都回答了前一层提出的问题，同时暴露了下一层需要的改进。

### 5. 外部评审
主动请了 Kimi 和 Gemini 做 code/doc review（04-27~28），发现了 t-test 误用、突变编号未清理、docking pose bias、盐浓度注释不符等致命/严重问题。逐条验证并修复。这种"自己花钱请人找茬"的做法在学术计算项目中很少见。

### 6. 知道什么时候停
S305 磷酸化 story 跑了 1119ns MD + MM-GBSA，发现了明确的解离信号，但意识到这是正交问题（不是 4mut 机制的核心），果断停止扩展（§51.3："停止所有新的磷酸化相关实验"）。Auto-launcher 系统写了但发现 race condition 和 CUDA context 问题后也及时停用（§53.5）。

---

## 三、做得不好的

### 1. 早期突变编号混乱
项目最初用 C463S/K479E/L495Y/K498T（基于错误的裸鼹鼠→人映射）。直到 external review 指出才纠正为 D431S。这导致：
- 早期 fasta 文件、脚本、部分文档使用了错误编号
- `af3_mutation_analysis.md` 中的 Hgal_rev RMSD 数据（13.04 Å）是基于错误位点的——后来证明纠正后 RMSD 只有 0.98 Å，完全不同的结论
- 清理不彻底——旧错误 fasta 长期保留在项目中

**教训**：项目第一天就应该做全局序列比对确认编号，而不是在 review 阶段才发现。

### 2. 首次四元复合物尝试的工程失误
第一次 build_quaternary_mvp.py 尝试（05-03）的问题：
- 用 Kabsch 刚体对齐做 4mut → 灾难性 clashes（cGAS 穿透 RING, 999+ contacts）
- tleap/pdb4amber 后坐标漂移 10 Å（K315 从 16.6→26.3 Å），至今不清楚原因
- 没有先做 frame-0 验证（"K315 距离检查会立刻发现坐标漂移"——项目自己的 lessons learned）

**教训**：复杂结构拼装要先做 clash check + 坐标验证，再投入 GPU 时间。永远检查 frame 0。

### 3. S305 磷酸化 story 消耗了不成比例的资源
- S305-phos 600ns + S305E 600ns = 1200ns，加上 MM-GBSA 和深度分析
- 结论在 200ns 时就已明确（完全解离）
- 对论文主线的贡献是 side note（"5.2.7 磷酸化作为正交开关"），不应该是 1200ns 的投入

**教训**：明确主线后，side projects 应该用最小可行实验回答，达到结论就停。

### 4. 论文草稿过早 + 质量差
04-26 左右就写了 paper.md 和 main.tex，但到 05-06 时：
- 7 条参考文献中 5 条是 [PLACEHOLDER]
- Table 7 一整行 [PLACEHOLDER]
- Hgal 数据写着"in progress"
- 突变编号还是错的（C463S）
- "Tight-but-Floppy"概念从未被明确定义却出现在 Discussion 中

**教训**：论文应该在核心数据齐备后写，否则早期草稿的框架会限制思维（"我们已经有这个 section 了，不舍得删"），或者变成不断修补僵尸框架。

### 5. lib/ 僵尸代码
§50.5 花了一小时整理了 scripts/lib/（paths.py, stats.py, mda_tools.py, plot_style.py），但直到 05-06 的组织整理中才发现——**没有一个生产脚本实际 import 它**。33 个脚本仍然各自硬编码 `/home/scroll/...` 路径。

**教训**：创建共享库的同时必须立即把它接入至少 2-3 个实际脚本，否则它永远是"下周我会用"的僵尸。

### 6. PDB 格式对齐的反复踩坑
Bio.PDB 的列对齐问题在项目中出现了至少 3 次：
- mutate.py 写出的 PDB 残基名变成 YS 而非 LYS（列错位）
- pdb4amber 输出的 chain ID 丢失
- tleap 后的 residue ID 与原始 PDB 不同

每次都用 30-60 分钟 debug。如果一开始就写一个经过列验证的 PDB writer，这个问题只需要处理一次。

---

## 四、关键转折点

| 转折 | 日期 | 影响 |
|------|------|------|
| 界面不在 C 端 | 04-23 | 推翻原始假设，确立变构方向 |
| AF3 突变序列验证 | 04-23 | 推翻"突变驱动几何改变"假说 |
| CMAP bug 发现 | 05-01 | GROMACS 数据从"不可用"变"可用" |
| 四系统 MD 对比 | 05-05 | 跨物种证据：Hgal 敏感、Hsap 稳健 |
| FULL 四元构建 | 05-06 | SPRY 锚定→K315 稳定、零漂移 |
| US PMF 漂移方向 | 05-07 | 定量自由能证据、催化准备度 4.3× |

---

## 五、如果要重来一次会怎么做

1. **D1 做全局序列比对**，确认突变编号，避免后续清理
2. **先跑 10ns 测试所有体系**，确认 stable 再扩展到 200ns×3
3. **四元复合物直接用 FULL 模型**（+RING+SPRY+CC），跳过 v1 和 v2 的弯路
4. **S305 停在 200ns**，不扩展到 600ns+600ns
5. **论文等数据齐了再写**，先写 Methods 和 Introduction，Results 最后
6. **lib/ 创建后立即接入** 3 个最常用脚本，设 CI 检查
7. **写一个经过验证的 PDB writer** 避免 Bio.PDB 列对齐反复踩坑

---

*最后更新：2026-05-07*
*维护者：Claude Code CLI*
