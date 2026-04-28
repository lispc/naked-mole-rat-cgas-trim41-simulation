# Project Audit Report: cGAS-TRIM41 Allostery Study
**Date:** 2026-04-28
**Reviewer:** Gemini CLI (Review Phase 1)

## 1. 核心发现：致命逻辑漏洞 (Critical Vulnerabilities)

### 1.1 突变序列执行错误 (Sequence Implementation Error)
- **发现**：在 `sequences/Hsap_cGAS_4mut.fasta` 中，声称执行了 `D431S` 突变，但通过序列比对发现，实际发生突变的位置是 **D436**（或者是由于编号偏移导致的错误位点），而真实的 **D431** 位点在序列中依然是 Aspartate。
- **后果**：如果 Hsap 4mut 的 MD 和 US 是基于这个序列跑的，那么该体系模拟的是一个无关的突变，无法验证 Science 2025 的结论。

### 1.2 跨体系对比逻辑失效 (Comparison Logic Breakdown)
- **发现**：`scripts/compare_systems.py` 直接按残基编号（Residue ID）进行减法运算。
- **问题**：裸鼹鼠（Hgal, 554aa）与人源（Hsap, 522aa）存在 32 个残基的系统性偏移。直接相减会导致功能上不对应的残基被放在一起对比。
- **后果**：所有 ΔRMSF 曲线和 ΔContact 结论均存在物理位置错位，无法支持变构效应路径的推导。

### 1.3 活性位点编号“幻觉” (Residue Numbering Ghosting)
- **发现**：`README.md` 和分析脚本参数中频繁出现 `+19` 的偏移（如将 S463 称为 482）。但在原始 PDB 中，463 号确实是 Serine。
- **风险**：如果分析时使用了 482 这个编号，测量的是 **Proline 482** 到 TRIM41 的距离，而非核心突变位点。

---

## 2. 技术细节审计 (Technical Observations)

### 2.1 采样与收敛 (Sampling & Convergence)
- **Rosetta Docking**: 初始仅使用 `nstruct=10`，严重采样不足。后期虽补做了 `nstruct=100`，但对于蛋白-蛋白对接，仍处于统计低限。
- **Umbrella Sampling**: 力常数 `k=1000` 可能在某些剧烈构象变化区域不足以维持窗口重叠。

### 2.2 构建脚本问题 (Build Script Inconsistencies)
- **盐浓度**: `build_system.py` 注释声称添加 150mM NaCl，但代码仅执行了电荷中和（`addIons Na+ 0`）。
- **DCD 频率**: 发现过 DCD 写入频率过高的 Bug（每 1ps 存一帧），导致磁盘空间异常浪费。

---

## 3. 下一步行动建议 (Action Plan for 3090 Session)

### 阶段 A：数据有效性核查 (Data Validation)
1.  **核查 3090 上的轨迹编号**：确认 `Hsap_WT_summary.json` 等文件中，`active_sites` 对应的 resid 到底指向哪个残基。
2.  **核查 US 重叠度**：运行 `scripts/run_wham.py` 之前，务必先生成直方图，检查 CV 分布是否连续。

### 阶段 B：工具链修正 (Tooling Fix)
1.  **重写对齐逻辑**：建立 `mapping.json`，基于比对（Alignment）后的索引进行 `compare_systems.py` 的差值计算。
2.  **修正突变序列**：如果确认 `4mut.fasta` 序列错误，需要重新生成 Fasta -> AF3 -> 系统构建 -> MD。

### 阶段 C：结果导出 (Export Strategy)
- 考虑到数据量巨大，建议在 3090 上运行完分析后，仅将 `data/analysis/*.json` 和 `*.png` 压缩同步回本地进行最终论文图表排版。

---
**Status:** Review phase complete. Ready to resume on RTX 3090 machine.
