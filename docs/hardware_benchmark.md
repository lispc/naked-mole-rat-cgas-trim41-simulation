# cGAS-TRIM41 MD 研究项目日志

> 本文档记录项目执行过程中的关键决策、数据、估算和推理，供后续论文撰写和复盘使用。
> 
> 最后更新：2026-04-22

---


## 四、硬件资源评估与性能基准

### 4.1 硬件规格（原 Apple M3 Pro）

| 项目 | 规格 |
|------|------|
| CPU | Apple M3 Pro (12 核) |
| GPU | 18 核 Apple Silicon GPU (通过 OpenCL 访问) |
| 内存 | 36 GB 统一内存（无独立显存） |
| 操作系统 | macOS 15.3 (Darwin 25.3.0) |
| Python | 3.13.13 (conda-forge) |

### 4.1b 硬件规格（迁移后 Linux + RTX 3090）

| 项目 | 规格 |
|------|------|
| CPU | x86_64 (未详测) |
| GPU | 4× NVIDIA GeForce RTX 3090 (24GB VRAM each) |
| 内存 | 系统内存充裕 |
| 操作系统 | Linux (CUDA 13.0) |
| Python | 3.11.15 (conda-forge) |

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

### 4.2b RTX 3090 实测基准（迁移后）

#### 实测 1：alanine dipeptide vacuum（验证 CUDA 后端）
- 10k steps @ 2fs = 20 ps
- **CUDA**: 20 ps in 0.54s → **18,612 steps/s**
- 仅为理想上限参考

#### 实测 2：Hgal_domain 生产 MD（116,710 atoms, PME, HBonds, 2fs）
- 采样窗口: 120s 内跑了 0.2ns
- **实测速度: ~152 ns/day**
- 包含 PME、溶剂、离子、DCD I/O、checkpoint 等全部真实开销

**与 M3 Pro 对比**:

| 指标 | M3 Pro (OpenCL) | RTX 3090 (CUDA) | 加速比 |
|------|----------------|-----------------|--------|
| NVT Heating | ~46 ns/day | ~150-200 ns/day | **3-4x** |
| NPT Equil | ~21 ns/day | ~100-150 ns/day | **5-7x** |
| NVT Production | ~40-50 ns/day | **~152 ns/day** | **3-4x** |
| 单条 200ns 耗时 | ~4-5 天 | **~1.3 天** | **3-4x** |
| 3 重复总计 | ~12-15 天 | **~1.3 天** | **~10x** (并行优势) |

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

