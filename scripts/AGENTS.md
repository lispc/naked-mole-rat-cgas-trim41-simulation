# Scripts Agent Instructions

本目录包含项目全部 Python / Shell / PyMOL 脚本。

## 规范

### Python
- 命名：`snake_case`，函数名要有描述性
- 字符串：优先用 f-string
- 路径：不要硬编码绝对路径（如 `/home/user/...`），用 `pathlib.Path` 或相对路径
- 错误处理：IO 操作和外部命令调用需有 try/except
- Docstring：新脚本开头必须说明用途和用法示例

### 脚本状态标记
在 `scripts/README.md` 中维护状态：
- **Active** — 当前可用，近期维护
- **推荐** — 同功能多版本中的首选
- **Broken** — 已知问题，需修复
- **Deprecated** — 已废弃，不再维护

### 环境
- 所有 Python 脚本在 `cgas-md` conda 环境中运行
- GROMACS 相关脚本在 `gmx` conda 环境中运行

### 边界
- ✅ **Always**: 修改脚本后更新 `scripts/README.md` 中的描述
- ⚠️ **Ask first**: 删除脚本、重命名已有脚本
- 🚫 **Never**: 在脚本中内嵌密码或 API key

## 完整索引
见 [`README.md`](README.md)。
