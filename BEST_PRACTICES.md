# singlecell_factory Best Practices (Codex + Claude)

## 1) 先看作用域
- Codex 项目规则：`AGENTS.md`（本项目根目录）
- Claude 项目规则：`CLAUDE.md`（本项目根目录）
- 全局规则仍生效，但项目文件优先级更高。

## 2) 推荐任务起手式
把下面模板发给代理（Codex/Claude 都可）：

```text
$ralplan
目标：<要完成的分析或改造>
约束：不破坏 mandatory pipeline chain，不新增依赖（除非明确批准）
验收标准：给出测试结果 + 运行产物证据
执行要求：实现后必须 verifier 复核
```

## 3) Skills 最佳路由
- 跑流水线/失败恢复：`$execute-and-recover-pipeline`
- 新增/改造模块：`$develop-and-integrate-module`
- 性能优化：`$optimize-and-guard-performance`
- 大任务编排：`$orchestrate-large-task`
- 最终质量门：`$review-and-validate-quality`

建议固定链路：
- 开发类：`$ralplan -> $develop-and-integrate-module -> $review-and-validate-quality`
- 运行类：`$execute-and-recover-pipeline -> $review-and-validate-quality`

## 4) Agents 使用建议
- `explore`：快速定位模块/符号
- `executor`：实现改动
- `test-engineer`：补齐测试策略
- `verifier`：最终证据闭环
- `security-reviewer`：仅在高风险输入/外部交互变更时启用

## 5) 验证标准（最少要做）
- 快速回归：`pytest -q tests/test_modular.py tests/test_modular_optimizations.py`
- 全量回归（高风险变更）：`pytest -q`
- 流水线改动必须附带：
  - `module_status.csv`
  - `run_manifest.json`

## 6) Hook 使用建议
- 默认全局 hooks 已自动生效（无需手动触发）。
- 如果你用 tmux 持续运行模式（ralph/team/ultrawork），先配置 `.omx/tmux-hook.json` 的 pane id，否则会是 `invalid_config`。

## 7) 绝对不要做的事
- 不要破坏 mandatory chain：`cellranger -> qc -> doublet_detection`
- 不要跳过文档同步（`README.md` / `PROTOCOL.md`）
- 不要在无验证证据下宣称“完成”

## 8) Densify Policy（Phase 7A.2+，强制）

所有在 `workflow/modular/modules/` 下的 `.toarray()` / `.todense()` 调用必须遵守以下协议：

### 规则

1. **新代码必须走 `safe_densify()`**（来自 `_sparse_utils`）或主动调用 `plan_densify()`（来自 `_densify_policy`）判断 budget，拒绝后 raise `MemoryGuardError` 或改为 chunked 路径。

2. **已有代码（存量 allowlist）**：每个 `.toarray()` / `.todense()` 必须在同一行或紧邻上一行添加注释：
   ```python
   # densify-allowed: <reason explaining why this densification is bounded/safe>
   ```

3. **CI grep-ban**：`tests/test_densify_audit.py` 在每次 commit 时扫描 `modules/`，任何缺少标记的调用导致测试失败。

### 判断标准（写 reason 时参考）

| 场景 | 可接受理由 |
|---|---|
| 单列切片（n_cells × 1） | "single-gene column vector; trivially small" |
| 有界子集（top-k 基因） | "subset of top_k genes × n_cells; bounded by n_top_genes" |
| 引擎 flag 保护 | "only reached when SC_XX_ENGINE != 'chunked'; caller controls RAM" |
| obsm 低维矩阵 | "already reduced dimensionality; n_cells × n_cnv_bins" |
| 整矩阵但有上游 flag | 必须说明上游 flag + 说明无 chunked 替代 |

### 添加新模块时

```python
from workflow.modular._densify_policy import DensifyDecision, plan_densify
from workflow.modular._sparse_utils import MemoryGuardError

decision = plan_densify(X.shape, np.float32, reason="my reason")
if decision == DensifyDecision.ABORT:
    raise MemoryGuardError("...")
elif decision == DensifyDecision.CHUNK:
    # use chunked_row_densify or implement chunked path
    ...
else:
    dense = X.toarray()  # densify-allowed: covered by plan_densify GO decision
```

### 环境变量（覆盖默认 cap）

| 变量 | 默认 | 说明 |
|---|---|---|
| `SC_DENSIFY_SOFT_CAP_BYTES` | 4 GiB | 超过此值返回 CHUNK |
| `SC_DENSIFY_HARD_CAP_BYTES` | 12 GiB | 超过此值返回 ABORT |

## 9) 交付输出建议
每次交付固定 4 项：
1. 改了哪些文件
2. 跑了哪些验证命令
3. 产物证据（尤其 manifest/status）
4. 剩余风险与下一步

## 10) R Contract Test Discipline (Phase 7D)

`tests/test_r_bundle_contract.py` contains subprocess-based round-trip tests that export a synthetic v2 bundle from Python and read it back via `io_bundle.R` using `Rscript -e`. These tests guard against silent breakage of the Python exporter / R reader contract.

- All tests are tagged `@pytest.mark.r_contract` and excluded from the default CI run via `pyproject.toml` `addopts`.
- Run explicitly with: `pytest -m r_contract tests/test_r_bundle_contract.py`
- Tests auto-skip when Rscript is absent or required R packages (`arrow`, `jsonlite`, `Matrix`) are not installed — they never fail hard due to missing R environment.
- R output is parsed by scanning stdout for `KEY=VALUE` lines emitted via `sprintf` + `cat`. Never mix informational R output with these sentinel lines.
- The `rscript_path` session fixture in `conftest.py` handles both binary detection and R package probing.
- R source under test lives in `multiomics_r_factory/R_bundle/io_bundle.R` — never edit it from this repo (bridge symlink rule applies).
