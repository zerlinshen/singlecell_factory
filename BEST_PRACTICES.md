# singlecell_factory Best Practices (Codex + Claude)

## 1) 先看作用域
- Codex 项目规则：`AGENTS.md`（本项目根目录）
- Claude 项目规则：`CLAUDE.md`（本项目根目录）
- 全局规则仍生效，但项目文件优先级更高。

## 1.5) 两种操作方式都可行
- 直接在远端 `/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory` 操作：适合重计算、
  pipeline execution、run triage、remote R/reporting。
- 在 Mac 上通过 SSH 编排远端：适合协调、审阅、handoff、report packaging。

无论哪种方式，远端 artifact 仍是 run truth：优先看
`run_manifest.json`、`module_status.csv`、`ops/before_every_run/LATEST.md`、
`ops/run_ledger/`。Mac 侧 summary 负责组织和交付，不替代远端证据。

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
- 远端执行前记忆：`$before-every-run`
- Mac/remote 桥接与 handoff：`$singlecell-remote-workflow`
- 报告包整理：`$remote-run-report-bridge`
- 文档/治理漂移：`$reproduction-workspace-governance` + `$readme-sync-enforcer`
- 失败/过期结果治理：`$post-run-failure-cleanup`、`$reproduce-run-retention`

建议固定链路：
- 开发类：`$ralplan -> $develop-and-integrate-module -> $review-and-validate-quality`
- 运行类：`$execute-and-recover-pipeline -> $review-and-validate-quality`
- 远端运行类：`$before-every-run -> $singlecell-remote-workflow -> $execute-and-recover-pipeline -> $review-and-validate-quality`
- handoff/report 类：`$remote-run-report-bridge -> $reproduction-workspace-governance -> $review-and-validate-quality`

Skill 目录：
- Codex global: `~/.codex/skills/`
- Codex project: `.codex/skills/`
- Claude global: `~/.claude/skills/`
- Claude project: `.claude/skills/`

`codex_skills/` 只作为历史 project-local skills 的兼容镜像；新的或刷新后的
Codex project skills 放在 `.codex/skills/`。

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
- 不要把 `2026-04-25` 的 sparse-exact probe 目录当成 canonical success，
  除非看到 `final_adata.h5ad`、`run_manifest.json`、`module_status.csv`

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
- R source under test lives in `r_multiomics_factory/R_bundle/io_bundle.R` — never edit it from this repo (bridge symlink rule applies). R plot helpers live in `plotting_factory/r/` (introduced 2026-05-18; bridged via `bridges/local_plot_pipeline/`); never edit them inside `r_multiomics_factory/R/` or `singlecell_factory/workflow/`.

## 11) Staircase Testing Discipline (Phase 7C)

为防止 5k synthetic 测试无法触发的真实数据 bug（典型例子：doublet auto-threshold 在 100k 才生效），所有改动按"最低 tier"原则跑测试：

| 改动类别 | 必经 gate |
|---|---|
| 纯数值算法（DE / parity / 单 module 改进） | nano (默认 `pytest`) |
| 内存 / 密度 / chunked 路径 | nano + `pytest -m small_real` |
| 论文方法对齐 / capability flag 改动 | nano + small_real + `pytest -m medium_real` |
| 发表前最终复现 | nano + small_real + medium_real + `pytest -m full_real` + paper-aligned 完整 launcher |

工具链:
- `scripts/build_staircase_fixtures.py` 从 NC2024 prepared zarr 生成 `tests/data/staircase/{small_real,medium_real}.zarr` 子集（按 sample 分层抽样，seed=42 确定性）。运行一次即可，结果不入仓 (`tests/data/.gitignore` 排除大文件)。
- `tests/test_staircase_smoke.py` 是入口测试集合，验证 grouped Scrublet 阈值触发、disease 列覆盖等真实数据特性。
- 每次添加新 module 或修改内存敏感路径时必须扩展 staircase smoke 测试。
- pytest markers 注册在 `pyproject.toml`；默认 `addopts` 排除 `small_real / medium_real / full_real`，CI 速度不受影响。
- nano 是 synthetic（5k cells × 2k genes，density 0.1，CSR），永远在默认 suite 中跑——保证最快回归信号。
- full_real 不复制数据，写 `tests/data/staircase/full_real.path` 路径指针指向原 zarr。

参见 `tests/data/staircase/README.md` 获取每个 tier 的详细规模、用途与生成命令。
