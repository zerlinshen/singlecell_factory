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

## 8) 交付输出建议
每次交付固定 4 项：
1. 改了哪些文件
2. 跑了哪些验证命令
3. 产物证据（尤其 manifest/status）
4. 剩余风险与下一步
