# Migration Guide

## 命名变更
- 旧模块名：`Phase1 Optional Analysis Workspace`
- 新模块名：`RNA Velocity and Pseudotime Analysis`

## 路径映射
- `phase1_optional/` -> `rna_velocity_pseudotime_analysis/`
- `phase1_optional/scripts/run_phase1_optional.sh` -> `rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh`
- `phase1_optional/scripts/cleanup_phase1_artifacts.sh` -> `rna_velocity_pseudotime_analysis/scripts/cleanup_rna_velocity_pseudotime_artifacts.sh`
- `phase1_optional/configs/phase1.optional.yaml` -> `rna_velocity_pseudotime_analysis/configs/rna_velocity_pseudotime.yaml`
- `phase1_optional/envs/phase1.environment.yml` -> `rna_velocity_pseudotime_analysis/envs/rna_velocity_pseudotime.environment.yml`

## 向后兼容
- 旧目录 `phase1_optional/` 保留兼容入口。
- 旧运行脚本与清理脚本会转发到新脚本。
- 旧配置路径保留为兼容别名。

## 推荐迁移动作
1. 将文档和自动化命令中的旧路径替换为新路径。
2. 在 CI/CD 中排除 `rna_velocity_pseudotime_analysis/**`。
3. 保留旧路径一段过渡期后再移除。
