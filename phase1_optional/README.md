# Deprecated Alias: phase1_optional

该目录为兼容别名入口，已迁移到 `RNA Velocity and Pseudotime Analysis`。

## 新模块位置
- [README.md](file:///home/zerlinshen/singlecell_factory/rna_velocity_pseudotime_analysis/README.md)
- [USAGE.md](file:///home/zerlinshen/singlecell_factory/rna_velocity_pseudotime_analysis/docs/USAGE.md)
- [MIGRATION.md](file:///home/zerlinshen/singlecell_factory/rna_velocity_pseudotime_analysis/docs/MIGRATION.md)

## 兼容入口
- `phase1_optional/scripts/run_phase1_optional.sh` 会转发到新运行脚本
- `phase1_optional/scripts/cleanup_phase1_artifacts.sh` 会转发到新清理脚本

## 推荐使用
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3
```
