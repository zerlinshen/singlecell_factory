# RNA Velocity and Pseudotime Analysis Usage Guide

## 定位声明
本目录中的能力为可选分析步骤，非生产必需，不属于主流程交付范围。

## 运行前准备
1. 创建独立环境：
```bash
conda env create -f rna_velocity_pseudotime_analysis/envs/rna_velocity_pseudotime.environment.yml
```
2. 激活独立环境：
```bash
conda activate sc_rna_velocity_pseudotime
```

## 执行方式
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3
```
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3 --prepare-env
```
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3 --env-name sc_rna_velocity_pseudotime
```

## 从BAM生成loom
```bash
bash rna_velocity_pseudotime_analysis/scripts/generate_velocity_loom.sh \
  --samplefolder /home/zerlinshen/singlecell_factory/pbmc_1k_v3_count \
  --gtf /home/zerlinshen/singlecell_factory/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
  --threads 4
```
生成后将loom路径写入配置的 `input.loom`。

## 输出位置
- 日志：`rna_velocity_pseudotime_analysis/runtime/logs/`
- 报告：`rna_velocity_pseudotime_analysis/runtime/reports/<run_id>/`
- 结果：`rna_velocity_pseudotime_analysis/runtime/results/<run_id>/`
- 缓存：`rna_velocity_pseudotime_analysis/runtime/cache/<run_id>/`
- 图像：`rna_velocity_pseudotime_analysis/runtime/results/<run_id>/figures/`
- 结果表：`rna_velocity_pseudotime_analysis/runtime/results/<run_id>/tables/`

## 清理策略
```bash
bash rna_velocity_pseudotime_analysis/scripts/cleanup_rna_velocity_pseudotime_artifacts.sh --keep-days=7
```
或
```bash
bash rna_velocity_pseudotime_analysis/scripts/cleanup_rna_velocity_pseudotime_artifacts.sh --all
```

## 边界约束
- 不修改 `run_end_to_end.sh`
- 不修改 `downstream_scanpy/` 现有执行路径
- 不输出到仓库根目录 `output/`
