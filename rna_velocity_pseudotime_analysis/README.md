# RNA Velocity and Pseudotime Analysis

本模块用于 RNA velocity 与 pseudotime 的探索性分析，属于可选分析步骤，非生产必需。

## Scope
- 不接入现有端到端执行路径
- 不修改主流程配置
- 不写入正式交付目录
- 仅在本目录内生成运行产物

## 从哪个阶段开始导入数据
- 推荐入口阶段：`Scanpy QC 后、聚类前` 的表达矩阵（等价于 `after_qc.h5ad` 语义）
- 次优入口阶段：`Cell Ranger filtered_feature_bc_matrix`（需要在本模块内补充标准化/邻接图等前处理）
- 不推荐直接入口：`processed.h5ad` 最终文件（若已做强降维/子集化，可能影响 velocity 动力学拟合）

## 为什么推荐 QC 后数据
- RNA velocity 需要尽可能保留真实转录动态信号，过度筛选或过早固定低维结构会损失信息
- pseudotime 需要稳定的细胞拓扑关系，QC 后数据噪声更可控
- 在该阶段进入，既能避免低质量细胞干扰，又不会丢失关键过渡态细胞

## 数据阶段与用途映射
| 阶段 | 典型来源 | 适合分析 | 说明 |
|---|---|---|---|
| Stage A: Cell Ranger 输出 | `<run_id>/outs/filtered_feature_bc_matrix` | 基础重建 | 需要在本模块内完成标准化与图构建 |
| Stage B: Scanpy QC 后数据 | `output/<project>/data/after_qc.h5ad` | RNA velocity + pseudotime（推荐） | 质量控制完成，保留动力学信息较好 |
| Stage C: Scanpy 最终数据 | `output/<project>/data/processed.h5ad` | 复核与可视化补充 | 若已高度处理，可能不利于动力学建模 |

## 最小启动流程
1. 先完成主流程到可复用数据阶段（建议拿到 `after_qc.h5ad` 或 `filtered_feature_bc_matrix`）
2. 创建并激活独立环境
3. 在本模块目录中运行独立脚本
4. 检查 `runtime/results/<run_id>/` 与 `runtime/reports/<run_id>/` 输出

## 建议输入准备清单
- 必需：细胞 x 基因表达矩阵（10x MTX 或 h5ad）
- velocity 必需：spliced / unspliced 信息（通常来自 loom 或等价层）
- 建议：批次信息、样本分组、已知起始细胞群标签（用于 pseudotime root 设定）
- 可选：先验 marker 列表（用于结果解释与生物学验证）

## 配置输入项（推荐）
- `input.h5ad`：优先使用，建议指向 `after_qc.h5ad`
- `input.tenx_dir`：备用入口，指向 `filtered_feature_bc_matrix`
- `input.loom`：可选，提供 spliced/unspliced 用于 velocity
- `input.preferred_stage`：建议保持 `after_qc`

## 从BAM生成loom（spliced/unspliced）
```bash
bash rna_velocity_pseudotime_analysis/scripts/generate_velocity_loom.sh \
  --samplefolder /home/zerlinshen/singlecell_factory/pbmc_1k_v3_count \
  --gtf /home/zerlinshen/singlecell_factory/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
  --threads 4
```
- 生成loom默认路径：`<samplefolder>/velocyto/<samplefolder_basename>.loom`
- 校验参数文件：`rna_velocity_pseudotime_analysis/runtime/velocity_inputs/loom_generation_<tag>.json`
- 运行日志文件：`rna_velocity_pseudotime_analysis/runtime/logs/loom_generation_<tag>.log`

## Directory Contract
- `configs/`: 分析配置模板
- `envs/`: 独立运行环境定义
- `scripts/`: 独立执行与清理脚本
- `runtime/`: 运行日志、报告、缓存、结果
- `ci/`: CI/CD 排除策略模板
- `src/`: 分析代码
- `docs/`: 使用与迁移文档

## Isolation Rules
- 所有输出必须写入 `rna_velocity_pseudotime_analysis/runtime/`
- 所有缓存必须写入 `rna_velocity_pseudotime_analysis/runtime/cache/`
- 所有日志必须写入 `rna_velocity_pseudotime_analysis/runtime/logs/`
- 不调用 `run_end_to_end.sh`
- 不修改 `downstream_scanpy/` 下现有主流程脚本

## 如何开始（开箱即用）
1. 准备独立环境
```bash
conda env create -f rna_velocity_pseudotime_analysis/envs/rna_velocity_pseudotime.environment.yml
conda activate sc_rna_velocity_pseudotime
```
2. 复制并修改配置
```bash
cp rna_velocity_pseudotime_analysis/configs/rna_velocity_pseudotime.yaml /tmp/rvpt.yaml
```
3. 执行分析
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3 --config /tmp/rvpt.yaml
```
或自动准备独立环境后再运行：
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3 --prepare-env
```
4. 查看结果
- 日志：`rna_velocity_pseudotime_analysis/runtime/logs/`
- 报告：`rna_velocity_pseudotime_analysis/runtime/reports/<run_id>/`
- 结果：`rna_velocity_pseudotime_analysis/runtime/results/<run_id>/`

## 结果呈现层输出清单
- `figures/umap_leiden.png`：基础细胞图谱
- `figures/pseudotime_dpt_umap.png`：pseudotime 着色图
- `figures/velocity_stream_umap.png`：velocity 流线图（需 spliced/unspliced）
- `figures/velocity_grid_umap.png`：velocity 网格图（需 spliced/unspliced）
- `figures/cellrank_velocity_projection_umap.png`：CellRank 投影图（可选）
- `tables/module_status.csv`：各模块执行状态
- `tables/pseudotime_dpt_by_cell.csv`：单细胞 pseudotime
- `tables/pseudotime_dpt_by_cluster.csv`：按 cluster 汇总的 pseudotime
- `tables/velocity_speed_by_cluster.csv`：按 cluster 汇总的 velocity speed
- `integrated.h5ad`：合并后的分析对象

## 参数与后端选择建议
- `modules.velocity.backend=scvelo`：Python 生态主流 velocity 路径
- `modules.pseudotime.backend=cellrank`：与 scVelo/Scanpy 对接自然
- `modules.r_bridge.backend=monocle3`：需要 R 侧轨迹对照时启用
- `modules.r_bridge.enabled=false`：默认关闭跨语言桥接以减少环境复杂度

## 结果解读最小标准
- velocity：检查流场方向是否与已知发育方向一致
- pseudotime：检查排序是否符合 marker 基因表达梯度
- 一致性：对同一细胞群比较 CellRank 与 Monocle3 趋势是否同向
- 质量门槛：若出现大面积逆向流或断裂轨迹，优先回溯输入阶段与 root 设定

## Run
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3
```
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh --dataset pbmc_1k_v3 --env-name sc_rna_velocity_pseudotime
```
```bash
bash rna_velocity_pseudotime_analysis/scripts/run_rna_velocity_pseudotime_analysis.sh \
  --dataset pbmc_1k_v3 \
  --config rna_velocity_pseudotime_analysis/configs/rna_velocity_pseudotime.pbmc_loom.yaml \
  --python-bin /home/zerlinshen/.conda/envs/sc_rna_velocity_pseudotime/bin/python
```

## Cleanup
```bash
bash rna_velocity_pseudotime_analysis/scripts/cleanup_rna_velocity_pseudotime_artifacts.sh --all
```

## CI/CD Exclusion
请将 `rna_velocity_pseudotime_analysis/**` 设为流水线排除路径，参考 `ci/ci_exclude_paths.txt` 与 `ci/github_paths_ignore_example.yml`。

## Backward Compatibility
兼容入口保留在 `phase1_optional/`，详见 [MIGRATION.md](file:///home/zerlinshen/singlecell_factory/rna_velocity_pseudotime_analysis/docs/MIGRATION.md)。
