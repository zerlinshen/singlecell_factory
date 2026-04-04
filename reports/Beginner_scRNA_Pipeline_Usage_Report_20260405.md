# 单细胞分析流程使用报告（初学者版）

日期：2026-04-05  
项目路径：`/home/zerlinshen/singlecell_factory`

## 1. 这套流程做什么

这是一套面向 10x 单细胞转录组的自动化分析流程，核心目标是：

1. 从 Cell Ranger 矩阵开始完成标准 scRNA 分析
2. 自动完成质量控制、聚类、注释、差异分析、轨迹等模块
3. 在可选条件下完成 RNA velocity（真实剪接动力学）
4. 生成标准化图表和表格，便于复现和汇报

你可以把它理解为“开箱可跑”的模块化分析框架。

## 2. 你需要准备什么

## 2.1 必需输入

1. `sample-root` 目录，至少包含：
- `outs/filtered_feature_bc_matrix`（Cell Ranger count 输出）

2. Python 环境中已安装项目依赖

## 2.2 RNA velocity 额外输入（可选模块）

二选一：

1. `--velocity-loom`（已有 loom）
2. `--velocity-bam` + reference（推荐）

现在流程已支持自动找 GTF：

- 你给 `--transcriptome-dir` 后，会自动从 reference 中查找 `genes.gtf(.gz)`
- 典型位置：`<transcriptome-dir>/genes/genes.gtf.gz`

## 3. 与 reference 的正确关系（非常关键）

结论：**Cell Ranger 与 RNA velocity 必须使用同一套 reference。**

1. `cellranger count --transcriptome` 用哪套 reference
2. RNA velocity 就从同一 reference 取 GTF（自动或手动）

这样可以保证：

1. 染色体命名一致
2. 基因注释版本一致
3. 坐标体系一致

## 4. 如何运行

## 4.1 最小运行（先看数据质量和聚类）

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. python -m workflow.modular.cli \
  --project DEMO_MINI \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering
```

## 4.2 推荐全流程（不含 RNA velocity）

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. python -m workflow.modular.cli \
  --project DEMO_FULL \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,differential_expression,annotation,trajectory,pseudo_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution \
  --parallel-workers 4
```

## 4.3 加入 RNA velocity（BAM + 自动 GTF）

```bash
NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. python -m workflow.modular.cli \
  --project DEMO_VELOCITY \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,rna_velocity \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir /home/zerlinshen/singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A
```

## 5. 结果先看哪些文件

每次运行都会生成一个时间戳目录，重点看：

1. `run_manifest.json`：全局元数据、模块顺序、关键统计
2. `module_status.csv`：每个模块是否成功
3. `final_adata.h5ad`：后续二次分析主文件

建议“先判是否可用，再做生物解释”：

1. 先看 `module_status.csv` 是否全是 `ok`
2. 再看 QC 图是否合理
3. 再看聚类、注释、DE、轨迹

## 6. 初学者最重要：如何看 QC 并调参数

QC 没有“万能固定阈值”，要按你的组织、测序深度、细胞类型分布调整。下面给你可执行规则。

## 6.1 常用 QC 参数（当前默认）

1. `--min-genes 200`
2. `--max-genes 7000`
3. `--min-counts 500`
4. `--max-counts 50000`
5. `--max-mito-pct 20`
6. `--max-ribo-pct 50`

## 6.2 如何根据图判断要不要改

看这 4 张图：

1. `qc_violin_pre_filter.png`
2. `qc_scatter_pre_filter.png`
3. `qc_violin_post_filter.png`
4. `qc_scatter_post_filter.png`

判断逻辑：

1. 如果大量细胞 `n_genes_by_counts` 很低并挤在底部：
- 提高 `--min-genes`（例如 200 -> 300）

2. 如果有明显“超高 UMI + 超高基因数”的长尾：
- 降低 `--max-counts` 或 `--max-genes`
- 这类常是潜在 doublet/multiplet

3. 如果线粒体比例整体偏高：
- 降低 `--max-mito-pct`（如 20 -> 15 或 10）
- 但肿瘤/坏死样本可能本来高，不要一刀切

4. 如果过滤后剩余细胞太少（例如只剩原始的 <40%）：
- 阈值过严，放宽 `min/max genes` 或 `max-mito-pct`

## 6.3 一个稳健调参顺序

1. 先固定 `min/max counts`
2. 再调 `min/max genes`
3. 最后调 `max-mito-pct`
4. 每次只改 1-2 个参数并重跑

## 7. Doublet 怎么看与怎么调

关注：

1. `doublet_scores.png`
2. `run_manifest.json` 中 `doublets_detected`、`doublet_rate_pct`

建议：

1. 肿瘤组织常见 doublet 率可高于 PBMC
2. 如果你怀疑检出过少，可提高 `--expected-doublet-rate`（如 0.06 -> 0.08）
3. 如果你担心误删，先加 `--no-remove-doublets` 做对照

当前流程对极小数据集有 fallback 保护：

1. Scrublet 不稳定时不再崩溃
2. 自动标记为 singlet 并记录 `doublet_method=fallback_all_singlets`

## 8. 聚类参数怎么调

关键参数：

1. `--n-pcs`
2. `--n-neighbors`
3. `--leiden-resolution`

经验规则：

1. cluster 太碎：降低 `--leiden-resolution`（0.8 -> 0.5）
2. cluster 太粗：提高 `--leiden-resolution`（0.8 -> 1.0/1.2）
3. 稳定性差：适度提高 `n_neighbors`（15 -> 20/30）
4. 结构复杂：增加 `n-pcs`（40 -> 50）

## 9. RNA velocity 参数怎么调

关键参数：

1. `--velocity-mode stochastic|dynamical`
2. `--velocity-min-shared-counts`
3. `--velocity-n-pcs`
4. `--velocity-n-neighbors`
5. `--velocity-n-jobs`

建议：

1. 新手先用 `stochastic`，更稳更快
2. 动力学足够好再尝试 `dynamical`
3. 若基因过滤过多，降低 `--velocity-min-shared-counts`（20 -> 10）
4. 大样本可增加 `--velocity-n-jobs`

如何判断 velocity 质量：

1. `velocity_confidence.csv` 均值不能过低
2. `velocity_stream_umap.png` 流向是否与已知生物过程一致
3. `velocity_length_distribution.png` 是否过于集中在接近 0
4. dynamical 模式下 `velocity_latent_time_umap.png` 是否有合理连续梯度

## 10. 人类与小鼠都能跑吗

可以，但要区分“技术可运行”和“生物注释准确性”。

1. 技术层面：
- 只要 reference 一致，Cell Ranger + clustering + trajectory + velocity 都能跑

2. 生物注释层面：
- 部分 marker/注释策略偏人类资源
- 小鼠建议替换 marker 集和物种相关资源后再做解释

## 11. 失败时先排查这几件事

1. `module_status.csv` 哪个模块 first fail
2. `run_manifest.json` 的 metadata 是否异常（细胞数骤降、cluster 数异常）
3. 输入路径是否真实存在
4. reference 与 BAM 是否同源
5. 是否在日志里出现“No spliced/unspliced layers found”

## 12. 给初学者的执行建议（强烈推荐）

1. 第一次先跑最小流程（QC + clustering）
2. QC 合理后再加 annotation + DE
3. 轨迹/Pseudo velocity 放在第三轮
4. RNA velocity 放最后，确保 reference/BAM/GTF一致
5. 每轮只改少量参数，保存命令和结果目录

---

如果你愿意，我下一步可以按你的数据目标（比如免疫微环境优先、肿瘤进化优先）给你出一套“参数模板 + 解释模板”，你可以直接复用到后续样本。
