# LUSC RNA Velocity 真实运行报告

日期：2026-04-05  
样本：`data/raw/lung_carcinoma_3k_count`  
最终成功运行目录：`/home/zerlinshen/singlecell_factory/results/LUSC_RNA_velocity_REAL_v5_20260405_20260405_020515`

## 1. 任务目标

本次任务是对 LUSC 数据真实执行 `rna_velocity` 模块，并输出：

1. 真实 velocity 图像
2. velocity 量化表格
3. 结果解读与参数调整建议

## 2. 运行配置

执行命令：

```bash
MPLCONFIGDIR=/tmp/matplotlib NUMBA_CACHE_DIR=/tmp/numba_cache PYTHONPATH=. python -m workflow.modular.cli \
  --project LUSC_RNA_velocity_REAL_v5_20260405 \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,rna_velocity \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir /home/zerlinshen/singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A \
  --velocity-n-jobs 8 \
  --parallel-workers 2
```

关键说明：

1. `--velocity-gtf` 未显式提供，系统已自动从 `--transcriptome-dir` 解析到：
   `/home/zerlinshen/singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz`
2. 模块状态全通过（`rna_velocity = ok`）。

## 3. 产物清单

`rna_velocity/` 目录下共 5 个文件：

1. `velocity_stream_umap.png`
2. `velocity_grid_umap.png`
3. `velocity_length_distribution.png`
4. `velocity_confidence.csv`
5. `velocity_top_genes.csv`

## 4. Figure 展示与解读

## 4.1 Stream 图（主图）

![velocity_stream_umap](/home/zerlinshen/singlecell_factory/results/LUSC_RNA_velocity_REAL_v5_20260405_20260405_020515/rna_velocity/velocity_stream_umap.png)

解读：

1. 主体簇群（0/1/6/8 附近）内部存在连续流向，说明局部状态转换较平滑。
2. 右侧簇群（2/3/11）内部与相邻簇之间可见方向性，提示该区域存在活跃转录动态。
3. 孤立簇（如 5、12）主要表现为簇内流向，跨簇联系较弱，可能代表相对独立状态。

## 4.2 Grid 图（局部方向细节）

![velocity_grid_umap](/home/zerlinshen/singlecell_factory/results/LUSC_RNA_velocity_REAL_v5_20260405_20260405_020515/rna_velocity/velocity_grid_umap.png)

解读：

1. 与 stream 图一致，主体区箭头方向存在一致性，不是纯随机噪声。
2. 少数小簇方向较分散，建议结合注释与 DE 结果判断是否为稀有群或技术噪声。

## 4.3 Velocity Length 分布（QC）

![velocity_length_distribution](/home/zerlinshen/singlecell_factory/results/LUSC_RNA_velocity_REAL_v5_20260405_20260405_020515/rna_velocity/velocity_length_distribution.png)

解读：

1. 分布右偏，主峰在低到中等长度区间，符合常见 velocity 长度分布。
2. 中位数约 `4.79`，并有少量长尾细胞（>15），提示存在高动态子群。

## 5. 量化结果摘要

来自 `run_manifest.json` + `velocity_confidence.csv`：

1. 细胞数：`2377`
2. `velocity_mean_confidence`：`0.8612`
3. `confidence median`：`0.8808`
4. `confidence p10 / p90`：`0.7454 / 0.9636`
5. `confidence < 0.5`：`1.60%`
6. `velocity length mean / median`：`5.5049 / 4.79`
7. `length p10 / p90`：`3.01 / 8.284`
8. `corr(confidence, length)`：`0.136`（弱正相关）

整体判断：本次 velocity 结果质量较好，置信度整体偏高，低置信细胞占比低。

## 6. 按簇的重点观察

按簇平均 `velocity_confidence`：

高置信 Top 5：

1. Leiden 3：`0.959`（n=92）
2. Leiden 0：`0.926`（n=406）
3. Leiden 10：`0.921`（n=31）
4. Leiden 11：`0.916`（n=51）
5. Leiden 2：`0.887`（n=328）

低置信 Bottom 5：

1. Leiden 13：`0.724`（n=48）
2. Leiden 9：`0.765`（n=67）
3. Leiden 12：`0.807`（n=81）
4. Leiden 5：`0.818`（n=155）
5. Leiden 4：`0.819`（n=145）

建议优先复核低置信簇（13/9/12）：结合 marker、DE、细胞周期和双细胞概率判断是否混合群。

## 7. Top Velocity Genes（摘要）

`velocity_top_genes.csv` 中重复出现较多的基因包括：

1. `ST6GAL1`
2. `ZFYVE16`
3. `HDAC9`
4. `ITGAX`
5. `PLCG2`
6. `LEF1`

建议：这些基因先作为“候选动态基因”，再与 cluster 注释和已知通路交叉验证，不直接作为最终生物结论。

## 8. 针对本次结果的调参建议

如果你要进一步优化本样本的 velocity，可按以下顺序迭代：

1. 保持 `stochastic` 先做稳定版本（当前质量已可用）。
2. 对低置信簇（13/9/12）尝试：
   - `--velocity-min-shared-counts 20 -> 15`（保留更多基因）
   - 或 `--velocity-n-neighbors 30 -> 40`（增强局部平滑）
3. 若要更强时间结构解释，再跑 `--velocity-mode dynamical`，并新增查看：
   - `velocity_latent_time_umap.png`
   - `velocity_phase_portraits.png`

## 9. 结论

本次 LUSC RNA velocity 真实运行已成功完成，且结果具备可解释性：

1. 模块执行成功，输入链路（BAM + 自动 GTF）打通。
2. velocity 图存在稳定方向结构，不是纯噪声。
3. 置信度指标整体较高，可进入后续生物学解释阶段。

