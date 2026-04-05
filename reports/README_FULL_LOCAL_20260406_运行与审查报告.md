# README_FULL_LOCAL_20260406 运行与审查报告

## 1) 执行元信息
- 操作者: Codex (GPT-5) 按用户指令执行
- 执行日期: 2026-04-06 (Asia/Shanghai)
- 目标: 按 README 最新“Full local analysis”配置完成真实全流程运行，并验证模块状态、关键产物与可复现性
- 代码仓库: `/home/zerlinshen/singlecell_factory`

## 2) 实际执行命令
```bash
PYTHONPATH=. MPLCONFIGDIR=$PWD/.mplconfig NUMBA_CACHE_DIR=/tmp/numba_cache \
python -m workflow.modular.cli \
  --project README_FULL_LOCAL_20260406 \
  --sample-root data/raw/lung_carcinoma_3k_count \
  --optional-modules clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,rna_velocity,cnv_inference,pathway_analysis,cell_communication,gene_regulatory_network,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,evolution,pseudobulk_de,cell_fate,composition,metacell \
  --velocity-bam data/raw/lung_carcinoma_3k_count/outs/possorted_genome_bam.bam \
  --transcriptome-dir ref/reference/refdata-gex-GRCh38-2024-A
```

## 3) 运行产物位置
- 运行目录: `results/README_FULL_LOCAL_20260406_20260406_012249`
- 状态文件: `results/README_FULL_LOCAL_20260406_20260406_012249/module_status.csv`
- 清单文件: `results/README_FULL_LOCAL_20260406_20260406_012249/run_manifest.json`

## 4) 结果概览
- 模块状态: **23/23 ok**（mandatory 3 + optional 20）
- Pipeline 总耗时: **40.962 s**
- 细胞变化: `2588 -> 2382 (QC) -> 2377 (doublet removal)`
- 聚类数: `n_clusters = 18`
- GPU后端: `clustering_backend=cpu`, `de_backend=cpu`（当前主机未检测到可用 CUDA 后端）
- RNA velocity: **成功**
  - `velocity_gtf_resolved = ref/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz`
  - `velocity_extract_n_jobs = 8`
  - `velocity_extract_seconds = 0.093`
  - `velocity_mean_confidence = 0.86117`
- 单样本相关预期行为:
  - `batch_correction_status = skipped_missing_batch_key`
  - `pseudobulk_de_status = skipped_missing_grouping_columns`

## 5) 关键输出核验
- 聚类 UMAP: `clustering/umap_leiden.png`
- RNA velocity 主图: `rna_velocity/velocity_stream_umap.png`
- 进化模块克隆 marker: `evolution/evolution_clone_markers.csv`
- 细胞命运热图: `cell_fate/fate_heatmap.png`
- 组成差异统计: `composition/composition_test_results.csv`
- 元细胞对象: `metacell/metacells.h5ad`

## 6) 本次审查结论
- README 中更新的本地全流程命令可直接复现并跑通。
- `rna_velocity` 输入链路（BAM + transcriptome-dir 自动发现 GTF）验证通过。
- 该数据为单样本，`batch_correction` 与 `pseudobulk_de` 的 skip 属设计内行为，不是失败。

## 7) 建议的下游分析
1. 多样本扩展后启用 `pseudobulk_de` 做条件比较（治疗前后/亚型间），避免单细胞层面伪重复。
2. 针对 `evolution` 的 3 个 clone 做通路差异富集，优先比较 Clone_1 vs Clone_3 的 EMT/IFN/代谢轴。
3. 联动 `cell_fate` 与 `rna_velocity`：筛选 fate 高概率且 velocity 长度高的过渡群，定位驱动基因。
4. 在 `immune_phenotyping` 中追踪 `CD8_exhausted` 与 checkpoint 表达（PDCD1/LAG3/TIGIT）的空间与轨迹位置。

## 8) 可能的生物学故事（LUSC 当前数据）
- 故事线A（克隆进化）: CNV 高分群体构成 3 个克隆分支，其中主干克隆（Clone_1）占比最高，可能代表优势扩增亚群；次级克隆（Clone_2/3）对应不同微环境适应策略。
- 故事线B（免疫压力与逃逸）: IFN-gamma/TIS 通路富集与 exhaustion 亚群共现，提示“免疫激活-免疫抑制并存”的动态平衡。
- 故事线C（转移潜能）: EMT 与 invasion signature 在部分群体升高，并与 pseudotime 末端状态重叠，支持“上皮-间质转化驱动进展”的假设。

> 说明: 上述故事为可检验假设，不等同于最终结论；建议结合临床分组、多样本重复与外部队列验证。
