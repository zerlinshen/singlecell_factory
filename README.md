# singlecell_factory v1.0.0

本仓库提供两条可独立运行的单细胞分析流程：

- `workflow/standard.py`：Cell Ranger `outs/filtered_feature_bc_matrix` 到 Scanpy 聚类与UMAP。
- `workflow/velocity.py`：loom 生成、RNA velocity 估计与速度场可视化。

两条流程共享 `data/` 与 `ref/`，互不依赖。

## 目录规范

```text
singlecell_factory/
├── data/
│   └── raw/
│       ├── fastq/
│       ├── *_count/
│       │   └── outs/
│       └── cellranger_runs/
├── ref/
│   └── reference/
├── workflow/
│   ├── standard.py
│   └── velocity.py
├── rna_velocity_pseudotime_analysis/
├── tests/
├── environment.yml
├── pyproject.toml
└── .pre-commit-config.yaml
```

## 一键环境

```bash
cd /home/zerlinshen/singlecell_factory
conda env create -f environment.yml
conda activate sc10x
pre-commit install
```

## 测试数据下载

PBMC 1k:

```bash
mkdir -p /home/zerlinshen/singlecell_factory/data/raw/fastq
cd /home/zerlinshen/singlecell_factory/data/raw/fastq
wget -O pbmc_1k_v3_fastqs.tar https://cf.10xgenomics.com/samples/cell-exp/7.1.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
```

Lung carcinoma:

```bash
cd /home/zerlinshen/singlecell_factory/data/raw/fastq
wget -O Chromium_5p_GEX_Human_Lung_Carcinoma_DTC_fastqs.tar https://cf.10xgenomics.com/samples/cell-vdj/7.0.1/SC5pv2_GEX_Human_Lung_Carcinoma_DTC/SC5pv2_GEX_Human_Lung_Carcinoma_DTC_fastqs.tar
tar -xvf Chromium_5p_GEX_Human_Lung_Carcinoma_DTC_fastqs.tar
```

## 标准流程（copy-paste）

```bash
python -m workflow.standard \
  --project pbmc_standard_v1 \
  --tenx-dir /home/zerlinshen/singlecell_factory/data/raw/pbmc_1k_v3_count/outs/filtered_feature_bc_matrix \
  --output-dir /home/zerlinshen/singlecell_factory/output/workflow_standard
```

预计运行时长：

- 1k cells: 2–5 分钟
- 3k cells: 5–12 分钟
- 50k cells: 35–90 分钟（依赖CPU/磁盘）

## Velocity流程（copy-paste）

```bash
python -m workflow.velocity \
  --samplefolder /home/zerlinshen/singlecell_factory/data/raw/lung_carcinoma_3k_count \
  --gtf-path /home/zerlinshen/singlecell_factory/ref/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
  --tenx-dir /home/zerlinshen/singlecell_factory/data/raw/lung_carcinoma_3k_count/outs/filtered_feature_bc_matrix \
  --run-id lung_velocity_v1 \
  --config-path /home/zerlinshen/singlecell_factory/rna_velocity_pseudotime_analysis/configs/rna_velocity_pseudotime.lung_loom.yaml \
  --python-bin /home/zerlinshen/.conda/envs/sc_rna_velocity_pseudotime/bin/python
```

预计运行时长：

- loom 生成：20–90 分钟
- velocity估计与可视化：5–20 分钟

## 质量保障

```bash
pre-commit run --all-files
pytest
```

覆盖率门槛配置在 `pyproject.toml`：`workflow` 包 `>=90%`。

## 性能目标

- 总时长下降目标：`>=50%`
- 峰值内存下降目标：`>=30%`

建议先跑旧入口（`downstream_scanpy`）记录基线，再跑 `workflow` 入口做对比。
