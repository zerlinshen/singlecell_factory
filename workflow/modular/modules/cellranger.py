from __future__ import annotations

from pathlib import Path
import subprocess

import anndata as ad
import pandas as pd

from ._scanpy_compat import has_api, import_scanpy_or_stub, scanpy_import_error
from ._counts_contract import (
    stamp_factory_counts_provenance,
    validate_counts_layer_contract,
)
from .._gene_symbols import normalize_var_to_symbols

sc = import_scanpy_or_stub()

from ..context import PipelineContext


__references__ = {
    "Zheng_10x_2017": {
        "title": "Massively parallel digital transcriptional profiling of single cells",
        "authors": "Zheng et al.",
        "journal": "Nature Communications",
        "year": "2017",
        "doi": "10.1038/ncomms14049",
        "description": "10x Chromium Single Cell 3' chemistry \u2014 describes the output format this module ingests.",
    },
    "CellRanger_software": {
        "title": "Cell Ranger Single Cell Software (v7+)",
        "authors": "10x Genomics",
        "journal": "Software documentation",
        "year": "2024",
        "doi": "https://support.10xgenomics.com/single-cell-gene-expression/software",
        "description": "Authoritative spec for the barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz layout consumed here.",
    },
}



class CellRangerModule:
    """Mandatory module: validate/run Cell Ranger and load matrix."""

    name = "cellranger"
    required = True

    FLEX_16_PROBE_MAP = {
        "ACTTTAGG": "BC001",
        "AACGGGAA": "BC002",
        "AGTAGGCT": "BC003",
        "ATGTTGAC": "BC004",
        "ACAGACCT": "BC005",
        "ATCCCAAC": "BC006",
        "AAGTAGAG": "BC007",
        "AGCTGTGA": "BC008",
        "ACAGTCTG": "BC009",
        "AGTGAGTG": "BC010",
        "AGAGGCAA": "BC011",
        "ACTACTCA": "BC012",
        "ATACGTCA": "BC013",
        "ATCATGTG": "BC014",
        "AACGCCGA": "BC015",
        "ATTCGGTT": "BC016",
    }

    @staticmethod
    def _needs_counts_layer(ctx: PipelineContext) -> bool:
        return "pseudobulk_de" in set(ctx.cfg.optional_modules)

    @staticmethod
    def _should_lazy_read(zarr_path: Path, ctx: PipelineContext) -> bool:
        """Decide whether to use lazy zarr loading.

        Priority: SC_LAZY_READ env var > cfg.lazy_read > auto (file size > 5 GB).
        """
        import os as _os
        flag = (
            _os.environ.get("SC_LAZY_READ", "").strip().lower()
            or getattr(ctx.cfg, "lazy_read", "auto")
        )
        if flag == "true":
            return True
        if flag == "false":
            return False
        # "auto": trigger lazy read when the zarr store is larger than 5 GB
        try:
            total = sum(f.stat().st_size for f in zarr_path.rglob("*") if f.is_file())
            return total > 5 * 1024 ** 3
        except Exception:
            return False

    def _annotate_flex_probe_groups(self, adata, sample_root: Path, ctx: PipelineContext) -> None:
        if adata.n_obs == 0:
            return
        if "sample" in adata.obs.columns and pd.Series(adata.obs["sample"]).nunique() > 1:
            return
        suffixes = pd.Index(adata.obs_names.astype(str)).str.split("-").str[0].str[-8:]
        mapped = suffixes.map(self.FLEX_16_PROBE_MAP)
        nunique = mapped.nunique(dropna=True)
        if nunique < 2:
            return
        adata.obs["probe_barcode_seq"] = suffixes.values
        adata.obs["probe_barcode_id"] = mapped.fillna("UNKNOWN").values
        adata.obs["sample"] = adata.obs["probe_barcode_id"].astype(str).values
        adata.obs["sample_label_strategy"] = "probe_barcode_16plex_proxy"
        agg = sample_root / "16plex_900k_32_NSCLC_multiplex_aggregation.csv"
        if agg.exists():
            try:
                meta = pd.read_csv(agg)
                if {"description", "sample_id"}.issubset(meta.columns):
                    meta["probe_barcode_id"] = meta["description"].astype(str).str.extract(r"_(BC\d+)$", expand=False)
                    per_bc = meta.dropna(subset=["probe_barcode_id"]).groupby("probe_barcode_id")["description"].apply(lambda s: "|".join(sorted(set(map(str, s)))))
                    adata.obs["probe_barcode_description_pool"] = adata.obs["probe_barcode_id"].map(per_bc).fillna("")
            except Exception as exc:
                ctx.metadata["flex_probe_annotation_warning"] = str(exc)
        ctx.metadata["sample_label_strategy"] = "probe_barcode_16plex_proxy"
        ctx.metadata["sample_label_nunique"] = int(pd.Series(adata.obs["sample"]).nunique())

    @staticmethod
    def _apply_cohort_subset(adata, ctx: PipelineContext):
        """Filter adata to cohort subset defined by cfg.cohort_subset predicates."""
        import logging as _logging
        _log = _logging.getLogger(__name__)
        cohort_subset = getattr(ctx.cfg, "cohort_subset", None)
        if not cohort_subset:
            return adata
        if isinstance(cohort_subset, str):
            cohort_subset = [cohort_subset]
        mask = None
        for spec in cohort_subset:
            if "=" not in spec:
                _log.warning("cohort_subset spec '%s' has no '='; skipping.", spec)
                continue
            col, vals_str = spec.split("=", 1)
            col = col.strip()
            vals = [v.strip() for v in vals_str.split(",") if v.strip()]
            if col not in adata.obs.columns:
                raise ValueError(f"cohort_subset column '{col}' not found in adata.obs.")
            submask = adata.obs[col].isin(vals).values
            mask = submask if mask is None else (mask & submask)
        if mask is None:
            return adata
        n_before = int(adata.n_obs)
        adata = adata[mask].copy()
        n_after = int(adata.n_obs)
        subset_str = "; ".join(cohort_subset) if isinstance(cohort_subset, list) else str(cohort_subset)
        _log.info(
            "Cohort subset applied: %s, n_obs reduced from %d to %d",
            subset_str, n_before, n_after,
        )
        ctx.metadata["cohort_subset_applied"] = subset_str
        ctx.metadata["cohort_subset_n_obs_before"] = n_before
        ctx.metadata["cohort_subset_n_obs_after"] = n_after
        return adata

    def run(self, ctx: PipelineContext) -> None:
        cfg = ctx.cfg.cellranger
        sample_root = cfg.sample_root
        direct_h5ad = cfg.input_h5ad
        prepared_h5ad = (
            Path(direct_h5ad)
            if direct_h5ad is not None
            else sample_root / "prepared_input.h5ad"
        )
        prepared_zarr = sample_root / "prepared_input.zarr"
        outs = cfg.outs_dir

        if direct_h5ad is not None:
            if prepared_h5ad.suffix.lower() != ".h5ad":
                raise ValueError(
                    f"direct input must use the .h5ad suffix: {prepared_h5ad}"
                )
            if not prepared_h5ad.is_file():
                raise FileNotFoundError(f"direct input h5ad not found: {prepared_h5ad}")

        if direct_h5ad is not None or prepared_h5ad.exists() or prepared_zarr.exists():
            if direct_h5ad is not None or prepared_h5ad.exists():
                adata = ad.read_h5ad(prepared_h5ad)
                ctx.metadata["prepared_input_source"] = str(prepared_h5ad)
                ctx.metadata["prepared_input_loading_mode"] = "eager_h5ad"
                ctx.metadata["direct_input_h5ad"] = direct_h5ad is not None
            else:
                if self._should_lazy_read(prepared_zarr, ctx) and hasattr(ad.experimental, "read_lazy"):
                    adata = ad.experimental.read_lazy(prepared_zarr)
                    ctx.metadata["prepared_input_loading_mode"] = "lazy_zarr"
                else:
                    adata = ad.read_zarr(prepared_zarr)
                    ctx.metadata["prepared_input_loading_mode"] = "eager_zarr"
                ctx.metadata["prepared_input_source"] = str(prepared_zarr)
            adata.var_names_make_unique()
            # Establish the same gene-identifier contract the Cell Ranger path
            # gets from read_10x_mtx(var_names="gene_symbols"). Public atlases
            # distributed under the CZ CELLxGENE schema are Ensembl-indexed, and
            # every symbol-keyed stage downstream (mito/ribo/hb QC flags, marker
            # annotation, signature scoring) silently matches nothing otherwise.
            ctx.metadata["gene_namespace"] = normalize_var_to_symbols(adata)
            self._annotate_flex_probe_groups(adata, sample_root, ctx)
            adata = self._apply_cohort_subset(adata, ctx)
            if cfg.sample_id and cfg.sample_id != "lusc":
                if "sample" not in adata.obs.columns:
                    raise ValueError(
                        f"--sample-id={cfg.sample_id!r} given but adata.obs has no 'sample' column; "
                        f"use --cohort-subset instead"
                    )
                sample_ids = [s.strip() for s in cfg.sample_id.split(",")]
                mask = adata.obs["sample"].astype(str).isin(sample_ids)
                n_before = adata.n_obs
                adata = adata[mask].copy()
                ctx.metadata["sample_id_filter_applied"] = sample_ids
                ctx.metadata["sample_id_filter_n_before"] = int(n_before)
                ctx.metadata["sample_id_filter_n_after"] = int(adata.n_obs)
            if self._needs_counts_layer(ctx):
                if "counts" in adata.layers:
                    contract = validate_counts_layer_contract(adata)
                    ctx.metadata["counts_layer_contract"] = contract
                    ctx.metadata["counts_layer_preserved"] = True
                    ctx.metadata["counts_provenance_status"] = (
                        "validated_explicit_counts_layer"
                    )
                else:
                    # An arbitrary/prepared AnnData X may be log-normalized. Never
                    # manufacture a raw-count claim by copying it into a layer whose
                    # name downstream code trusts. Confirmatory pseudobulk will fail
                    # closed with the missing-count contract; exploratory runs skip.
                    ctx.metadata["counts_layer_preserved"] = False
                    ctx.metadata["counts_provenance_status"] = (
                        "missing_explicit_counts_layer"
                    )
            else:
                ctx.metadata["counts_layer_preserved"] = "counts" in adata.layers
                ctx.metadata["counts_provenance_status"] = "not_requested"
            ctx.adata = adata
            ctx.metadata["prepared_input_h5ad"] = (
                str(prepared_h5ad)
                if direct_h5ad is not None or prepared_h5ad.exists()
                else str(prepared_zarr)
            )
            ctx.metadata["raw_cells"] = int(adata.n_obs)
            ctx.metadata["raw_genes"] = int(adata.n_vars)
            return

        # If user requests force run or no previous Cell Ranger output exists, run `cellranger count`.
        if cfg.force_run or not outs.exists():
            if not cfg.run_if_missing:
                raise FileNotFoundError(f"Cell Ranger output not found: {outs}")
            self._run_cellranger(
                cfg.sample_root,
                cfg.fastq_dir,
                cfg.transcriptome_dir,
                cfg.sample_id,
                cfg.localcores,
                cfg.localmem,
            )

        if not outs.exists():
            raise FileNotFoundError(f"Cell Ranger output still missing after attempt: {outs}")

        if not has_api(sc, "read_10x_mtx"):
            err = scanpy_import_error(sc)
            raise ImportError(f"scanpy.read_10x_mtx is unavailable: {err}") from err
        adata = sc.read_10x_mtx(str(outs), var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata = self._apply_cohort_subset(adata, ctx)
        # Preserve raw UMI counts only when downstream modules require them
        # (notably pseudobulk_de). This avoids a large duplicate matrix for
        # massive clustering-first runs.
        if self._needs_counts_layer(ctx):
            adata.layers["counts"] = adata.X.copy()
            stamp_factory_counts_provenance(
                adata, source="cellranger_raw_feature_bc_matrix"
            )
            ctx.metadata["counts_layer_contract"] = validate_counts_layer_contract(
                adata
            )
            ctx.metadata["counts_provenance_status"] = "validated_factory_reader"
        else:
            ctx.metadata["counts_provenance_status"] = "not_requested"
        ctx.metadata["counts_layer_preserved"] = bool(
            self._needs_counts_layer(ctx) and "counts" in adata.layers
        )
        ctx.adata = adata
        ctx.metadata["cellranger_outs"] = str(outs)
        ctx.metadata["raw_cells"] = int(adata.n_obs)
        ctx.metadata["raw_genes"] = int(adata.n_vars)

    @staticmethod
    def _run_cellranger(
        sample_root: Path,
        fastq_dir: Path | None,
        transcriptome_dir: Path | None,
        sample_id: str,
        localcores: int,
        localmem: int,
    ) -> None:
        if fastq_dir is None or transcriptome_dir is None:
            raise ValueError("Running Cell Ranger requires fastq_dir and transcriptome_dir.")
        cmd = [
            "cellranger",
            "count",
            f"--id={sample_id}",
            "--create-bam=true",
            f"--fastqs={fastq_dir}",
            f"--transcriptome={transcriptome_dir}",
            f"--localcores={localcores}",
            f"--localmem={localmem}",
        ]
        subprocess.run(cmd, check=True, cwd=sample_root)
