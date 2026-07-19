from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from ..context import PipelineContext


__references__ = {
    "Squair_pseudobulk_2021": {
        "title": "Confronting false discoveries in single-cell differential expression",
        "authors": "Squair et al.",
        "journal": "Nature Communications",
        "year": "2021",
        "doi": "10.1038/s41467-021-25960-2",
        "description": "Benchmark demonstrating pseudobulk DE outperforms per-cell DE \u2014 methodology applied here.",
    },
}


logger = logging.getLogger(__name__)
MIN_CELLS_PER_SAMPLE = 3
MIN_SAMPLES_PER_CONDITION = 2


class PseudobulkInferenceContractError(ValueError):
    """Machine-classified failure of the biological replicate contract."""

    def __init__(self, inference_status: str, message: str):
        super().__init__(message)
        self.inference_status = inference_status


class PseudobulkDEModule:
    """Pseudobulk DE with an explicit contrast contract.

    Supported modes:
    - confirmatory: explicit contrast(s) across a contrast column
    - exploratory: optional group-vs-rest summaries when explicitly enabled
    """

    name = "pseudobulk_de"

    @staticmethod
    def _record_inference(
        ctx: PipelineContext,
        *,
        inference_class: str,
        inference_status: str,
        claimable: bool,
        **details,
    ) -> None:
        payload = {
            "inference_class": inference_class,
            "inference_status": inference_status,
            "claimable": bool(claimable),
            **details,
        }
        contrast_contract = ctx.metadata.get("pseudobulk_de_contrast_contract")
        if isinstance(contrast_contract, dict):
            payload["contrast_contract"] = json.loads(json.dumps(contrast_contract))
        ctx.metadata["pseudobulk_de_inference_class"] = inference_class
        ctx.metadata["pseudobulk_de_inference_status"] = inference_status
        ctx.metadata["pseudobulk_de_claimable"] = bool(claimable)
        if ctx.adata is not None:
            ctx.adata.uns["pseudobulk_de"] = payload

    def run(self, ctx: PipelineContext) -> None:
        import os
        if os.environ.get("SC_MEM_GUARD", "").lower() == "on":
            from .._mem_guard import MemoryGuard, MemoryGuardError
            try:
                with MemoryGuard(ctx, self.name) as mg:
                    mg.check("entry")
                    return self._run_impl(ctx)
            except MemoryGuardError as exc:
                logger.warning("MemoryGuard: %s", exc)
                ctx.metadata.setdefault("mem_warnings", []).append(
                    {"module": self.name, "error": str(exc)}
                )
        return self._run_impl(ctx)

    def _run_impl(self, ctx: PipelineContext) -> None:
        adata = ctx.adata
        if adata is None:
            raise ValueError("Pseudobulk DE requires AnnData.")

        cfg = getattr(ctx.cfg, "pseudobulk", None)
        min_cells = int(getattr(cfg, "min_cells_per_sample", MIN_CELLS_PER_SAMPLE))
        min_samples = int(getattr(cfg, "min_samples_per_condition", MIN_SAMPLES_PER_CONDITION))
        self._record_inference(
            ctx,
            inference_class="not_tested",
            inference_status="not_started",
            claimable=False,
        )

        confirmatory_requested = bool(
            getattr(cfg, "contrast_json", None)
            or getattr(cfg, "contrast_a", None)
            or getattr(cfg, "contrast_b", None)
        )
        try:
            explicit_contrasts, contrast_col, contract_message = (
                self._resolve_confirmatory_contract(adata, ctx)
            )
        except PseudobulkInferenceContractError as exc:
            ctx.metadata["pseudobulk_de_status"] = "failed_invalid_contrast_contract"
            ctx.metadata["pseudobulk_de_mode"] = "failed"
            self._record_inference(
                ctx,
                inference_class="replicate_aware_pseudobulk",
                inference_status=exc.inference_status,
                claimable=False,
            )
            raise
        except ValueError as exc:
            ctx.metadata["pseudobulk_de_status"] = "failed_invalid_contrast_contract"
            ctx.metadata["pseudobulk_de_mode"] = "failed"
            self._record_inference(
                ctx,
                inference_class="replicate_aware_pseudobulk",
                inference_status="not_testable_invalid_contrast_contract",
                claimable=False,
            )
            raise ValueError(str(exc)) from exc
        if confirmatory_requested and not explicit_contrasts:
            ctx.metadata["pseudobulk_de_status"] = "failed_invalid_contrast_contract"
            ctx.metadata["pseudobulk_de_mode"] = "failed"
            self._record_inference(
                ctx,
                inference_class="replicate_aware_pseudobulk",
                inference_status="not_testable_invalid_contrast_contract",
                claimable=False,
            )
            raise ValueError(contract_message or "Invalid confirmatory pseudobulk contract.")
        exploratory = bool(getattr(cfg, "exploratory_group_vs_rest", False))

        configured_sample_col = getattr(cfg, "sample_col", None)
        if explicit_contrasts and (
            not configured_sample_col or configured_sample_col not in adata.obs.columns
        ):
            msg = (
                "Confirmatory condition contrasts require an explicit biological "
                "--pseudobulk-sample-col that exists in adata.obs; technical batch "
                "auto-detection is not claim-safe."
            )
            ctx.metadata["pseudobulk_de_status"] = "failed_missing_biological_sample_col"
            ctx.metadata["pseudobulk_de_mode"] = "failed"
            self._record_inference(
                ctx,
                inference_class="replicate_aware_pseudobulk",
                inference_status="not_testable_missing_biological_sample_col",
                claimable=False,
                requested_biological_sample_col=configured_sample_col,
            )
            raise ValueError(msg)

        try:
            sample_col = self._resolve_sample_col(
                adata,
                getattr(cfg, "sample_col", None) or ctx.cfg.batch.batch_key,
            )
            group_col = self._resolve_group_col(
                adata,
                getattr(cfg, "group_col", "cell_type"),
            )
        except ValueError as exc:
            logger.warning("%s; skipping pseudobulk DE.", exc)
            ctx.metadata["pseudobulk_de_status"] = (
                "failed_missing_grouping_columns"
                if confirmatory_requested
                else "skipped_missing_grouping_columns"
            )
            ctx.metadata["pseudobulk_de_mode"] = (
                "failed" if confirmatory_requested else "skipped"
            )
            self._record_inference(
                ctx,
                inference_class=(
                    "replicate_aware_pseudobulk" if explicit_contrasts else "not_tested"
                ),
                inference_status="not_testable_missing_grouping_columns",
                claimable=False,
            )
            if confirmatory_requested:
                raise ValueError(str(exc)) from exc
            ctx.status(self.name, "skipped", str(exc))
            return

        if len(adata.obs[sample_col].unique()) < 2:
            msg = f"Only 1 sample in '{sample_col}'; pseudobulk DE needs >= 2. Skipping."
            logger.warning(msg)
            ctx.metadata["pseudobulk_de_status"] = (
                "failed_single_sample" if confirmatory_requested else "skipped_single_sample"
            )
            ctx.metadata["pseudobulk_de_mode"] = (
                "failed" if confirmatory_requested else "skipped"
            )
            self._record_inference(
                ctx,
                inference_class=(
                    "replicate_aware_pseudobulk" if explicit_contrasts else "not_tested"
                ),
                inference_status="not_testable_insufficient_biological_replicates",
                claimable=False,
                biological_sample_col=sample_col,
            )
            if confirmatory_requested:
                raise ValueError(msg)
            ctx.status(self.name, "skipped", msg)
            return

        if "counts" not in adata.layers:
            msg = "Missing adata.layers['counts']; pseudobulk DE requires raw UMI counts."
            logger.warning("%s; skipping pseudobulk DE.", msg)
            ctx.metadata["pseudobulk_de_status"] = (
                "failed_missing_counts_layer"
                if confirmatory_requested
                else "skipped_missing_counts_layer"
            )
            ctx.metadata["pseudobulk_de_mode"] = (
                "failed" if confirmatory_requested else "skipped"
            )
            self._record_inference(
                ctx,
                inference_class=(
                    "replicate_aware_pseudobulk" if explicit_contrasts else "not_tested"
                ),
                inference_status="not_testable_missing_raw_counts",
                claimable=False,
                biological_sample_col=sample_col,
            )
            if confirmatory_requested:
                raise ValueError(msg)
            ctx.status(self.name, "skipped", msg)
            return

        try:
            if explicit_contrasts:
                replicate_contract = self._validate_biological_replicates(
                    adata,
                    sample_col=sample_col,
                    contrast_col=contrast_col,
                    contrasts=explicit_contrasts,
                    min_samples_per_condition=min_samples,
                )
                pb_counts, pb_meta = self._aggregate(
                    adata,
                    [sample_col, group_col, contrast_col],
                    min_cells_per_sample=min_cells,
                )
                results = self._run_confirmatory(
                    pb_counts,
                    pb_meta,
                    sample_col,
                    group_col,
                    contrast_col,
                    explicit_contrasts,
                    min_samples_per_condition=min_samples,
                )
                mode = "confirmatory"
                ctx.metadata["pseudobulk_de_contrast_contract"] = {
                    "mode": mode,
                    "sample_col": sample_col,
                    "group_col": group_col,
                    "contrast_col": contrast_col,
                    "contrasts": explicit_contrasts,
                    "biological_replicates": replicate_contract,
                }
            elif exploratory:
                pb_counts, pb_meta = self._aggregate(
                    adata,
                    [sample_col, group_col],
                    min_cells_per_sample=min_cells,
                )
                results = self._run_exploratory(
                    pb_counts,
                    pb_meta,
                    group_col,
                    min_samples_per_condition=min_samples,
                )
                mode = "exploratory"
                ctx.metadata["pseudobulk_de_contrast_contract"] = {
                    "mode": mode,
                    "sample_col": sample_col,
                    "group_col": group_col,
                    "contrast_policy": "group_vs_rest",
                }
            else:
                msg = contract_message or (
                    "No explicit pseudobulk contrast contract provided. Pass "
                    "--pseudobulk-contrast-col/--pseudobulk-contrast-a/--pseudobulk-contrast-b "
                    "for confirmatory pseudobulk or enable --pseudobulk-exploratory-group-vs-rest."
                )
                ctx.metadata["pseudobulk_de_status"] = "skipped_missing_contrast_contract"
                ctx.metadata["pseudobulk_de_mode"] = "skipped"
                self._record_inference(
                    ctx,
                    inference_class="not_tested",
                    inference_status="not_testable_missing_contrast_contract",
                    claimable=False,
                )
                ctx.status(self.name, "skipped", msg)
                return
        except PseudobulkInferenceContractError as exc:
            logger.warning("%s; skipping pseudobulk DE.", exc)
            ctx.metadata["pseudobulk_de_status"] = f"failed_{exc.inference_status}"
            ctx.metadata["pseudobulk_de_mode"] = "failed"
            self._record_inference(
                ctx,
                inference_class="replicate_aware_pseudobulk",
                inference_status=exc.inference_status,
                claimable=False,
                biological_sample_col=sample_col,
                contrast_col=contrast_col,
            )
            raise
        except ValueError as exc:
            logger.warning("%s; skipping pseudobulk DE.", exc)
            ctx.metadata["pseudobulk_de_status"] = (
                "failed_invalid_counts_or_grouping"
                if confirmatory_requested
                else "skipped_missing_counts_layer"
            )
            ctx.metadata["pseudobulk_de_mode"] = (
                "failed" if confirmatory_requested else "skipped"
            )
            self._record_inference(
                ctx,
                inference_class=(
                    "replicate_aware_pseudobulk"
                    if explicit_contrasts
                    else "exploratory_pseudobulk"
                ),
                inference_status="not_testable_invalid_counts_or_grouping",
                claimable=False,
            )
            if confirmatory_requested:
                raise
            ctx.status(self.name, "skipped", str(exc))
            return

        if pb_counts.empty or pb_meta.empty:
            ctx.metadata["pseudobulk_de_status"] = (
                "failed_no_valid_pseudobulk_groups"
                if confirmatory_requested
                else "skipped_no_valid_pseudobulk_groups"
            )
            ctx.metadata["pseudobulk_de_mode"] = (
                "failed" if confirmatory_requested else "skipped"
            )
            self._record_inference(
                ctx,
                inference_class=(
                    "replicate_aware_pseudobulk"
                    if explicit_contrasts
                    else "exploratory_pseudobulk"
                ),
                inference_status="not_testable_no_valid_pseudobulk_groups",
                claimable=False,
            )
            msg = "No pseudobulk groups passed minimum cell threshold."
            if confirmatory_requested:
                raise ValueError(msg)
            ctx.status(self.name, "skipped", msg)
            return

        if results is None or results.empty:
            ctx.metadata["pseudobulk_de_status"] = (
                "failed_no_testable_contrasts"
                if confirmatory_requested
                else "skipped_no_testable_contrasts"
            )
            ctx.metadata["pseudobulk_de_mode"] = (
                "failed" if confirmatory_requested else mode
            )
            self._record_inference(
                ctx,
                inference_class=(
                    "replicate_aware_pseudobulk"
                    if mode == "confirmatory"
                    else "exploratory_pseudobulk"
                ),
                inference_status="not_testable_insufficient_biological_replicates",
                claimable=False,
            )
            self._write_output_tables(pb_counts, pb_meta, ctx)
            msg = "No pseudobulk contrasts met minimum sample requirements."
            if confirmatory_requested:
                raise ValueError(msg)
            ctx.status(self.name, "skipped", msg)
            return

        test_used = set(results.get("test_used", pd.Series(dtype=str)).dropna().astype(str))
        if mode == "confirmatory" and test_used == {"pydeseq2"}:
            inference_class = "replicate_aware_pseudobulk"
            inference_status = "supported_confirmatory"
            claimable = True
        elif mode == "confirmatory":
            inference_class = "replicate_aware_pseudobulk_backend_fallback"
            inference_status = "exploratory_nonclaimable_backend_fallback"
            claimable = False
        else:
            inference_class = "exploratory_pseudobulk_group_vs_rest"
            inference_status = "exploratory_nonclaimable"
            claimable = False

        ctx.metadata["pseudobulk_de_contrast_contract"].update({
            "inference_class": inference_class,
            "inference_status": inference_status,
            "claimable": claimable,
        })
        self._record_inference(
            ctx,
            inference_class=inference_class,
            inference_status=inference_status,
            claimable=claimable,
            biological_sample_col=sample_col,
            contrast_col=contrast_col,
            test_used=sorted(test_used),
        )
        results["inference_class"] = inference_class
        results["inference_status"] = inference_status
        results["claimable"] = claimable
        results["biological_sample_col"] = sample_col
        results.to_csv(ctx.table_dir / "pseudobulk_de_results.csv", index=False)
        self._write_output_tables(pb_counts, pb_meta, ctx)
        n_sig = int((results["padj"] < 0.05).sum()) if "padj" in results.columns else 0
        ctx.metadata["pseudobulk_de_significant_genes"] = n_sig
        if mode == "confirmatory" and claimable:
            module_status = "completed_confirmatory"
        elif mode == "confirmatory":
            module_status = "completed_nonclaimable_backend_fallback"
        else:
            module_status = "completed_exploratory_nonclaimable"
        ctx.metadata["pseudobulk_de_status"] = module_status
        ctx.metadata["pseudobulk_de_mode"] = mode
        # H-3 audit fix (2026-05-22): manifest-level summary of which DE
        # test produced the rows. Squair 2021 (the cited pseudobulk paper)
        # presumes a DESeq2-class model; reviewers must be able to gate
        # publication claims by reading the manifest, not the per-row table.
        if "test_used" in results.columns:
            test_counts = results["test_used"].value_counts().to_dict()
            ctx.metadata["pseudobulk_de_test_used_counts"] = {
                str(k): int(v) for k, v in test_counts.items()
            }
            ctx.metadata["pseudobulk_de_test_actually_used"] = (
                "pydeseq2"
                if test_counts.get("pydeseq2", 0) > 0
                and test_counts.get("pydeseq2", 0) == len(results)
                else "mixed_or_rank_fallback"
            )
        self._plot_volcano(results, ctx)
        self._plot_heatmap(results, pb_counts, pb_meta, sample_col, ctx)

    # -- contract helpers ------------------------------------------------
    @staticmethod
    def _write_output_tables(
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        ctx: PipelineContext,
    ) -> None:
        counts.to_csv(ctx.table_dir / "pseudobulk_counts.csv", index=False)
        annotated_metadata = metadata.reset_index(drop=True).copy()
        annotated_metadata["inference_class"] = ctx.metadata[
            "pseudobulk_de_inference_class"
        ]
        annotated_metadata["inference_status"] = ctx.metadata[
            "pseudobulk_de_inference_status"
        ]
        annotated_metadata["claimable"] = ctx.metadata["pseudobulk_de_claimable"]
        annotated_metadata.to_csv(
            ctx.table_dir / "pseudobulk_metadata.csv",
            index=False,
        )

    @staticmethod
    def _validate_biological_replicates(
        adata,
        *,
        sample_col: str,
        contrast_col: str,
        contrasts: list[dict[str, str]],
        min_samples_per_condition: int,
    ) -> list[dict[str, object]]:
        identities = adata.obs[[sample_col, contrast_col]].copy()
        invalid_identity = identities.isna() | identities.apply(
            lambda column: column.astype(str).str.strip().eq("")
        )
        if invalid_identity.any(axis=None):
            raise PseudobulkInferenceContractError(
                "not_testable_invalid_biological_identifiers",
                "Biological sample and condition identities must be non-null and non-blank.",
            )
        mapping = identities.astype(str).drop_duplicates()
        conditions_per_sample = mapping.groupby(sample_col, observed=True)[
            contrast_col
        ].nunique()
        ambiguous = conditions_per_sample[conditions_per_sample != 1].index.astype(str).tolist()
        if ambiguous:
            raise PseudobulkInferenceContractError(
                "not_testable_sample_condition_nonunique",
                "Each biological sample must map to exactly one condition; ambiguous "
                f"samples: {ambiguous[:5]}",
            )

        evidence: list[dict[str, object]] = []
        for spec in contrasts:
            contrast_a = str(spec["contrast_a"])
            contrast_b = str(spec["contrast_b"])
            samples_a = set(
                mapping.loc[mapping[contrast_col].astype(str) == contrast_a, sample_col]
                .astype(str)
            )
            samples_b = set(
                mapping.loc[mapping[contrast_col].astype(str) == contrast_b, sample_col]
                .astype(str)
            )
            if not samples_a or not samples_b:
                raise PseudobulkInferenceContractError(
                    "not_testable_invalid_contrast_labels",
                    f"Contrast labels {contrast_a!r} and {contrast_b!r} must both "
                    f"exist in adata.obs[{contrast_col!r}].",
                )
            overlap = sorted(samples_a & samples_b)
            if overlap:
                raise PseudobulkInferenceContractError(
                    "not_testable_overlapping_biological_replicates",
                    f"Contrast {contrast_a} vs {contrast_b} reuses samples: {overlap[:5]}",
                )
            if (
                len(samples_a) < min_samples_per_condition
                or len(samples_b) < min_samples_per_condition
            ):
                raise PseudobulkInferenceContractError(
                    "not_testable_insufficient_biological_replicates",
                    f"Contrast {contrast_a} vs {contrast_b} has {len(samples_a)} vs "
                    f"{len(samples_b)} biological replicates; requires at least "
                    f"{min_samples_per_condition} per condition.",
                )
            evidence.append(
                {
                    "name": str(spec.get("name", f"{contrast_a}_vs_{contrast_b}")),
                    "contrast_a": contrast_a,
                    "contrast_b": contrast_b,
                    "n_biological_replicates_a": len(samples_a),
                    "n_biological_replicates_b": len(samples_b),
                }
            )
        return evidence

    @classmethod
    def _resolve_sample_col(cls, adata, preferred: str) -> str:
        return cls._resolve_col(
            adata,
            preferred,
            ("sample", "batch", "donor", "patient"),
            "sample/batch",
        )

    @classmethod
    def _resolve_group_col(cls, adata, preferred: str) -> str:
        if preferred:
            return cls._resolve_col(adata, preferred, ("cell_type", "leiden"), "grouping")
        return cls._resolve_col(adata, "cell_type", ("leiden",), "grouping")

    @classmethod
    def _resolve_optional_contrast_col(cls, adata, preferred: str | None) -> str | None:
        try:
            return cls._resolve_col(
                adata,
                preferred or "condition",
                ("condition", "disease", "group", "cohort"),
                "contrast",
            )
        except ValueError:
            return None

    def _resolve_confirmatory_contract(self, adata, ctx: PipelineContext):
        cfg = getattr(ctx.cfg, "pseudobulk", None)
        if cfg is None:
            return [], None, "No pseudobulk configuration found."

        contrast_specs: list[dict[str, str]] = []
        if getattr(cfg, "contrast_json", None):
            contrast_specs = self._load_contrast_json(Path(cfg.contrast_json))
        elif getattr(cfg, "contrast_a", None) and getattr(cfg, "contrast_b", None):
            contrast_specs = [
                {
                    "name": f"{cfg.contrast_a}_vs_{cfg.contrast_b}",
                    "contrast_a": str(cfg.contrast_a),
                    "contrast_b": str(cfg.contrast_b),
                }
            ]
        else:
            return [], None, "No explicit confirmatory pseudobulk contrast was provided."

        json_columns = {
            spec["contrast_col"] for spec in contrast_specs if spec.get("contrast_col")
        }
        configured_contrast_col = getattr(cfg, "contrast_col", None)
        requested_columns = set(json_columns)
        if configured_contrast_col:
            requested_columns.add(str(configured_contrast_col))
        if len(requested_columns) > 1:
            raise PseudobulkInferenceContractError(
                "not_testable_mixed_contrast_columns",
                "A pseudobulk contrast contract must use one consistent contrast_col; "
                f"received {sorted(requested_columns)}.",
            )
        requested_contrast_col = next(iter(requested_columns), None)
        if requested_contrast_col:
            contrast_col = (
                requested_contrast_col
                if requested_contrast_col in adata.obs.columns
                else None
            )
        else:
            contrast_col = self._resolve_optional_contrast_col(adata, None)
        if not contrast_col:
            return [], None, "Requested pseudobulk contrasts but no usable contrast column was found in adata.obs."

        for spec in contrast_specs:
            spec.setdefault("contrast_col", contrast_col)
        return contrast_specs, contrast_col, ""

    @staticmethod
    def _load_contrast_json(path: Path) -> list[dict[str, str]]:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, list):
            raise ValueError("pseudobulk contrast JSON must be a list of contrast specs.")
        if not payload:
            raise ValueError("pseudobulk contrast JSON must contain at least one contrast spec.")
        specs: list[dict[str, str]] = []
        for index, item in enumerate(payload):
            if not isinstance(item, dict):
                raise ValueError(f"pseudobulk contrast JSON item {index} must be an object.")
            a = item.get("contrast_a")
            b = item.get("contrast_b")
            if not a or not b:
                raise ValueError(
                    f"pseudobulk contrast JSON item {index} requires contrast_a and contrast_b."
                )
            spec = {
                "name": str(item.get("name") or f"{a}_vs_{b}"),
                "contrast_a": str(a),
                "contrast_b": str(b),
            }
            if item.get("contrast_col"):
                spec["contrast_col"] = str(item["contrast_col"])
            specs.append(spec)
        return specs

    # -- column helpers --------------------------------------------------
    @staticmethod
    def _resolve_col(adata, preferred: str, fallbacks: tuple, label: str) -> str:
        if preferred in adata.obs.columns:
            return preferred
        for c in fallbacks:
            if c in adata.obs.columns:
                logger.info("Using '%s' as %s column ('%s' not found).", c, label, preferred)
                return c
        raise ValueError(
            f"No {label} column found in adata.obs (tried '{preferred}' + {fallbacks})."
        )

    # -- aggregation -----------------------------------------------------
    @staticmethod
    def _aggregate(
        adata,
        group_cols: list[str] | str,
        group_col_legacy: str | None = None,
        min_cells_per_sample: int = MIN_CELLS_PER_SAMPLE,
    ):
        if isinstance(group_cols, str):
            cols = [group_cols]
            if group_col_legacy is not None:
                cols.append(group_col_legacy)
        else:
            cols = list(group_cols)
        raw = adata.layers["counts"]
        if hasattr(raw, "compute"):
            raw = raw.compute()
        elif hasattr(raw, "to_memory"):
            raw = raw.to_memory()
        if sparse.issparse(raw):
            raw = raw.tocsr()
        else:
            raw = sparse.csr_matrix(raw)
        if raw.shape != adata.shape:
            raise ValueError(
                f"counts layer shape {raw.shape} does not match adata shape {adata.shape}."
            )
        genes = adata.var_names.tolist()
        records, meta_rows = [], []
        for keys, idx in adata.obs.groupby(cols, observed=True).groups.items():
            if len(idx) < min_cells_per_sample:
                continue
            sub = raw[adata.obs.index.get_indexer(idx), :]
            sums = np.asarray(sub.sum(axis=0)).ravel()
            records.append(sums)
            if not isinstance(keys, tuple):
                keys = (keys,)
            row = {col: key for col, key in zip(cols, keys)}
            row["n_cells"] = len(idx)
            meta_rows.append(row)
        return pd.DataFrame(records, columns=genes), pd.DataFrame(meta_rows)

    def _run_confirmatory(
        self,
        pb_counts: pd.DataFrame,
        pb_meta: pd.DataFrame,
        sample_col: str,
        group_col: str,
        contrast_col: str,
        contrasts: list[dict[str, str]],
        min_samples_per_condition: int,
    ) -> pd.DataFrame | None:
        results: list[pd.DataFrame] = []
        for group in pb_meta[group_col].unique():
            mask = pb_meta[group_col] == group
            group_meta = pb_meta.loc[mask].reset_index(drop=True)
            group_meta[contrast_col] = group_meta[contrast_col].astype(str)
            group_counts = pb_counts.loc[mask].reset_index(drop=True)
            for spec in contrasts:
                a = spec["contrast_a"]
                b = spec["contrast_b"]
                samples_a = sorted(
                    group_meta.loc[group_meta[contrast_col] == a, sample_col]
                    .astype(str)
                    .unique()
                )
                samples_b = sorted(
                    group_meta.loc[group_meta[contrast_col] == b, sample_col]
                    .astype(str)
                    .unique()
                )
                res = self._run_de(
                    group_counts,
                    group_meta,
                    contrast_col,
                    a,
                    b,
                    min_samples_per_condition=min_samples_per_condition,
                )
                if res is None:
                    continue
                res["group"] = f"{group}|{spec['name']}"
                res["contrast_col"] = contrast_col
                res["contrast_a"] = a
                res["contrast_b"] = b
                res["n_biological_replicates_a"] = len(samples_a)
                res["n_biological_replicates_b"] = len(samples_b)
                res["biological_replicate_ids_a"] = ";".join(samples_a)
                res["biological_replicate_ids_b"] = ";".join(samples_b)
                results.append(res)
        if not results:
            return None
        return pd.concat(results, ignore_index=True)

    def _run_exploratory(
        self,
        pb_counts: pd.DataFrame,
        pb_meta: pd.DataFrame,
        group_col: str,
        min_samples_per_condition: int,
    ) -> pd.DataFrame | None:
        groups = pb_meta[group_col].unique()
        if len(groups) < 2:
            return None

        all_results: list[pd.DataFrame] = []
        if len(groups) == 2:
            res = self._run_de(
                pb_counts,
                pb_meta,
                group_col,
                groups[0],
                groups[1],
                min_samples_per_condition=min_samples_per_condition,
            )
            if res is not None:
                res["group"] = f"{groups[0]}_vs_{groups[1]}"
                all_results.append(res)
        else:
            for grp in groups:
                m = pb_meta.copy()
                m["condition"] = np.where(m[group_col] == grp, grp, "rest")
                res = self._run_de(
                    pb_counts,
                    m,
                    "condition",
                    grp,
                    "rest",
                    min_samples_per_condition=min_samples_per_condition,
                )
                if res is not None:
                    res["group"] = f"{grp}_vs_rest"
                    all_results.append(res)
        if not all_results:
            return None
        return pd.concat(all_results, ignore_index=True)

    # -- DE dispatch -----------------------------------------------------
    def _run_de(
        self,
        counts,
        meta,
        cond_col,
        ca,
        cb,
        min_samples_per_condition: int = MIN_SAMPLES_PER_CONDITION,
    ):
        na, nb = (meta[cond_col] == ca).sum(), (meta[cond_col] == cb).sum()
        if na < min_samples_per_condition or nb < min_samples_per_condition:
            logger.info("Skipping %s vs %s: samples %d vs %d.", ca, cb, na, nb)
            return None
        # H-3 audit fix (2026-05-22): record which test actually produced
        # the p-values. Squair 2021 (the cited pseudobulk paper) presumes a
        # DESeq2-class model; if we silently fall through to rank tests the
        # result rows must say so explicitly so reviewers cannot mis-cite.
        try:
            res = self._de_pydeseq2(counts, meta, cond_col, ca, cb)
            if res is not None:
                res["test_used"] = "pydeseq2"
            return res
        except ImportError:
            logger.info("pydeseq2 not available, falling back to Mann-Whitney U.")
        except Exception as exc:
            logger.warning("pydeseq2 failed (%s), falling back.", exc)
        try:
            from scipy.stats import mannwhitneyu

            res = self._de_ranktest(
                counts,
                meta,
                cond_col,
                ca,
                cb,
                mannwhitneyu,
                two_sided_kw={"alternative": "two-sided"},
            )
            if res is not None:
                res["test_used"] = "mannwhitneyu"
            return res
        except Exception as exc:
            logger.warning("Mann-Whitney failed (%s), falling back to Wilcoxon.", exc)
        from scipy.stats import ranksums

        res = self._de_ranktest(counts, meta, cond_col, ca, cb, ranksums)
        if res is not None:
            res["test_used"] = "ranksums"
        return res

    # -- backend: pydeseq2 ----------------------------------------------
    @staticmethod
    def _de_pydeseq2(counts, meta, cond_col, ca, cb):
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        ci = counts.round().astype(int)
        md = meta[[cond_col]].copy()
        md.index = ci.index = pd.RangeIndex(len(md))
        md[cond_col] = md[cond_col].astype(str)
        dds = DeseqDataSet(counts=ci, metadata=md, design=f"~{cond_col}")
        dds.deseq2()
        # Never silently drop the requested contrast: an older/incompatible
        # PyDESeq2 API must fall through to the explicitly non-claimable rank
        # backend rather than producing a default coefficient labeled claimable.
        sr = DeseqStats(dds, contrast=[cond_col, str(ca), str(cb)])
        sr.summary()
        r = sr.results_df.reset_index().rename(
            columns={"index": "gene", "log2FoldChange": "log2FC"}
        )
        return r[["gene", "log2FC", "pvalue", "padj"]].copy()

    # -- backend: generic rank test (MWU / ranksums) ---------------------
    @staticmethod
    def _de_ranktest(counts, meta, cond_col, ca, cb, test_fn, two_sided_kw=None):
        from .._mem_guard import MemoryGuard, MemoryAbortError
        kw = two_sided_kw or {}
        ia = meta[meta[cond_col] == ca].index
        ib = meta[meta[cond_col] == cb].index
        ma, mb = counts.loc[ia].values, counts.loc[ib].values

        # MEDIUM-4 fix (2026-05-19): the previous implementation computed
        # lfc as log2((mean(a)+1)/(mean(b)+1)) on raw pseudobulk counts,
        # which is biased against samples with deeper sequencing depth.
        # Normalize each pseudobulk to a per-sample CPM (counts per million)
        # before averaging so lfc is library-size-corrected. The rank test
        # itself is unaffected by per-sample monotone transforms, so the
        # p-value path stays equivalent.
        eps = 1.0
        def _cpm(matrix: np.ndarray) -> np.ndarray:
            totals = matrix.sum(axis=1, keepdims=True)
            totals = np.where(totals > 0, totals, 1.0)
            return matrix * 1e6 / totals

        ma_cpm = _cpm(ma.astype(np.float64))
        mb_cpm = _cpm(mb.astype(np.float64))

        genes, pvals, lfcs = [], [], []
        for j, g in enumerate(counts.columns):
            if MemoryGuard.abort_requested():
                raise MemoryAbortError("watchdog abort during pseudobulk DE gene loop")
            va, vb = ma[:, j], mb[:, j]
            va_cpm, vb_cpm = ma_cpm[:, j], mb_cpm[:, j]
            lfc = np.log2((va_cpm.mean() + eps) / (vb_cpm.mean() + eps))
            try:
                _, p = test_fn(va, vb, **kw)
            except ValueError:
                p = 1.0
            genes.append(g)
            pvals.append(p)
            lfcs.append(lfc)
        res = pd.DataFrame({"gene": genes, "log2FC": lfcs, "pvalue": pvals})
        res["padj"] = _bh_adjust(res["pvalue"].values)
        return res

    # -- plots -----------------------------------------------------------
    @staticmethod
    def _plot_volcano(results, ctx):
        if results.empty or "padj" not in results.columns:
            return
        fig, ax = plt.subplots(figsize=(10, 7))
        lp = -np.log10(results["padj"].clip(lower=1e-50))
        lfc = results["log2FC"]
        colors = np.where((results["padj"] < 0.05) & (lfc.abs() > 0.5), "tab:red", "grey")
        ax.scatter(lfc, lp, c=colors, s=6, alpha=0.6, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="grey", ls="--", lw=0.8)
        for v in (-0.5, 0.5):
            ax.axvline(v, color="grey", ls="--", lw=0.8)
        claim_label = PseudobulkDEModule._plot_claim_label(results)
        ax.set(
            xlabel="Log2 Fold Change",
            ylabel="-log10(adjusted p-value)",
            title=f"Pseudobulk DE -- Volcano Plot [{claim_label}]",
        )
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudobulk_volcano.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_heatmap(results, pb_counts, pb_meta, sample_col, ctx):
        if results.empty or "padj" not in results.columns:
            return
        top = (
            results[results["padj"] < 0.05]
            .sort_values("padj")
            .drop_duplicates("gene")
            .head(30)["gene"]
            .tolist()
        )
        avail = [g for g in top if g in pb_counts.columns]
        if len(avail) < 2:
            return
        mat = np.log2(pb_counts[avail] + 1)
        labels = pb_meta[sample_col].astype(str).values if sample_col in pb_meta.columns else np.arange(len(mat))
        fig, ax = plt.subplots(figsize=(max(8, len(avail) * 0.35), max(4, len(mat) * 0.4)))
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(avail)))
        ax.set_xticklabels(avail, rotation=90, fontsize=7)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        claim_label = PseudobulkDEModule._plot_claim_label(results)
        ax.set(
            ylabel="Sample",
            title=f"Top DE Genes -- Pseudobulk (log2 counts) [{claim_label}]",
        )
        plt.colorbar(im, ax=ax, label="log2(count + 1)")
        plt.tight_layout()
        plt.savefig(ctx.figure_dir / "pseudobulk_heatmap.png", dpi=160, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_claim_label(results: pd.DataFrame) -> str:
        claimable = bool(
            "claimable" in results.columns
            and not results.empty
            and results["claimable"].fillna(False).astype(bool).all()
        )
        if claimable:
            return "CLAIMABLE"
        if "inference_status" in results.columns and not results.empty:
            status = str(results["inference_status"].iloc[0])
        else:
            status = "unknown_inference_status"
        return f"NON-CLAIMABLE: {status}"


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return pvals.copy()
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    adj = pvals * n / ranks
    sorted_adj = adj[order]
    cummin = np.minimum.accumulate(sorted_adj[::-1])[::-1]
    adj[order] = cummin
    return np.clip(adj, 0.0, 1.0)
