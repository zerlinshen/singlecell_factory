from __future__ import annotations

import gzip
import hashlib
import json
import logging
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import jsonschema
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse

from ...factory_paths import R_MULTIOMICS_FACTORY_ROOT, SINGLECELL_FACTORY_ROOT, SUITE_ROOT
from ..context import PipelineContext
from ..manifest_writer import factory_git_state

logger = logging.getLogger(__name__)

__references__ = {
    "Squair_pseudobulk_2021": {
        "title": "Confronting false discoveries in single-cell differential expression",
        "authors": "Squair et al.",
        "journal": "Nature Communications",
        "year": "2021",
        "doi": "10.1038/s41467-021-25960-2",
        "description": "Sample-level pseudobulk avoids cell-level pseudoreplication.",
    },
    "Granja_ArchR_2021": {
        "title": "ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis",
        "authors": "Granja et al.",
        "journal": "Nature Genetics",
        "year": "2021",
        "doi": "10.1038/s41588-021-00790-6",
        "description": "Sample-aware pseudobulk replication for single-cell ATAC-seq accessibility modeling.",
    },
}

_SCHEMA_VERSION = "2.0.0"
_PEAK_AXIS_SCHEMA_VERSION = "atac_peak_axis/v1"
_INT64_MAX = int(np.iinfo(np.int64).max)
_R_INT_MAX = int(np.iinfo(np.int32).max)
_NA_LIKE = frozenset({"", "na", "nan", "none", "null", "<na>"})
_ENV_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


class ScatacContractError(ValueError):
    """Failure of the scATAC data, axis, design, or configuration contract."""


class ScatacInferenceError(RuntimeError):
    """An R-engine failure with enough state for a machine-readable failure record."""

    def __init__(
        self,
        message: str,
        *,
        completed_group_ids: Iterable[str] = (),
        failed_group_id: str | None = None,
    ) -> None:
        super().__init__(message)
        self.completed_group_ids = list(completed_group_ids)
        self.failed_group_id = failed_group_id


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(65536):
            h.update(chunk)
    return h.hexdigest()


def _nul_sha256(values: Iterable[str]) -> str:
    return hashlib.sha256(b"\0".join(str(value).encode("utf-8") for value in values)).hexdigest()


def _ordered_peak_id_sha256(peak_ids: Iterable[str]) -> str:
    return _nul_sha256(peak_ids)


def _frame_sha256(frame: pd.DataFrame) -> str:
    """Hash a deterministic TSV representation without changing biological values."""
    return hashlib.sha256(frame.to_csv(sep="\t", index=True, lineterminator="\n").encode("utf-8")).hexdigest()


def _sparse_matrix_sha256(matrix: sparse.spmatrix) -> str:
    """Hash a CSR matrix structurally; this never materializes a dense matrix."""
    csr = matrix.tocsr(copy=True)
    csr.sum_duplicates()
    csr.sort_indices()
    h = hashlib.sha256()
    h.update(f"{csr.shape[0]}x{csr.shape[1]}\0{csr.dtype}\0".encode("ascii"))
    for values in (csr.indptr, csr.indices, csr.data):
        contiguous = np.ascontiguousarray(values)
        h.update(str(contiguous.dtype).encode("ascii"))
        h.update(contiguous.tobytes())
    return h.hexdigest()


def _integer_vector_sha256(values: np.ndarray) -> str:
    return _nul_sha256(str(int(value)) for value in np.asarray(values).ravel())


def _exact_sparse_total(matrix: sparse.spmatrix) -> int:
    """Sum sparse stored nonnegative integers in Python-int space."""
    values = np.asarray(matrix.data)
    total = 0
    chunk_size = 65_536
    for start in range(0, values.size, chunk_size):
        total += sum(int(value) for value in values[start:start + chunk_size])
    return total


def _validate_sparse_stored_values(matrix: sparse.spmatrix) -> None:
    """Reject values that cannot participate in checked signed-int64 aggregation."""
    values = np.asarray(matrix.data)
    chunk_size = 65_536
    for start in range(0, values.size, chunk_size):
        chunk = values[start:start + chunk_size]
        if not np.all(np.isfinite(chunk)):
            raise ScatacContractError("atac_peaks contains non-finite (NaN/Inf) stored values.")
        if np.any(chunk < 0):
            raise ScatacContractError("atac_peaks contains negative stored counts.")
        if not np.all(chunk == np.floor(chunk)):
            raise ScatacContractError("atac_peaks contains fractional/non-integer stored counts.")
        if np.any(chunk > _INT64_MAX):
            raise ScatacContractError("atac_peaks contains a stored count outside signed int64 range.")


def _validated_labels(values: pd.Series, field_name: str) -> pd.Series:
    """Reject missing/blank/NA-like labels before any string conversion occurs."""
    null_mask = values.isna()
    if bool(null_mask.any()):
        raise ScatacContractError(
            f"{field_name} contains null values for {int(null_mask.sum())} analyzed cell(s)."
        )
    labels: list[str] = []
    for value in values.tolist():
        label = str(value).strip()
        if label.lower() in _NA_LIKE:
            raise ScatacContractError(
                f"{field_name} contains blank or NA-like values among analyzed cells."
            )
        labels.append(label)
    return pd.Series(labels, index=values.index, dtype="object")


def _validated_config_label(value: object, field_name: str) -> str:
    if value is None:
        raise ScatacContractError(f"{field_name} must be explicitly configured.")
    label = str(value).strip()
    if label.lower() in _NA_LIKE:
        raise ScatacContractError(f"{field_name} must be a non-blank, non-NA-like value.")
    return label


def _factory_states() -> dict[str, dict[str, Any]]:
    return {
        "suite": factory_git_state(SUITE_ROOT),
        "singlecell_factory": factory_git_state(SINGLECELL_FACTORY_ROOT),
        "r_multiomics_factory": factory_git_state(R_MULTIOMICS_FACTORY_ROOT),
    }


def _loaded_gpu_modules() -> list[str]:
    prefixes = ("cupy", "cudf", "rapids", "torch.cuda", "jax")
    return sorted(name for name in sys.modules if name.startswith(prefixes))


def _sanitize_argv(argv: list[str], out_dir: Path) -> list[str]:
    """Keep argv auditable while avoiding machine-specific absolute paths."""
    sanitized: list[str] = []
    out_root = out_dir.resolve()
    suite_root = SUITE_ROOT.resolve()
    for value in argv:
        path = Path(value)
        if not path.is_absolute():
            sanitized.append(value)
            continue
        candidate = path.resolve(strict=False)
        try:
            sanitized.append(str(candidate.relative_to(out_root)))
            continue
        except ValueError:
            pass
        try:
            sanitized.append(f"$SUITE_ROOT/{candidate.relative_to(suite_root)}")
            continue
        except ValueError:
            pass
        sanitized.append(f"$ABSOLUTE/{candidate.name}")
    return sanitized


def _schema_definition(definition: str) -> dict[str, Any]:
    schema_path = SINGLECELL_FACTORY_ROOT / "contracts" / "scatac_pseudobulk_da.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    # Definitions reference sibling definitions.  Preserve the full namespace
    # while selecting one document shape as the validation root.
    return {"definitions": schema["definitions"], **schema["definitions"][definition]}


def _write_validated_json(path: Path, payload: dict[str, Any], definition: str) -> Path:
    jsonschema.validate(payload, _schema_definition(definition))
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # Re-open after writing so a corrupt or partial file cannot be accepted only
    # because the in-memory object was valid.
    jsonschema.validate(json.loads(path.read_text(encoding="utf-8")), _schema_definition(definition))
    return path


def _resolve_r_adapter_path() -> Path:
    adapter = R_MULTIOMICS_FACTORY_ROOT / "bulk" / "rna_diff" / "scatac_pseudobulk_da.R"
    if not adapter.is_file():
        raise FileNotFoundError(
            "Canonical scATAC R adapter is missing at "
            f"{adapter}; resolve the r_multiomics_factory sibling through factory_paths."
        )
    return adapter


def _resolve_conda_executable() -> Path:
    """Resolve conda from the current runtime, never from CONDA_PREFIX."""
    candidates: list[Path] = []
    on_path = shutil.which("conda")
    if on_path:
        candidates.append(Path(on_path))
    executable = Path(sys.executable).resolve()
    # A conda Python conventionally lives at <conda>/envs/<env>/bin/python.
    if len(executable.parents) >= 4:
        candidates.append(executable.parents[3] / "bin" / "conda")
    for candidate in candidates:
        if candidate.is_file() and candidate.exists():
            return candidate
    raise FileNotFoundError(
        "conda executable is required to invoke the declared r_multiomics environment; "
        "no executable was found on PATH or relative to the active Python runtime."
    )


class ScatacPseudobulkDAModule:
    """CPU-only, sample-replicated scATAC pseudobulk DA adapter.

    The module is deliberately thin at the Python/R boundary: Python owns the
    AnnData-axis certificate, exact sparse aggregation, and provenance; the
    shared R helper owns DESeq2 primary inference and the mandatory edgeR QL
    cross-check.  It remains technically staged and never marks an arbitrary
    operator input as a biological claim.
    """

    name = "scatac_pseudobulk_da"
    fail_pipeline_on_error = True
    force_sequential = True
    parallel_safe = False
    mutates_structure = False
    requires_keys = {"obsm": ["atac_peaks"], "uns": ["atac_var", "atac_peak_axis"]}

    def run(self, ctx: PipelineContext) -> None:
        if ctx.adata is None:
            raise ScatacContractError("AnnData object required for scatac_pseudobulk_da.")

        out_dir = self._resolve_output_dir(ctx)
        cfg = getattr(ctx.cfg, "scatac_pseudobulk_da", None)
        if cfg is None:
            raise ScatacContractError("Missing scatac_pseudobulk_da configuration.")

        stage = "validation"
        try:
            val_info = self.validate_scatac_contract(ctx.adata, cfg)
            stage = "aggregation"
            agg_info = self.aggregate_sparse_counts(ctx.adata, cfg, val_info, out_dir)
            mode = getattr(cfg, "mode", "confirmatory_da")
            if mode == "aggregation_only":
                claim_path = self._write_claimability_manifest(
                    out_dir=out_dir,
                    cfg=cfg,
                    mode=mode,
                    technical_status="aggregation_completed",
                    group_records={},
                    agg_info=agg_info,
                )
                self._record_context_success(ctx, mode, claim_path, ())
                return

            stage = "inference"
            group_records = self.run_inference(cfg, val_info, agg_info, out_dir)
            claim_path = self._write_claimability_manifest(
                out_dir=out_dir,
                cfg=cfg,
                mode=mode,
                technical_status="dual_engine_inference_completed",
                group_records=group_records,
                agg_info=agg_info,
            )
            self._record_context_success(ctx, mode, claim_path, tuple(group_records))
        except Exception as exc:
            completed_group_ids = getattr(exc, "completed_group_ids", [])
            failed_group_id = getattr(exc, "failed_group_id", None)
            self._write_failure_manifest(
                out_dir,
                failed_stage=stage,
                error=exc,
                completed_group_ids=completed_group_ids,
                failed_group_id=failed_group_id,
            )
            ctx.metadata["scatac_pseudobulk_da_status"] = "failed"
            ctx.metadata["scatac_pseudobulk_da_error"] = str(exc)
            raise

    @staticmethod
    def _resolve_output_dir(ctx: PipelineContext) -> Path:
        lookup_module_dir = getattr(ctx, "module_output_dir", None)
        registered_dir = (
            lookup_module_dir(ScatacPseudobulkDAModule.name)
            if callable(lookup_module_dir)
            else None
        )
        if registered_dir is not None:
            out_dir = registered_dir
        elif getattr(ctx, "run_dir", None):
            out_dir = ctx.run_dir / ScatacPseudobulkDAModule.name
        elif hasattr(ctx, "output_dir"):
            out_dir = ctx.output_dir / ScatacPseudobulkDAModule.name
        else:
            out_dir = Path(getattr(ctx.cfg, "output_dir", ".")) / ScatacPseudobulkDAModule.name
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir

    @staticmethod
    def _record_context_success(
        ctx: PipelineContext,
        mode: str,
        claim_path: Path,
        group_ids: tuple[str, ...],
    ) -> None:
        ctx.metadata["scatac_pseudobulk_da_status"] = "completed"
        ctx.metadata["scatac_pseudobulk_da_mode"] = mode
        ctx.metadata["scatac_pseudobulk_da_claimable"] = False
        ctx.metadata["scatac_pseudobulk_da_claimability_manifest"] = str(claim_path)
        ctx.adata.uns["scatac_pseudobulk_da"] = {
            "claimable": False,
            "capability_release_status": "staged_nonproduction_missing_representative_real_data",
            "groups_evaluated": list(group_ids),
        }

    @staticmethod
    def _validate_config(cfg: Any) -> dict[str, Any]:
        mode = getattr(cfg, "mode", None)
        if mode not in {"confirmatory_da", "aggregation_only"}:
            raise ScatacContractError("mode must be exactly 'confirmatory_da' or 'aggregation_only'.")
        sample_col = _validated_config_label(getattr(cfg, "sample_col", None), "sample_col")
        group_col = _validated_config_label(getattr(cfg, "group_col", None), "group_col")
        peak_id_col = _validated_config_label(getattr(cfg, "peak_id_col", None), "peak_id_col")
        min_samples = int(getattr(cfg, "min_samples_per_condition", 0))
        if min_samples < 2:
            raise ScatacContractError(
                "min_samples_per_condition must be >= 2 biological samples per condition."
            )
        min_total_count = int(getattr(cfg, "min_total_count", -1))
        if min_total_count < 0:
            raise ScatacContractError("min_total_count must be a nonnegative integer.")
        fdr_threshold = float(getattr(cfg, "fdr_threshold", 0.0))
        if not 0.0 < fdr_threshold <= 1.0:
            raise ScatacContractError("fdr_threshold must be in the interval (0, 1].")
        abs_log2fc_threshold = float(getattr(cfg, "abs_log2fc_threshold", -1.0))
        if not np.isfinite(abs_log2fc_threshold) or abs_log2fc_threshold < 0.0:
            raise ScatacContractError("abs_log2fc_threshold must be finite and nonnegative.")
        r_conda_env = _validated_config_label(getattr(cfg, "r_conda_env", None), "r_conda_env")
        if not _ENV_NAME_RE.fullmatch(r_conda_env) or r_conda_env != "r_multiomics":
            raise ScatacContractError(
                "r_conda_env must be the declared 'r_multiomics' environment for this release."
            )
        timeout = int(getattr(cfg, "subprocess_timeout", 0))
        if timeout <= 0:
            raise ScatacContractError("subprocess_timeout must be a positive number of seconds.")
        aggregation_backend = _validated_config_label(
            getattr(cfg, "aggregation_backend", None), "aggregation_backend"
        )
        if aggregation_backend != "cpu":
            raise ScatacContractError(
                "Only the exact CPU sparse aggregation backend is supported in this staged release."
            )

        normalized: dict[str, Any] = {
            "mode": mode,
            "sample_col": sample_col,
            "group_col": group_col,
            "peak_id_col": peak_id_col,
            "min_samples": min_samples,
            "min_total_count": min_total_count,
            "fdr_threshold": fdr_threshold,
            "abs_log2fc_threshold": abs_log2fc_threshold,
            "r_conda_env": r_conda_env,
            "timeout": timeout,
            "aggregation_backend": aggregation_backend,
        }
        if mode == "confirmatory_da":
            test_level = _validated_config_label(getattr(cfg, "test_level", None), "test_level")
            reference_level = _validated_config_label(
                getattr(cfg, "reference_level", None), "reference_level"
            )
            if test_level == reference_level:
                raise ScatacContractError("test_level and reference_level must be distinct.")
            normalized.update(
                {
                    "condition_col": _validated_config_label(
                        getattr(cfg, "condition_col", None), "condition_col"
                    ),
                    "test_level": test_level,
                    "reference_level": reference_level,
                }
            )
        else:
            normalized.update({"condition_col": "", "test_level": "", "reference_level": ""})

        configured_groups = tuple(getattr(cfg, "groups", ()) or ())
        labels = tuple(_validated_config_label(value, "groups") for value in configured_groups)
        if len(labels) != len(set(labels)):
            raise ScatacContractError("groups must not contain duplicate labels.")
        normalized["groups"] = labels
        return normalized

    @classmethod
    def validate_scatac_contract(cls, adata, cfg) -> dict[str, Any]:
        config = cls._validate_config(cfg)
        if "atac_peaks" not in adata.obsm:
            raise ScatacContractError("Missing adata.obsm['atac_peaks'].")
        atac_peaks = adata.obsm["atac_peaks"]
        if isinstance(atac_peaks, np.ndarray) or not sparse.issparse(atac_peaks):
            raise ScatacContractError(
                "Dense matrix rejected: adata.obsm['atac_peaks'] must be a scipy sparse matrix."
            )
        if atac_peaks.shape[0] != adata.n_obs:
            raise ScatacContractError(
                f"atac_peaks row count ({atac_peaks.shape[0]}) does not match adata.n_obs ({adata.n_obs})."
            )
        _validate_sparse_stored_values(atac_peaks)

        if "atac_var" not in adata.uns:
            raise ScatacContractError("Missing adata.uns['atac_var'].")
        raw_atac_var = adata.uns["atac_var"]
        if isinstance(raw_atac_var, pd.DataFrame):
            atac_var = raw_atac_var.copy()
        elif isinstance(raw_atac_var, dict):
            atac_var = pd.DataFrame(raw_atac_var)
        else:
            raise ScatacContractError("adata.uns['atac_var'] must be a DataFrame or dict.")
        n_peaks = int(atac_peaks.shape[1])
        if len(atac_var) != n_peaks:
            raise ScatacContractError(
                f"atac_var row count ({len(atac_var)}) does not match atac_peaks column count ({n_peaks})."
            )
        peak_id_col = config["peak_id_col"]
        if peak_id_col not in atac_var.columns:
            raise ScatacContractError(
                f"Configured peak_id_col '{peak_id_col}' is absent from atac_var columns: "
                f"{list(atac_var.columns)}. Index fallback is prohibited."
            )
        peak_id_series = _validated_labels(atac_var[peak_id_col], f"atac_var['{peak_id_col}']")
        peak_ids = peak_id_series.tolist()
        if len(peak_ids) != len(set(peak_ids)):
            raise ScatacContractError("Peak IDs in atac_var must be strictly unique.")
        atac_var[peak_id_col] = peak_id_series
        peak_id_hash = _ordered_peak_id_sha256(peak_ids)

        certificate = adata.uns.get("atac_peak_axis")
        if not isinstance(certificate, dict):
            raise ScatacContractError(
                "Missing producer-emitted adata.uns['atac_peak_axis'] certificate; "
                "length and uniqueness cannot prove matrix-axis order."
            )
        required_certificate = {
            "schema_version", "producer", "source_status", "peak_id_col", "n_peaks", "ordered_peak_id_sha256"
        }
        missing_certificate = sorted(required_certificate - set(certificate))
        if missing_certificate:
            raise ScatacContractError(
                f"atac_peak_axis certificate is incomplete: missing {missing_certificate}."
            )
        if certificate["schema_version"] != _PEAK_AXIS_SCHEMA_VERSION:
            raise ScatacContractError("atac_peak_axis certificate has an unsupported schema version.")
        if certificate["producer"] != "atac_ingest":
            raise ScatacContractError("atac_peak_axis certificate was not emitted by atac_ingest.")
        if certificate["peak_id_col"] != peak_id_col:
            raise ScatacContractError(
                "Configured peak_id_col does not match the producer-emitted atac_peak_axis certificate."
            )
        if int(certificate["n_peaks"]) != n_peaks:
            raise ScatacContractError("atac_peak_axis certificate peak count does not match atac_peaks.")
        if certificate["ordered_peak_id_sha256"] != peak_id_hash:
            raise ScatacContractError(
                "Ordered peak-ID digest does not match the producer-emitted atac_peak_axis certificate."
            )

        obs = adata.obs
        for field in (config["sample_col"], config["group_col"]):
            if field not in obs.columns:
                raise ScatacContractError(f"Configured obs column '{field}' is absent.")
        # Group labels are validated before allowlisting so unseen missing values
        # cannot become implicit groups after later string conversion.
        all_group_labels = _validated_labels(obs[config["group_col"]], config["group_col"])
        configured_groups = config["groups"]
        if configured_groups:
            available_groups = set(all_group_labels.tolist())
            unknown_groups = sorted(set(configured_groups) - available_groups)
            if unknown_groups:
                raise ScatacContractError(
                    f"Configured groups are absent from obs['{config['group_col']}']: {unknown_groups}."
                )
            group_mask = all_group_labels.isin(configured_groups)
        else:
            group_mask = pd.Series(True, index=obs.index)
        if not bool(group_mask.any()):
            raise ScatacContractError("No cells remain after applying the configured group allowlist.")

        selected_obs = obs.loc[group_mask].copy()
        sample_labels = _validated_labels(selected_obs[config["sample_col"]], config["sample_col"])
        group_labels = _validated_labels(selected_obs[config["group_col"]], config["group_col"])
        if config["mode"] == "confirmatory_da":
            condition_col = config["condition_col"]
            if condition_col not in selected_obs.columns:
                raise ScatacContractError(f"Configured obs column '{condition_col}' is absent.")
            condition_labels = _validated_labels(selected_obs[condition_col], condition_col)
            observed_levels = set(condition_labels.tolist())
            configured_levels = {config["reference_level"], config["test_level"]}
            if observed_levels != configured_levels:
                raise ScatacContractError(
                    "Confirmatory DA requires exactly the configured two condition levels after group allowlisting; "
                    f"observed={sorted(observed_levels)}, configured={sorted(configured_levels)}."
                )
            design = pd.DataFrame(
                {"sample": sample_labels, "group": group_labels, "condition": condition_labels},
                index=selected_obs.index,
            )
            sample_condition_counts = design.groupby("sample", sort=False)["condition"].nunique(dropna=False)
            inconsistent = sample_condition_counts[sample_condition_counts != 1]
            if not inconsistent.empty:
                raise ScatacContractError(
                    "Sample-to-condition integrity error: biological sample(s) mapped to multiple conditions: "
                    f"{inconsistent.index.tolist()}."
                )
            groups = list(configured_groups) if configured_groups else sorted(group_labels.unique().tolist())
            underpowered: dict[str, dict[str, int]] = {}
            for group in groups:
                group_design = design.loc[design["group"] == group]
                counts = group_design.groupby("condition")["sample"].nunique()
                n_reference = int(counts.get(config["reference_level"], 0))
                n_test = int(counts.get(config["test_level"], 0))
                if n_reference < config["min_samples"] or n_test < config["min_samples"]:
                    underpowered[group] = {
                        config["reference_level"]: n_reference,
                        config["test_level"]: n_test,
                    }
            if underpowered:
                raise ScatacContractError(
                    f"Confirmatory DA aborts: group(s) have fewer than {config['min_samples']} "
                    f"biological samples per condition: {underpowered}. Groups are never silently skipped."
                )
        else:
            design = pd.DataFrame(
                {"sample": sample_labels, "group": group_labels, "condition": ""},
                index=selected_obs.index,
            )
            groups = list(configured_groups) if configured_groups else sorted(group_labels.unique().tolist())

        selected_indices = np.flatnonzero(group_mask.to_numpy())
        return {
            "atac_peaks": atac_peaks,
            "atac_var": atac_var,
            "peak_ids": peak_ids,
            "peak_id_hash": peak_id_hash,
            "peak_axis_certificate": dict(certificate),
            "config": config,
            "groups": groups,
            "selected_indices": selected_indices,
            "selected_design": design,
            "selected_design_sha256": _frame_sha256(design),
        }

    @staticmethod
    def aggregate_sparse_counts(adata, cfg, val_info: dict[str, Any], out_dir: Path) -> dict[str, Any]:
        del adata, cfg  # Contract-normalized values in val_info are the only accepted input.
        agg_dir = out_dir / "aggregation"
        agg_dir.mkdir(parents=True, exist_ok=True)
        selected_indices = val_info["selected_indices"]
        design = val_info["selected_design"].copy()
        peak_ids = val_info["peak_ids"]
        atac_peaks = val_info["atac_peaks"]
        config = val_info["config"]

        pb_units = (
            design.assign(_row_position=np.arange(len(design), dtype=np.int64))
            .groupby(["group", "sample", "condition"], observed=True, sort=True)
            .size()
            .reset_index(name="cell_count")
            .sort_values(["group", "sample", "condition"], kind="mergesort")
            .reset_index(drop=True)
        )
        pb_units["pb_id"] = [f"pb_{index:04d}" for index in range(len(pb_units))]
        unit_to_idx = {
            (row.group, row.sample, row.condition): index
            for index, row in enumerate(pb_units.itertuples(index=False))
        }
        cell_pb_idx = np.fromiter(
            (unit_to_idx[(row.group, row.sample, row.condition)] for row in design.itertuples(index=False)),
            dtype=np.int64,
            count=len(design),
        )
        n_cells = int(len(selected_indices))
        n_pb = int(len(pb_units))
        n_peaks = int(len(peak_ids))
        if n_cells == 0 or n_pb == 0:
            raise ScatacContractError("Aggregation requires at least one selected cell and pseudobulk sample.")

        x_sub = atac_peaks[selected_indices].tocsr(copy=True)
        x_sub.sum_duplicates()
        x_sub.sort_indices()
        selected_total_exact = _exact_sparse_total(x_sub)
        if selected_total_exact > _INT64_MAX:
            raise ScatacContractError(
                "Checked int64 aggregation overflow: exact selected input total "
                f"{selected_total_exact} exceeds signed int64 maximum {_INT64_MAX}."
            )
        x_int = x_sub.astype(np.int64, copy=False)
        membership = sparse.csr_matrix(
            (
                np.ones(n_cells, dtype=np.int64),
                (np.arange(n_cells, dtype=np.int64), cell_pb_idx),
            ),
            shape=(n_cells, n_pb),
            dtype=np.int64,
        )
        aggregate = (x_int.T.tocsr() @ membership).tocsr()
        aggregate.sum_duplicates()
        aggregate.sort_indices()

        aggregate_total_exact = _exact_sparse_total(aggregate)
        if selected_total_exact != aggregate_total_exact:
            raise ScatacContractError(
                "Count invariant violation: exact input total "
                f"{selected_total_exact} != aggregate total {aggregate_total_exact}."
            )
        input_per_peak = np.asarray(x_int.sum(axis=0)).ravel().astype(np.int64, copy=False)
        output_per_peak = np.asarray(aggregate.sum(axis=1)).ravel().astype(np.int64, copy=False)
        if not np.array_equal(input_per_peak, output_per_peak):
            raise ScatacContractError("Per-peak count invariant violation after sparse aggregation.")
        cell_totals = np.asarray(x_int.sum(axis=1)).ravel().astype(np.int64, copy=False)
        input_per_pseudobulk = np.zeros(n_pb, dtype=np.int64)
        for pseudobulk_index, cell_total in zip(cell_pb_idx, cell_totals, strict=True):
            input_per_pseudobulk[int(pseudobulk_index)] += int(cell_total)
        output_per_pseudobulk = np.asarray(aggregate.sum(axis=0)).ravel().astype(np.int64, copy=False)
        if not np.array_equal(input_per_pseudobulk, output_per_pseudobulk):
            raise ScatacContractError("Per-pseudobulk-sample count invariant violation after aggregation.")

        peaks_file = agg_dir / "peaks.tsv.gz"
        peaks_df = val_info["atac_var"].copy()
        peaks_df.insert(0, "peak_index", np.arange(n_peaks, dtype=np.int64))
        with gzip.open(peaks_file, "wt", encoding="utf-8") as handle:
            peaks_df.to_csv(handle, sep="\t", index=False)
        samples_file = agg_dir / "pseudobulk_samples.tsv"
        pb_samples_df = pd.DataFrame(
            {
                "pb_index": np.arange(n_pb, dtype=np.int64),
                "pb_id": pb_units["pb_id"],
                "group": pb_units["group"],
                "sample": pb_units["sample"],
                "condition": pb_units["condition"],
                "cell_count": pb_units["cell_count"].astype(np.int64),
            }
        )
        pb_samples_df.to_csv(samples_file, sep="\t", index=False)
        counts_mtx_file = agg_dir / "peak_by_pseudobulk_sample_counts.mtx.gz"
        with gzip.open(counts_mtx_file, "wb") as handle:
            spio.mmwrite(handle, aggregate)

        backend_observation = {
            "requested": config["aggregation_backend"],
            "resolved": "cpu",
            "matrix_implementation": f"{type(x_int).__module__}.{type(x_int).__name__}",
            "gpu_modules_loaded": _loaded_gpu_modules(),
            "cpu_only_observed": not bool(_loaded_gpu_modules()),
        }
        agg_manifest = {
            "schema_version": _SCHEMA_VERSION,
            "module_name": "scatac_pseudobulk_da",
            "timestamp_utc": _utc_now(),
            "n_cells": n_cells,
            "n_peaks": n_peaks,
            "n_pseudobulk_samples": n_pb,
            "peak_id_hash": val_info["peak_id_hash"],
            "peak_axis_certificate": val_info["peak_axis_certificate"],
            "aggregation_backend_requested": config["aggregation_backend"],
            "aggregation_backend_resolved": "cpu",
            "backend_observation": backend_observation,
            "integer_parity_certificate": None,
            "production_dual_run": False,
            "sample_col": config["sample_col"],
            "group_col": config["group_col"],
            "condition_col": config["condition_col"],
            "peak_id_col": config["peak_id_col"],
            "count_invariants": {
                "input_total_exact": str(selected_total_exact),
                "output_total_exact": str(aggregate_total_exact),
                "per_peak_sha256": _integer_vector_sha256(output_per_peak),
                "per_pseudobulk_sample_sha256": _integer_vector_sha256(output_per_pseudobulk),
            },
            "input_hashes": {
                "selected_sparse_counts_sha256": _sparse_matrix_sha256(x_int),
                "selected_design_sha256": val_info["selected_design_sha256"],
            },
            "repository_state": _factory_states(),
            "files": {
                "peaks_file": str(peaks_file.relative_to(out_dir)),
                "pseudobulk_samples_file": str(samples_file.relative_to(out_dir)),
                "counts_matrix_file": str(counts_mtx_file.relative_to(out_dir)),
            },
            "artifact_hashes": {
                "peaks_file": _compute_sha256(peaks_file),
                "pseudobulk_samples_file": _compute_sha256(samples_file),
                "counts_matrix_file": _compute_sha256(counts_mtx_file),
            },
        }
        agg_manifest_path = _write_validated_json(
            agg_dir / "aggregation_manifest.json", agg_manifest, "aggregation_manifest"
        )
        return {
            "aggregate": aggregate,
            "pb_samples_df": pb_samples_df,
            "peak_ids": peak_ids,
            "agg_dir": agg_dir,
            "agg_manifest": agg_manifest,
            "agg_manifest_path": agg_manifest_path,
            "peaks_file": peaks_file,
            "samples_file": samples_file,
            "counts_mtx_file": counts_mtx_file,
        }

    def run_inference(
        self,
        cfg: Any,
        val_info: dict[str, Any],
        agg_info: dict[str, Any],
        out_dir: Path,
    ) -> dict[str, dict[str, Any]]:
        del cfg
        groups_dir = out_dir / "groups"
        groups_dir.mkdir(parents=True, exist_ok=True)
        config = val_info["config"]
        r_adapter = _resolve_r_adapter_path()
        conda = _resolve_conda_executable()
        aggregate = agg_info["aggregate"]
        pb_samples_df = agg_info["pb_samples_df"]
        peak_ids = agg_info["peak_ids"]
        group_records: dict[str, dict[str, Any]] = {}

        for group_label in val_info["groups"]:
            group_id = "grp_" + hashlib.sha256(group_label.encode("utf-8")).hexdigest()[:12]
            group_dir = groups_dir / group_id
            group_dir.mkdir(parents=True, exist_ok=True)
            try:
                group_samples = pb_samples_df.loc[pb_samples_df["group"] == group_label].copy()
                if group_samples.empty:
                    raise ScatacContractError(f"Configured group '{group_label}' has no pseudobulk samples.")
                group_indices = group_samples["pb_index"].to_numpy(dtype=np.int64)
                group_counts = aggregate[:, group_indices].tocsr()
                if group_counts.data.size and int(np.max(group_counts.data)) > _R_INT_MAX:
                    raise ScatacContractError(
                        "R integer-range guard: a pseudobulk peak count exceeds "
                        f".Machine$integer.max ({_R_INT_MAX}) before R coercion."
                    )
                counts_tsv = group_dir / "counts.tsv"
                # Materialization is allowed only at the narrow R boundary; the
                # cell-by-peak source matrix remains sparse throughout aggregation.
                counts_df = pd.DataFrame(
                    group_counts.toarray(), index=peak_ids, columns=group_samples["sample"].tolist()
                )
                counts_df.index.name = "peak_id"
                counts_df.to_csv(counts_tsv, sep="\t")
                design_tsv = group_dir / "design.tsv"
                group_samples[["sample", "condition", "group", "cell_count"]].to_csv(
                    design_tsv, sep="\t", index=False
                )
                cmd = [
                    str(conda), "run", "-n", config["r_conda_env"], "Rscript", str(r_adapter),
                    "--counts", str(counts_tsv),
                    "--design", str(design_tsv),
                    "--test", config["test_level"],
                    "--ref", config["reference_level"],
                    "--out-dir", str(group_dir),
                    "--group-id", group_id,
                    "--group-label", group_label,
                    "--padj", str(config["fdr_threshold"]),
                    "--lfc", str(config["abs_log2fc_threshold"]),
                    "--min-total-count", str(config["min_total_count"]),
                ]
                try:
                    proc = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=config["timeout"],
                        check=False,
                    )
                except subprocess.TimeoutExpired as exc:
                    stdout = exc.stdout if isinstance(exc.stdout, str) else ""
                    stderr = exc.stderr if isinstance(exc.stderr, str) else ""
                    (group_dir / "r_stdout.log").write_text(stdout, encoding="utf-8")
                    (group_dir / "r_stderr.log").write_text(stderr, encoding="utf-8")
                    raise ScatacInferenceError(
                        f"scATAC DA inference timed out after {config['timeout']} seconds for group '{group_label}'.",
                        completed_group_ids=group_records,
                        failed_group_id=group_id,
                    ) from exc
                stdout_path = group_dir / "r_stdout.log"
                stderr_path = group_dir / "r_stderr.log"
                stdout_path.write_text(proc.stdout, encoding="utf-8")
                stderr_path.write_text(proc.stderr, encoding="utf-8")
                if proc.returncode != 0:
                    raise ScatacInferenceError(
                        f"scATAC DA inference failed for group '{group_label}' (exit code {proc.returncode}): "
                        f"{proc.stderr[:1000]}",
                        completed_group_ids=group_records,
                        failed_group_id=group_id,
                    )
                result_path = group_dir / "da_results.tsv.gz"
                r_manifest_path = group_dir / "inference_manifest.json"
                if not result_path.is_file() or not r_manifest_path.is_file():
                    raise ScatacInferenceError(
                        f"R adapter did not emit required result and manifest artifacts for group '{group_label}'.",
                        completed_group_ids=group_records,
                        failed_group_id=group_id,
                    )
                with gzip.open(result_path, "rt", encoding="utf-8") as handle:
                    result = pd.read_csv(handle, sep="\t")
                if len(result) != len(peak_ids) or set(result.get("peak_id", [])) != set(peak_ids):
                    raise ScatacInferenceError(
                        f"R adapter did not retain every input peak for group '{group_label}'.",
                        completed_group_ids=group_records,
                        failed_group_id=group_id,
                    )
                manifest = json.loads(r_manifest_path.read_text(encoding="utf-8"))
                # R's raw commandArgs can expose host-specific launch paths.
                # The governed record below is the only argv representation we
                # retain: sanitized arguments plus a hash of the actual vector.
                manifest.pop("command_argv", None)
                manifest.update(
                    {
                        "schema_version": _SCHEMA_VERSION,
                        "technical_inference_status": "dual_engine_inference_completed",
                        "adapter": {
                            "path": "$SUITE_ROOT/r_multiomics_factory/bulk/rna_diff/scatac_pseudobulk_da.R",
                            "sha256": _compute_sha256(r_adapter),
                            "version": manifest.get("adapter_version", "unknown"),
                        },
                        "argv": {
                            "sanitized_relative": _sanitize_argv(cmd, out_dir),
                            "nul_delimited_sha256": _nul_sha256(cmd),
                        },
                        "input_hashes": {
                            "counts_tsv_sha256": _compute_sha256(counts_tsv),
                            "design_tsv_sha256": _compute_sha256(design_tsv),
                            "peak_id_sha256": val_info["peak_id_hash"],
                        },
                        "artifact_hashes": {
                            "da_results_tsv_gz": _compute_sha256(result_path),
                            "counts_tsv": _compute_sha256(counts_tsv),
                            "design_tsv": _compute_sha256(design_tsv),
                            "r_stdout_log": _compute_sha256(stdout_path),
                            "r_stderr_log": _compute_sha256(stderr_path),
                        },
                        "repository_state": _factory_states(),
                        "cpu_backend": {
                            "aggregation_backend_resolved": "cpu",
                            "r_conda_env": config["r_conda_env"],
                            "production_dual_run": False,
                            "integer_parity_certificate": None,
                        },
                    }
                )
                _write_validated_json(r_manifest_path, manifest, "inference_manifest")
                group_records[group_id] = {
                    "group_label": group_label,
                    "manifest": manifest,
                    "directory": group_dir,
                    "result_path": result_path,
                    "manifest_path": r_manifest_path,
                    "stdout_path": stdout_path,
                    "stderr_path": stderr_path,
                }
            except ScatacInferenceError:
                raise
            except Exception as exc:
                raise ScatacInferenceError(
                    str(exc), completed_group_ids=group_records, failed_group_id=group_id
                ) from exc
        return group_records

    @staticmethod
    def _write_claimability_manifest(
        *,
        out_dir: Path,
        cfg: Any,
        mode: str,
        technical_status: str,
        group_records: dict[str, dict[str, Any]],
        agg_info: dict[str, Any],
    ) -> Path:
        per_group = {}
        group_artifacts: dict[str, dict[str, str]] = {}
        for group_id, record in group_records.items():
            manifest = record["manifest"]
            per_group[group_id] = {
                "biological_sample_counts_by_condition": manifest["biological_sample_counts_by_condition"],
                "cell_counts_by_condition": manifest["cell_counts_by_condition"],
                "n_deseq2_significant": manifest["n_deseq2_significant"],
                "n_edger_significant": manifest["n_edger_significant"],
                "n_cross_engine_supported": manifest["n_cross_engine_supported"],
                "n_concordance_denominator": manifest["n_concordance_denominator"],
                "concordance_fraction_edger": manifest["concordance_fraction_edger"],
            }
            group_artifacts[group_id] = {
                "result": _compute_sha256(record["result_path"]),
                "manifest": _compute_sha256(record["manifest_path"]),
                "stdout": _compute_sha256(record["stdout_path"]),
                "stderr": _compute_sha256(record["stderr_path"]),
            }
        config = agg_info["agg_manifest"]
        payload = {
            "schema_version": _SCHEMA_VERSION,
            "module_name": "scatac_pseudobulk_da",
            "execution_mode": mode,
            "technical_inference_status": technical_status,
            "scientific_release_status": "staged_nonproduction_missing_representative_real_data",
            "capability_release_status": "staged_nonproduction_missing_representative_real_data",
            "claimable": False,
            "biological_unit": "biological_sample",
            "groups_evaluated": list(group_records),
            "completed_group_ids": list(group_records),
            "failed_group_ids": [],
            "aggregation_backend": {
                "aggregation_backend_requested": config["aggregation_backend_requested"],
                "aggregation_backend_resolved": config["aggregation_backend_resolved"],
                "integer_parity_certificate": None,
                "production_dual_run": False,
            },
            "thresholds": {
                "fdr_threshold": float(config.get("fdr_threshold", getattr(cfg, "fdr_threshold", 0.05))),
                "abs_log2fc_threshold": float(getattr(cfg, "abs_log2fc_threshold", 1.0)),
                "min_total_count": int(getattr(cfg, "min_total_count", 10)),
                "min_samples_per_condition": int(getattr(cfg, "min_samples_per_condition", 2)),
            },
            "blockers": ["missing_representative_real_data_with_biological_samples_and_two_level_contrast"],
            "per_group_sample_denominators": per_group,
            "repository_state": _factory_states(),
            "software_versions": {
                "python": sys.version.split()[0],
                "r_multiomics": (
                    next(iter(group_records.values()))["manifest"]["software_versions"]
                    if group_records else {}
                ),
            },
            "artifact_hashes": {
                "aggregation_manifest": _compute_sha256(agg_info["agg_manifest_path"]),
                "peaks": _compute_sha256(agg_info["peaks_file"]),
                "pseudobulk_samples": _compute_sha256(agg_info["samples_file"]),
                "counts_matrix": _compute_sha256(agg_info["counts_mtx_file"]),
                "groups": group_artifacts,
            },
        }
        return _write_validated_json(out_dir / "claimability_manifest.json", payload, "claimability_manifest")

    @staticmethod
    def _write_failure_manifest(
        out_dir: Path,
        *,
        failed_stage: str,
        error: Exception,
        completed_group_ids: Iterable[str],
        failed_group_id: str | None,
    ) -> Path:
        partial_hashes: dict[str, str] = {}
        for candidate in sorted(out_dir.rglob("*")) if out_dir.exists() else []:
            if candidate.is_file() and candidate.name != "failure_manifest.json":
                try:
                    partial_hashes[str(candidate.relative_to(out_dir))] = _compute_sha256(candidate)
                except OSError:
                    continue
        payload = {
            "schema_version": _SCHEMA_VERSION,
            "module_name": "scatac_pseudobulk_da",
            "status": f"{failed_stage}_failed",
            "failed_stage": failed_stage,
            "error_type": type(error).__name__,
            "error_message": str(error),
            "timestamp_utc": _utc_now(),
            "completed_group_ids": list(completed_group_ids),
            "failed_group_id": failed_group_id,
            "capability_release_status": "staged_nonproduction_missing_representative_real_data",
            "claimable": False,
            "partial_artifact_hashes": partial_hashes,
        }
        return _write_validated_json(out_dir / "failure_manifest.json", payload, "failure_manifest")
