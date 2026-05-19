from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import anndata as ad
from anndata.experimental import concat_on_disk
import numpy as np
import pandas as pd
import zarr
from scipy import io

ROOT = Path("/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory/data/raw/nc2024_nsclc_emtab13526")
SDRF = ROOT / "E-MTAB-13526.sdrf.txt"
FULL_ROOT = ROOT / "full_cohort"
PARTS_ROOT = FULL_ROOT / "_prepared_parts_zarr"
PREPARED_ZARR = FULL_ROOT / "prepared_input.zarr"
SUMMARY_JSON = FULL_ROOT / "prepared_input.summary.json"
READY_SENTINEL = FULL_ROOT / "prepared_input.ready"
EXPECTED_N_SAMPLES = 81
CELL_CALLING_VERSION = "nc2024_hybrid_knee_qc_v1"
MIN_UMI_FLOOR = 1_000
MAX_UMI_CAP = 2_000
MIN_GENES_FLOOR = 500
MAX_TOP_GENE_FRACTION = 0.60
MAX_MITO_FRACTION = 0.35
MIN_RETAINED_BARCODES_PER_SAMPLE = 100
MAX_VALIDATED_RETENTION_FRACTION = 0.05

REQUIRED_OBS_COLUMNS = [
    "sample",
    "patient",
    "batch",
    "disease",
    "condition",
    "sorting",
    "sampling_site",
    "sex",
    "original_source_name",
    "tumor_type",
]


def read_features(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None)
    if df.shape[1] == 3:
        df.columns = ["gene_id", "gene_name", "feature_type"]
    elif df.shape[1] == 2:
        df.columns = ["gene_id", "gene_name"]
        df["feature_type"] = "Gene Expression"
    else:
        raise ValueError(f"Unexpected feature file format: {path}")
    return df


def read_barcodes(path: Path) -> pd.Series:
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", header=None)[0]


def read_matrix(path: Path):
    with gzip.open(path, "rb") as fh:
        return io.mmread(fh).tocsr().T.tocsr()


def build_sample_meta(sdrf_path: Path) -> pd.DataFrame:
    sdrf = pd.read_csv(sdrf_path, sep="\t")
    keep = [
        "Source Name",
        "Characteristics[individual]",
        "Characteristics[disease]",
        "Characteristics[FACS]",
        "Characteristics[sampling site]",
        "Characteristics[sex]",
        "Characteristics[original source name]",
    ]
    meta = sdrf[keep].drop_duplicates().rename(
        columns={
            "Source Name": "sample",
            "Characteristics[individual]": "patient",
            "Characteristics[disease]": "disease",
            "Characteristics[FACS]": "sorting",
            "Characteristics[sampling site]": "sampling_site",
            "Characteristics[sex]": "sex",
            "Characteristics[original source name]": "original_source_name",
        }
    )
    meta = meta.astype(str)
    meta["condition"] = np.where(
        meta["disease"].str.lower().eq("normal"),
        "healthy_background",
        "tumor",
    )
    meta["tumor_type"] = np.where(
        meta["condition"].eq("tumor"),
        "NSCLC",
        "non_involved",
    )
    meta["batch"] = meta["sample"]
    return meta.sort_values("sample").reset_index(drop=True)


def make_var_out(var_ref: pd.DataFrame) -> pd.DataFrame:
    var_out = var_ref.copy()
    var_out.index = pd.Index(var_out["gene_name"].astype(str), name="gene_symbol")
    var_out["gene_symbol"] = var_out["gene_name"].astype(str)
    var_out = var_out.drop(columns=["gene_name"])
    return var_out


def safe_fraction(numerator: float, denominator: float) -> float:
    if denominator <= 0:
        return 0.0
    return float(numerator) / float(denominator)


def compute_row_max(mat) -> np.ndarray:
    row_nnz = np.diff(mat.indptr)
    row_max = np.zeros(mat.shape[0], dtype=np.float32)
    nonempty = row_nnz > 0
    if np.any(nonempty):
        starts = mat.indptr[:-1][nonempty]
        row_max[nonempty] = np.maximum.reduceat(mat.data, starts)
    return row_max


def estimate_knee_umi_threshold(total_counts: np.ndarray) -> int:
    positive = np.asarray(total_counts[total_counts > 0], dtype=np.int64)
    if positive.size == 0:
        return MIN_UMI_FLOOR
    positive = np.sort(positive)[::-1]
    if positive.size < 32:
        return int(max(MIN_UMI_FLOOR, np.median(positive)))

    log_rank = np.log10(np.arange(1, positive.size + 1, dtype=np.float64))
    log_counts = np.log10(np.maximum(positive.astype(np.float64), 1.0))
    window = max(31, min(501, positive.size // 200))
    if window % 2 == 0:
        window += 1
    smooth_counts = (
        pd.Series(log_counts).rolling(window=window, center=True, min_periods=1).median().to_numpy()
    )
    first_derivative = np.gradient(smooth_counts, log_rank)
    second_derivative = np.gradient(first_derivative, log_rank)
    curvature = np.abs(second_derivative) / np.power(1.0 + np.square(first_derivative), 1.5)

    lo = max(10, int(positive.size * 0.0005))
    hi = min(positive.size - 1, max(lo + 1, int(positive.size * 0.05)))
    knee_idx = lo + int(np.argmax(curvature[lo:hi])) if hi > lo else int(np.argmax(curvature))
    return int(positive[knee_idx])


def build_called_barcode_mask(
    sample: str,
    mat,
    var_ref: pd.DataFrame,
    meta_row: pd.Series,
) -> tuple[np.ndarray, dict]:
    total_counts = np.asarray(mat.sum(axis=1)).ravel().astype(np.int64)
    detected_genes = np.diff(mat.indptr).astype(np.int32)
    positive_mask = total_counts > 0

    mito_mask = var_ref["gene_name"].astype(str).str.upper().str.startswith("MT-").to_numpy()
    row_max = compute_row_max(mat)
    top_gene_fraction = np.divide(
        row_max,
        total_counts,
        out=np.zeros(mat.shape[0], dtype=np.float32),
        where=total_counts > 0,
    )
    if np.any(mito_mask):
        mito_counts = np.asarray(mat[:, mito_mask].sum(axis=1)).ravel().astype(np.int64)
    else:
        mito_counts = np.zeros(mat.shape[0], dtype=np.int64)
    mito_fraction = np.divide(
        mito_counts,
        total_counts,
        out=np.zeros(mat.shape[0], dtype=np.float32),
        where=total_counts > 0,
    )

    knee_umi_threshold = estimate_knee_umi_threshold(total_counts[positive_mask])
    umi_threshold = int(np.clip(knee_umi_threshold, MIN_UMI_FLOOR, MAX_UMI_CAP))
    base_keep_mask = (
        positive_mask
        & (total_counts >= umi_threshold)
        & (detected_genes >= MIN_GENES_FLOOR)
    )
    qc_keep_mask = (
        base_keep_mask
        & (top_gene_fraction <= MAX_TOP_GENE_FRACTION)
        & (mito_fraction <= MAX_MITO_FRACTION)
    )

    qc_relaxed = False
    keep_mask = qc_keep_mask
    if int(np.count_nonzero(keep_mask)) < MIN_RETAINED_BARCODES_PER_SAMPLE:
        keep_mask = base_keep_mask
        qc_relaxed = True

    retained_total = int(np.count_nonzero(keep_mask))
    if retained_total == 0:
        raise ValueError(f"Cell-calling gate removed every barcode for sample {sample}")

    retained_counts = total_counts[keep_mask]
    retained_genes = detected_genes[keep_mask]
    retained_top_gene_fraction = top_gene_fraction[keep_mask]
    retained_mito_fraction = mito_fraction[keep_mask]

    summary = {
        "sample": sample,
        "patient": str(meta_row["patient"]),
        "condition": str(meta_row["condition"]),
        "disease": str(meta_row["disease"]),
        "raw_barcodes": int(mat.shape[0]),
        "positive_barcodes": int(np.count_nonzero(positive_mask)),
        "retained_barcodes": retained_total,
        "dropped_barcodes": int(mat.shape[0] - retained_total),
        "retention_fraction": safe_fraction(retained_total, int(mat.shape[0])),
        "barcodes_ge_100_umis": int(np.count_nonzero(total_counts >= 100)),
        "barcodes_ge_200_umis": int(np.count_nonzero(total_counts >= 200)),
        "barcodes_ge_500_umis": int(np.count_nonzero(total_counts >= 500)),
        "barcodes_ge_1000_umis": int(np.count_nonzero(total_counts >= 1_000)),
        "knee_umi_threshold": int(knee_umi_threshold),
        "umi_threshold_used": int(umi_threshold),
        "min_genes_threshold": MIN_GENES_FLOOR,
        "max_top_gene_fraction": MAX_TOP_GENE_FRACTION,
        "max_mito_fraction": MAX_MITO_FRACTION,
        "qc_relaxed": qc_relaxed,
        "median_umis_retained": float(np.median(retained_counts)),
        "median_genes_retained": float(np.median(retained_genes)),
        "median_top_gene_fraction_retained": float(np.median(retained_top_gene_fraction)),
        "median_mito_fraction_retained": float(np.median(retained_mito_fraction)),
        "mitochondrial_gene_count": int(np.count_nonzero(mito_mask)),
    }
    return keep_mask, summary


def build_prepare_summary(
    sample_summaries: list[dict],
    prepared_zarr: Path,
    expected_n_samples: int,
) -> dict:
    raw_total = sum(int(item["raw_barcodes"]) for item in sample_summaries)
    positive_total = sum(int(item["positive_barcodes"]) for item in sample_summaries)
    retained_total = sum(int(item["retained_barcodes"]) for item in sample_summaries)
    sample_summaries = sorted(sample_summaries, key=lambda item: item["sample"])
    return {
        "cell_calling_version": CELL_CALLING_VERSION,
        "method": "samplewise_hybrid_knee_plus_qc",
        "prepared_zarr": str(prepared_zarr),
        "expected_n_samples": int(expected_n_samples),
        "n_samples": int(len(sample_summaries)),
        "raw_barcodes_total": int(raw_total),
        "positive_barcodes_total": int(positive_total),
        "retained_barcodes_total": int(retained_total),
        "dropped_barcodes_total": int(raw_total - retained_total),
        "retention_fraction": safe_fraction(retained_total, raw_total),
        "gate_parameters": {
            "min_umi_floor": MIN_UMI_FLOOR,
            "max_umi_cap": MAX_UMI_CAP,
            "min_genes_floor": MIN_GENES_FLOOR,
            "max_top_gene_fraction": MAX_TOP_GENE_FRACTION,
            "max_mito_fraction": MAX_MITO_FRACTION,
            "min_retained_barcodes_per_sample": MIN_RETAINED_BARCODES_PER_SAMPLE,
            "max_validated_retention_fraction": MAX_VALIDATED_RETENTION_FRACTION,
        },
        "samples": sample_summaries,
    }


def validate_prepared_summary(
    summary_payload: dict,
    adata,
    expected_n_samples: int,
) -> dict:
    required_keys = {
        "cell_calling_version",
        "method",
        "expected_n_samples",
        "n_samples",
        "raw_barcodes_total",
        "retained_barcodes_total",
        "retention_fraction",
        "samples",
    }
    missing = sorted(required_keys - set(summary_payload))
    if missing:
        raise ValueError(f"Prepared summary missing required keys: {missing}")

    if int(summary_payload["expected_n_samples"]) != int(expected_n_samples):
        raise ValueError(
            "Prepared summary expected_n_samples mismatch: "
            f"expected {expected_n_samples}, got {summary_payload['expected_n_samples']}"
        )

    samples = list(summary_payload["samples"])
    if int(summary_payload["n_samples"]) != len(samples):
        raise ValueError("Prepared summary n_samples does not match samples payload length")
    if len(samples) != int(expected_n_samples):
        raise ValueError(
            f"Prepared summary sample count mismatch: expected {expected_n_samples}, got {len(samples)}"
        )

    raw_total = int(summary_payload["raw_barcodes_total"])
    retained_total = int(summary_payload["retained_barcodes_total"])
    if retained_total != int(getattr(adata, "n_obs", -1)):
        raise ValueError(
            "Prepared summary retained_barcodes_total mismatch: "
            f"summary={retained_total}, zarr={getattr(adata, 'n_obs', -1)}"
        )
    if raw_total <= retained_total:
        raise ValueError("Prepared summary indicates no effective barcode filtering")

    retention_fraction = float(summary_payload["retention_fraction"])
    if not (0.0 < retention_fraction <= MAX_VALIDATED_RETENTION_FRACTION):
        raise ValueError(
            "Prepared summary retention_fraction is implausible for a called-cell cohort: "
            f"{retention_fraction:.6f}"
        )

    obs_sample_counts = (
        pd.Series(adata.obs["sample"].astype(str)).value_counts().sort_index().to_dict()
    )
    for sample_summary in samples:
        sample = str(sample_summary["sample"])
        retained = int(sample_summary["retained_barcodes"])
        if retained <= 0:
            raise ValueError(f"Prepared summary retained zero barcodes for sample {sample}")
        observed = int(obs_sample_counts.get(sample, 0))
        if observed != retained:
            raise ValueError(
                f"Prepared summary sample mismatch for {sample}: summary={retained}, zarr={observed}"
            )

    return {
        "summary_cell_calling_version": str(summary_payload["cell_calling_version"]),
        "summary_method": str(summary_payload["method"]),
        "raw_barcodes_total": raw_total,
        "retained_barcodes_total": retained_total,
        "retention_fraction": retention_fraction,
    }


def write_part(
    sample: str,
    matrix_filename: str,
    mat,
    obs: pd.DataFrame,
    var_ref: pd.DataFrame,
    parts_root: Path,
) -> Path:
    var_out = make_var_out(var_ref)
    adata = ad.AnnData(X=mat.tocsr(), obs=obs, var=var_out)
    adata.var_names_make_unique()
    part_path = parts_root / f"{sample}.zarr"
    print(f"WRITE_PART {sample} -> {part_path}", flush=True)
    adata.write_zarr(part_path)
    return part_path


def open_prepared_input(prepared_zarr: Path):
    if hasattr(ad, "experimental") and hasattr(ad.experimental, "read_lazy"):
        return ad.experimental.read_lazy(prepared_zarr)
    return ad.read_zarr(prepared_zarr)


def validate_prepared_input(
    prepared_zarr: Path,
    expected_n_samples: int = EXPECTED_N_SAMPLES,
    summary_path: Path | None = None,
    require_summary: bool = False,
) -> dict:
    if not prepared_zarr.exists():
        raise FileNotFoundError(f"Prepared input missing: {prepared_zarr}")

    root = zarr.open_group(prepared_zarr, mode="r")
    if "X" not in root:
        raise ValueError(f"Prepared input missing X group: {prepared_zarr}")
    x_group = root["X"]
    x_attrs = dict(x_group.attrs)
    if x_attrs.get("encoding-type") != "csr_matrix":
        raise ValueError(
            f"Prepared input is not CSR-backed: encoding-type={x_attrs.get('encoding-type')}"
        )
    for key in ("data", "indices", "indptr"):
        if key not in x_group:
            raise ValueError(f"Prepared input missing X/{key}: {prepared_zarr}")

    adata = open_prepared_input(prepared_zarr)
    obs_columns = {str(col) for col in adata.obs.columns}
    missing = [col for col in REQUIRED_OBS_COLUMNS if col not in obs_columns]
    if missing:
        raise ValueError(f"Prepared input missing required obs columns: {missing}")

    n_samples = int(pd.Series(adata.obs["sample"].astype(str)).nunique())
    if n_samples != expected_n_samples:
        raise ValueError(
            f"Prepared input sample count mismatch: expected {expected_n_samples}, got {n_samples}"
        )

    summary_path = summary_path or (prepared_zarr.parent / SUMMARY_JSON.name)
    summary_validation = {}
    if summary_path.exists():
        summary_payload = json.loads(summary_path.read_text(encoding="utf-8"))
        summary_validation = validate_prepared_summary(
            summary_payload,
            adata,
            expected_n_samples=expected_n_samples,
        )
    elif require_summary:
        raise FileNotFoundError(f"Prepared summary missing: {summary_path}")

    payload = {
        "prepared_zarr": str(prepared_zarr),
        "shape": [int(getattr(adata, "n_obs", -1)), int(getattr(adata, "n_vars", -1))],
        "n_samples": n_samples,
        "required_obs_columns": REQUIRED_OBS_COLUMNS,
        "x_encoding": x_attrs.get("encoding-type"),
    }
    if summary_validation:
        payload["prepared_summary"] = {
            "path": str(summary_path),
            **summary_validation,
        }
    return payload


def prepare_full_cohort(
    root: Path = ROOT,
    expected_n_samples: int = EXPECTED_N_SAMPLES,
) -> dict:
    sdrf_path = root / SDRF.name
    sample_meta = build_sample_meta(sdrf_path).set_index("sample")
    matrix_files = sorted(root.glob("*-matrix.mtx.gz"))
    if not matrix_files:
        raise SystemExit("No matrix files found")
    if len(matrix_files) != expected_n_samples:
        raise ValueError(
            f"Expected {expected_n_samples} matrix files, found {len(matrix_files)}"
        )

    full_root = root / FULL_ROOT.name
    parts_root = full_root / PARTS_ROOT.name
    prepared_zarr = full_root / PREPARED_ZARR.name
    summary_json = full_root / SUMMARY_JSON.name
    ready_sentinel = full_root / READY_SENTINEL.name

    if parts_root.exists():
        shutil.rmtree(parts_root)
    full_root.mkdir(parents=True, exist_ok=True)
    parts_root.mkdir(parents=True, exist_ok=True)

    if prepared_zarr.exists():
        shutil.rmtree(prepared_zarr)
    if summary_json.exists():
        summary_json.unlink()
    if ready_sentinel.exists():
        ready_sentinel.unlink()

    var_ref = None
    part_paths: list[str] = []
    sample_summaries: list[dict] = []

    for matrix_path in matrix_files:
        sample = matrix_path.name.replace("-matrix.mtx.gz", "")
        print(f"READING {sample}", flush=True)
        feature_path = root / f"{sample}-features.tsv.gz"
        barcode_path = root / f"{sample}-barcodes.tsv.gz"

        if not feature_path.exists() or not barcode_path.exists():
            raise FileNotFoundError(f"Missing feature/barcode files for {sample}")
        if sample not in sample_meta.index:
            raise KeyError(f"Sample {sample} missing from SDRF metadata")

        var = read_features(feature_path)
        if var_ref is None:
            var_ref = var.copy()
        elif not var[["gene_id", "gene_name"]].equals(var_ref[["gene_id", "gene_name"]]):
            raise ValueError(f"Feature mismatch for sample {sample}")

        barcodes = read_barcodes(barcode_path).astype(str)
        mat = read_matrix(matrix_path)
        if mat.shape[0] != len(barcodes):
            raise ValueError(f"Barcode mismatch for {sample}: {mat.shape[0]} vs {len(barcodes)}")
        if mat.shape[1] != len(var_ref):
            raise ValueError(f"Feature mismatch in matrix for {sample}: {mat.shape[1]} vs {len(var_ref)}")

        meta_row = sample_meta.loc[sample]
        keep_mask, sample_summary = build_called_barcode_mask(sample, mat, var_ref, meta_row)
        barcodes = barcodes[keep_mask]
        mat = mat[keep_mask].tocsr()
        sample_summaries.append(sample_summary)

        cell_ids = pd.Index([f"{sample}:{bc}" for bc in barcodes], name="cell_id")
        obs = pd.DataFrame(index=cell_ids)
        obs["sample"] = sample
        obs["patient"] = meta_row["patient"]
        obs["batch"] = meta_row["batch"]
        obs["disease"] = meta_row["disease"]
        obs["condition"] = meta_row["condition"]
        obs["sorting"] = meta_row["sorting"]
        obs["sampling_site"] = meta_row["sampling_site"]
        obs["sex"] = meta_row["sex"]
        obs["original_source_name"] = meta_row["original_source_name"]
        obs["tumor_type"] = meta_row["tumor_type"]
        obs["donor"] = meta_row["patient"]
        obs["sample_id"] = sample
        obs["source_file"] = matrix_path.name
        obs["batch_join"] = meta_row["batch"]

        part_path = write_part(sample, matrix_path.name, mat, obs, var_ref, parts_root)
        part_paths.append(str(part_path))

    print(f"CONCAT_FULL_COHORT -> {prepared_zarr}", flush=True)
    concat_on_disk(part_paths, prepared_zarr, axis=0, join="outer", max_loaded_elems=50_000_000)
    # Ensure no spurious '_index' columns survive concat_on_disk (reserved by modern anndata).
    _merged = ad.read_zarr(prepared_zarr)
    _dirty = False
    for _df_attr in ("obs", "var"):
        _df = getattr(_merged, _df_attr)
        if "_index" in _df.columns:
            _df.drop(columns=["_index"], inplace=True)
            setattr(_merged, _df_attr, _df)
            _dirty = True
    if _merged.obs.index.name is None:
        _merged.obs.index.name = "cell_barcode"
        _dirty = True
    if _dirty:
        _merged.write_zarr(prepared_zarr)
    del _merged
    summary_payload = build_prepare_summary(
        sample_summaries,
        prepared_zarr=prepared_zarr,
        expected_n_samples=expected_n_samples,
    )
    summary_json.write_text(json.dumps(summary_payload, indent=2, ensure_ascii=False), encoding="utf-8")
    validation = validate_prepared_input(
        prepared_zarr,
        expected_n_samples=expected_n_samples,
        summary_path=summary_json,
        require_summary=True,
    )
    ready_sentinel.write_text("ok\n", encoding="utf-8")
    print("CONCAT_DONE_FULL_COHORT", flush=True)
    return validation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare and validate the NC2024 full-cohort zarr input."
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate the existing merged zarr without rebuilding it.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.validate_only:
        payload = validate_prepared_input(
            PREPARED_ZARR,
            expected_n_samples=EXPECTED_N_SAMPLES,
            summary_path=SUMMARY_JSON,
            require_summary=True,
        )
    else:
        payload = prepare_full_cohort(ROOT, expected_n_samples=EXPECTED_N_SAMPLES)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
