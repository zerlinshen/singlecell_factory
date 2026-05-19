#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <remote_result_dir> <remote_out_dir> [group_by] [cluster_by] [max_cells] [bundle_dir]" >&2
  exit 1
fi

REMOTE_RESULT_DIR="${1%/}"
REMOTE_OUT_DIR="$2"
GROUP_BY="${3:-cell_type}"
CLUSTER_BY="${4:-leiden}"
MAX_CELLS="${5:-200000}"
BUNDLE_DIR="${6:-${REMOTE_RESULT_DIR}/r_bundle}"
FINAL_ADATA="${REMOTE_RESULT_DIR}/final_adata.h5ad"
MANIFEST_TSV="${BUNDLE_DIR}/bundle_manifest.tsv"
MANIFEST_JSON="${BUNDLE_DIR}/bundle_manifest.json"

FACTORY_ROOT="${FACTORY_ROOT:-/home/zerlinshen/Bioinformatics Research Pipeline/singlecell_factory}"
CONDA_BIN="${CONDA_BIN:-/home/zerlinshen/conda/bin/conda}"
PY_ENV="${PY_ENV:-sc_gpu}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/home/zerlinshen/conda/envs/r_multiomics_arrow/bin/Rscript}"
PLOT_SCRIPT="${PLOT_SCRIPT:-/home/zerlinshen/Bioinformatics Research Pipeline/r_multiomics_factory/scripts/plot_remote_bundle_large.R}"
FORCE_R_BUNDLE_EXPORT="${FORCE_R_BUNDLE_EXPORT:-0}"
R_PLOT_THREADS="${R_PLOT_THREADS:-8}"
R_BUNDLE_MARKERS="${R_BUNDLE_MARKERS:-}"
R_BUNDLE_OBS_COLS="${R_BUNDLE_OBS_COLS:-}"
R_BUNDLE_OBSM="${R_BUNDLE_OBSM:-}"
R_REUSE_CHECK_V1='source("R/remote_bundle_manifest.R"); '
R_REUSE_CHECK_V1+='invisible(validate_remote_bundle(Sys.getenv("BUNDLE_DIR_FOR_R"), '
R_REUSE_CHECK_V1+='required_stems=c("obs","X_umap","X_pca","marker_expr"), '
R_REUSE_CHECK_V1+='allow_missing_manifest=FALSE)); cat("bundle-v1-reuse-integrity-ok\n")'
R_REUSE_CHECK_V2='source("R_bundle/io_bundle.R"); '
R_REUSE_CHECK_V2+='bundle <- read_bundle(Sys.getenv("BUNDLE_DIR_FOR_R")); '
R_REUSE_CHECK_V2+='missing <- setdiff(c("X_umap","X_pca"), names(bundle$obsm)); '
R_REUSE_CHECK_V2+='if (length(missing)) stop(sprintf("missing obsm: %s", paste(missing, collapse=","))); '
R_REUSE_CHECK_V2+='cat("bundle-v2-reuse-integrity-ok\n")'

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$R_PLOT_THREADS}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-$R_PLOT_THREADS}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-$R_PLOT_THREADS}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-$R_PLOT_THREADS}"
export R_DATATABLE_NUM_THREADS="${R_DATATABLE_NUM_THREADS:-$R_PLOT_THREADS}"

if [[ ! -x "$CONDA_BIN" ]]; then
  echo "Conda executable not found: $CONDA_BIN" >&2
  exit 1
fi

if [[ ! -x "$RSCRIPT_BIN" ]]; then
  echo "Remote Rscript not found: $RSCRIPT_BIN" >&2
  exit 1
fi

if [[ ! -f "$PLOT_SCRIPT" ]]; then
  echo "Remote R plot script not found: $PLOT_SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$FINAL_ADATA" ]]; then
  echo "Missing final_adata.h5ad in remote result dir: $REMOTE_RESULT_DIR" >&2
  exit 1
fi

manifest_value() {
  local key="$1"
  awk -F '\t' -v k="$key" '$1 == k { print $2; exit }' "$MANIFEST_TSV" | tr -d '\r'
}

manifest_json_scalar() {
  python3 - "$MANIFEST_JSON" "$@" <<'PY'
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as handle:
    value = json.load(handle)

for key in sys.argv[2:]:
    if not isinstance(value, dict):
        value = None
        break
    value = value.get(key)

if value is None:
    print("")
elif isinstance(value, float) and value.is_integer():
    print(int(value))
else:
    print(value)
PY
}

manifest_json_contains_all() {
  local requested="$1"
  shift
  [[ -z "$requested" ]] && return 0
  python3 - "$MANIFEST_JSON" "$requested" "$@" <<'PY'
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as handle:
    value = json.load(handle)

requested = [item.strip() for item in sys.argv[2].split(",") if item.strip()]
for key in sys.argv[3:]:
    if not isinstance(value, dict):
        value = []
        break
    value = value.get(key, [])

available = {str(item) for item in (value or [])}
missing = [item for item in requested if item not in available]
if missing:
    raise SystemExit("missing requested manifest values: " + ",".join(missing))
PY
}

current_final_adata_bytes() {
  stat -c '%s' "$FINAL_ADATA"
}

current_final_adata_mtime_epoch() {
  stat -c '%Y' "$FINAL_ADATA"
}

matches_optional_manifest_value() {
  local key="$1"
  local requested="$2"
  [[ -z "$requested" ]] && return 0
  [[ "$(manifest_value "$key")" == "$requested" ]]
}

csv_union() {
  python3 - "$@" <<'PY'
import sys

seen = set()
values = []
for raw in sys.argv[1:]:
    for item in raw.split(","):
        value = item.strip()
        if value and value not in seen:
            seen.add(value)
            values.append(value)
print(",".join(values))
PY
}

EFFECTIVE_R_BUNDLE_OBS_COLS="$R_BUNDLE_OBS_COLS"
if [[ -n "$R_BUNDLE_OBS_COLS" ]]; then
  EFFECTIVE_R_BUNDLE_OBS_COLS="$(csv_union "$GROUP_BY,$CLUSTER_BY" "$R_BUNDLE_OBS_COLS")"
fi

EFFECTIVE_R_BUNDLE_OBSM="$R_BUNDLE_OBSM"
if [[ -n "$R_BUNDLE_OBSM" ]]; then
  EFFECTIVE_R_BUNDLE_OBSM="$(csv_union "X_umap,X_pca" "$R_BUNDLE_OBSM")"
fi

can_reuse_v1_bundle() {
  [[ "$FORCE_R_BUNDLE_EXPORT" != "1" ]] || return 1
  [[ -f "$MANIFEST_TSV" ]] || return 1
  [[ "$MANIFEST_TSV" -nt "$FINAL_ADATA" ]] || return 1
  [[ "$(manifest_value input_h5ad)" == "$FINAL_ADATA" ]] || return 1
  [[ "$(manifest_value source_run_dir)" == "$REMOTE_RESULT_DIR" ]] || return 1
  [[ "$(manifest_value input_h5ad_bytes)" == "$(current_final_adata_bytes)" ]] || return 1
  [[ "$(manifest_value input_h5ad_mtime_epoch)" == "$(current_final_adata_mtime_epoch)" ]] || return 1
  matches_optional_manifest_value markers_requested "$R_BUNDLE_MARKERS" || return 1
  matches_optional_manifest_value obs_columns "$EFFECTIVE_R_BUNDLE_OBS_COLS" || return 1
  matches_optional_manifest_value obsm_keys "$EFFECTIVE_R_BUNDLE_OBSM" || return 1

  (
    cd "$FACTORY_ROOT/bridges/local_r_pipeline_macbook"
    BUNDLE_DIR_FOR_R="$BUNDLE_DIR" "$RSCRIPT_BIN" -e "$R_REUSE_CHECK_V1" >/dev/null
  )
}

manifest_json_files_present() {
  # C2 reader-side: refuse to reuse a partially written bundle. Walks every
  # entry in manifest$files, resolves the `path` against $BUNDLE_DIR, and
  # exits non-zero if anything the manifest references is absent on disk.
  # This is a presence sweep, deliberately cheaper than the per-file SHA256
  # check that runs later inside the R reuse probe.
  python3 - "$MANIFEST_JSON" "$BUNDLE_DIR" <<'PY'
import json
import os
import sys

manifest_path = sys.argv[1]
bundle_dir = sys.argv[2]
with open(manifest_path, "r", encoding="utf-8") as handle:
    manifest = json.load(handle)

files = manifest.get("files") or {}
if not isinstance(files, dict):
    raise SystemExit("manifest$files is not an object")

missing = []
for stem, rec in files.items():
    if not isinstance(rec, dict):
        continue
    rel = rec.get("path")
    if not rel:
        continue
    full = os.path.join(bundle_dir, rel)
    if not os.path.exists(full):
        missing.append(rel)

if missing:
    raise SystemExit("missing manifest files on disk: " + ",".join(missing))
PY
}

can_reuse_v2_bundle() {
  [[ "$FORCE_R_BUNDLE_EXPORT" != "1" ]] || return 1
  [[ -f "$MANIFEST_JSON" ]] || return 1
  [[ "$MANIFEST_JSON" -nt "$FINAL_ADATA" ]] || return 1
  [[ "$(manifest_json_scalar schema_version)" == "singlecell_r_bundle_v2" ]] || return 1
  [[ "$(manifest_json_scalar source input_h5ad)" == "$FINAL_ADATA" ]] || return 1
  [[ "$(manifest_json_scalar source input_h5ad_bytes)" == "$(current_final_adata_bytes)" ]] || return 1
  [[ "$(manifest_json_scalar source input_h5ad_mtime_epoch)" == "$(current_final_adata_mtime_epoch)" ]] || return 1
  manifest_json_contains_all "$R_BUNDLE_MARKERS" bundle markers_present || return 1
  manifest_json_contains_all "$EFFECTIVE_R_BUNDLE_OBS_COLS" bundle obs_columns_present || return 1
  manifest_json_contains_all "$EFFECTIVE_R_BUNDLE_OBSM" bundle obsm_keys || return 1
  # C2 reader-side guard: do a fast presence sweep BEFORE the per-file SHA256
  # check below so we treat partial bundles as "do not reuse, regenerate."
  manifest_json_files_present || return 1

  (
    cd "$FACTORY_ROOT/bridges/local_r_pipeline_macbook"
    BUNDLE_DIR_FOR_R="$BUNDLE_DIR" "$RSCRIPT_BIN" -e "$R_REUSE_CHECK_V2" >/dev/null
  )
}

can_reuse_bundle() {
  if [[ -f "$MANIFEST_JSON" ]]; then
    can_reuse_v2_bundle
    return
  fi
  can_reuse_v1_bundle
}

export_args=()
if [[ -n "$R_BUNDLE_MARKERS" ]]; then
  export_args+=(--markers "$R_BUNDLE_MARKERS")
fi
if [[ -n "$R_BUNDLE_OBS_COLS" ]]; then
  export_args+=(--obs-cols "$EFFECTIVE_R_BUNDLE_OBS_COLS")
fi
if [[ -n "$R_BUNDLE_OBSM" ]]; then
  export_args+=(--obsm "$EFFECTIVE_R_BUNDLE_OBSM")
fi

cd "$FACTORY_ROOT"
if can_reuse_bundle; then
  echo "Reusing validated R bundle: $BUNDLE_DIR"
else
  cd "$FACTORY_ROOT"
  "$CONDA_BIN" run -n "$PY_ENV" python scripts/export_singlecell_r_bundle.py \
    --input "$FINAL_ADATA" \
    --source-run-dir "$REMOTE_RESULT_DIR" \
    --output "$BUNDLE_DIR" \
    "${export_args[@]}"
fi

cd "$FACTORY_ROOT/bridges/local_r_pipeline_macbook"
"$RSCRIPT_BIN" "$PLOT_SCRIPT" \
  "$BUNDLE_DIR" \
  "$REMOTE_OUT_DIR" \
  "$GROUP_BY" \
  "$CLUSTER_BY" \
  "$MAX_CELLS"

echo "Remote R plots saved to: $REMOTE_OUT_DIR"
