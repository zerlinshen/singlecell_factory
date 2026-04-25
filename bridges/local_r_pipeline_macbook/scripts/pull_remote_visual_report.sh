#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <remote_result_dir> <local_tag>" >&2
  exit 1
fi
REMOTE_RESULT_DIR="$1"
LOCAL_TAG="$2"
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_DIR="$ROOT/output/${LOCAL_TAG}_remote_report_assets"
mkdir -p "$LOCAL_DIR"
for rel in \
  run_manifest.json \
  module_status.csv \
  clustering/umap_leiden.png \
  annotation/umap_cell_type.png \
  annotation/umap_annotation_confidence.png \
  annotation/cell_type_composition.png \
  composition/composition_barplot.png \
  composition/composition_boxplot.png \
  immune_phenotyping/umap_immune_subtype.png \
  immune_phenotyping/immune_signature_heatmap.png \
  tumor_microenvironment/tme_cyt_umap.png \
  tumor_microenvironment/tme_tis_umap.png \
  tumor_microenvironment/tme_signature_heatmap.png \
  qc/qc_violin_post_filter.png \
  qc/qc_scatter_post_filter.png
  do
  mkdir -p "$LOCAL_DIR/$(dirname "$rel")"
  scp "ubuntu-tail:${REMOTE_RESULT_DIR}/${rel}" "$LOCAL_DIR/${rel}" 2>/dev/null || true
 done
 echo "$LOCAL_DIR"
