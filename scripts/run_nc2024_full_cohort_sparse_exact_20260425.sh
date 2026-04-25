#!/usr/bin/env bash
# NC2024 NSCLC 900k full-cohort REAL run — sparse-exact engine path.
#
# Replaces the CSS approximation in clustering with sparse-exact PCA+Leiden,
# plus Phase 2B sparse Welch DE / sparse_groupby cell_communication / chunked CNV.
# Phase 3 MemoryGuard observes memory pressure but never aborts.
#
# Compared to scripts/run_nc2024_extended_full_cohort_20260424.sh:
#   - --scale-mode large  (was massive — CSS path)
#   - SC_CLUSTERING_ENGINE=sparse_exact, SC_DE_ENGINE=sparse,
#     SC_CNV_ENGINE=chunked, SC_CELLCOMM_ENGINE=sparse, SC_MEM_GUARD=on
#   - SC_DENSIFY_HARD_CAP_GB=12 (default 8 raised — workstation has 96GB)
#   - SCF_MASSIVE_CHECKPOINT_POLICY removed (no longer in massive mode)
#
# Failure recovery: --checkpoint allows resume. If sparse_exact clustering
# OOMs on the workstation, rerun with SC_CLUSTERING_ENGINE=  (empty) to fall
# back to massive/CSS while keeping the Phase 2B DE/CNV/cell_comm gains.
set -euo pipefail

REPO_ROOT="/home/zerlinshen/singlecell_factory"
RAW_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
SAMPLE_ROOT="${RAW_ROOT}/full_cohort"
RESULTS="${REPO_ROOT}/results"
PROJECT="NC2024_NSCLC_FULL_COHORT_SPARSE_EXACT_REAL_AUTO"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"
PAPER_SPEC="${REPO_ROOT}/ops/nc2024_sparse_exact/paper_repro_spec.json"
PROJECT_LOG="${RESULTS}/${PROJECT}.launch.log"
BEFORE_FILE="$(mktemp)"
AFTER_FILE="$(mktemp)"

MODULES="clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,pathway_analysis,cell_communication,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,pseudobulk_de,cell_fate,composition,metacell,paper_repro,validate_cbioportal"
CBIO_GENES="ELF3,EGFR,KRAS,TP53,CD274,PDCD1,CTLA4,TIGIT,HAVCR2,CD96,NECTIN1,CD80,CD86"

mkdir -p "${RESULTS}" "$(dirname "${PAPER_SPEC}")"
cat > "${PAPER_SPEC}" <<'JSON'
{
  "papers": [
    {
      "paper_id": "nc2024_nsclc_single_cell",
      "title": "Nature Communications 2024 NSCLC single-cell reproduction",
      "source": {
        "paper_path": "",
        "repo_url": "local-reproduction-workspace",
        "repo_commit": "workspace-artifact",
        "license": "research-reproduction"
      },
      "adaptation_targets": [
        "clustering",
        "annotation",
        "composition",
        "immune_phenotyping",
        "tumor_microenvironment",
        "cell_communication",
        "gene_signature_scoring",
        "pseudobulk_de",
        "metacell"
      ],
      "reproductions": [],
      "notes": "Sparse-exact engine path (Phase 2B). Clustering uses arpack PCA + scanpy neighbors + Leiden iGraph instead of CSS approximation; DE uses sparse Welch instead of toarray fallback; CNV uses chunked sparse smoothing; cell communication uses sparse_groupby_mean."
    }
  ]
}
JSON

: > "${PROJECT_LOG}"
find "${RESULTS}" -maxdepth 1 -type d -name "${PROJECT}_*" | sort > "${BEFORE_FILE}"

set +e
cd "${REPO_ROOT}"
SC_CLUSTERING_ENGINE=sparse_exact \
SC_DE_ENGINE=sparse \
SC_CNV_ENGINE=chunked \
SC_CELLCOMM_ENGINE=sparse \
SC_MEM_GUARD=on \
SC_MEM_WATCHDOG=on \
SC_DENSIFY_HARD_CAP_GB=12 \
SC_DENSIFY_SOFT_CAP_GB=6 \
CONDA_PREFIX=/home/zerlinshen/conda/envs/sc_gpu \
${CONDA_RUN} python -m workflow.modular.cli \
  --project "${PROJECT}" \
  --sample-root "${SAMPLE_ROOT}" \
  --output-dir "${RESULTS}" \
  --optional-modules "${MODULES}" \
  --batch-key sample \
  --checkpoint \
  --gpu-mode auto \
  --scale-mode large \
  --parallel-workers 1 \
  --pseudobulk-exploratory-group-vs-rest \
  --paper-spec-json "${PAPER_SPEC}" \
  --cbioportal-genes "${CBIO_GENES}" \
  --cbioportal-study lusc_tcga_pan_can_atlas_2018 \
  --cbioportal-top-n 10 \
  --no-run-cellranger-if-missing > "${PROJECT_LOG}" 2>&1
EC=$?
set -e

find "${RESULTS}" -maxdepth 1 -type d -name "${PROJECT}_*" | sort > "${AFTER_FILE}"
RUN_DIR="$(comm -13 "${BEFORE_FILE}" "${AFTER_FILE}" | tail -n 1 || true)"
if [[ -z "${RUN_DIR}" && ${EC} -eq 0 ]]; then
  RUN_DIR="$(tail -n 1 "${AFTER_FILE}" || true)"
fi
rm -f "${BEFORE_FILE}" "${AFTER_FILE}"
echo "RUN_DIR=${RUN_DIR}"
echo "EC=${EC}"
exit ${EC}
