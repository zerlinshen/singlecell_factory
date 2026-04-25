#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/zerlinshen/singlecell_factory"
RAW_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
SAMPLE_ROOT="${RAW_ROOT}/full_cohort"
RESULTS="${REPO_ROOT}/results"
PROJECT="NC2024_NSCLC_FULL_COHORT_RERUN_ALL_ELIGIBLE_AUTO"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"
PAPER_SPEC="${REPO_ROOT}/ops/nc2024_full_cohort_rerun_all_eligible/paper_repro_spec.json"
PROJECT_LOG="${RESULTS}/${PROJECT}.launch.log"
BEFORE_FILE="$(mktemp)"
AFTER_FILE="$(mktemp)"

MODULES="clustering,cell_cycle,batch_correction,differential_expression,annotation,trajectory,pseudo_velocity,pathway_analysis,cell_communication,immune_phenotyping,tumor_microenvironment,gene_signature_scoring,pseudobulk_de,cell_fate,composition,metacell,paper_repro"

mkdir -p "${RESULTS}" "$(dirname "${PAPER_SPEC}")"
cat > "${PAPER_SPEC}" <<'JSON'
{
  "papers": [
    {
      "paper_id": "nc2024_nsclc_single_cell_rerun_all_eligible",
      "title": "Nature Communications 2024 NSCLC single-cell rerun after pseudobulk repair",
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
        "metacell",
        "paper_repro"
      ],
      "reproductions": [],
      "notes": "Clean full-cohort rerun with all eligible modules only. validate_cbioportal excluded for determinism; rna_velocity/cnv_inference/evolution excluded for modality or scale-safety reasons. pseudobulk_de runs in explicit exploratory mode."
    }
  ]
}
JSON

: > "${PROJECT_LOG}"
find "${RESULTS}" -maxdepth 1 -type d -name "${PROJECT}_*" | sort > "${BEFORE_FILE}"

set +e
cd "${REPO_ROOT}"
SCF_MASSIVE_CHECKPOINT_POLICY=metadata_only \
CONDA_PREFIX=/home/zerlinshen/conda/envs/sc_gpu \
${CONDA_RUN} python -m workflow.modular.cli \
  --project "${PROJECT}" \
  --sample-root "${SAMPLE_ROOT}" \
  --output-dir "${RESULTS}" \
  --optional-modules "${MODULES}" \
  --batch-key sample \
  --checkpoint \
  --gpu-mode auto \
  --scale-mode massive \
  --parallel-workers 1 \
  --pseudobulk-exploratory-group-vs-rest \
  --paper-spec-json "${PAPER_SPEC}" \
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
exit "${EC}"
