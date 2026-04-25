#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="/home/zerlinshen/singlecell_factory"
RAW_ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
SAMPLE_ROOT="${RAW_ROOT}/full_cohort"
RESULTS="${REPO_ROOT}/results"
PROJECT="NC2024_NSCLC_FULL_COHORT_EXTENDED_MASSIVE_REAL_AUTO"
CONDA_RUN="/home/zerlinshen/conda/bin/conda run -n sc_gpu"
PAPER_SPEC="${REPO_ROOT}/ops/nc2024_extended_real_run/paper_repro_spec.json"
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
      "notes": "Non-strict provenance ledger for the 2026-04-24 NC2024 full-cohort extended real run. Figure-level comparison is packaged locally from remote outputs and existing paper boards."
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
exit "${EC}"
