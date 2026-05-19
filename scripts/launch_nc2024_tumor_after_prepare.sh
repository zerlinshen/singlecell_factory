#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
RESULTS="${REPO_ROOT}/results"
LOG="${RESULTS}/NC2024_NSCLC_TUMOR_PIPELINE.launch.log"
while true; do
  if [ -f "$RAW/tumor/prepared_input.ready" ]; then
    break
  fi
  sleep 30
done
cd "${REPO_ROOT}"
/home/zerlinshen/conda/bin/conda run -n sc_gpu python -m workflow.modular.cli \
  --project NC2024_NSCLC_TUMOR_REPRO_AUTO \
  --sample-root "$RAW/tumor" \
  --optional-modules clustering,batch_correction,differential_expression,annotation,composition,immune_phenotyping,tumor_microenvironment,trajectory,cnv_inference \
  --batch-key sample \
  --checkpoint \
  --gpu-mode auto \
  --scale-mode large > "$LOG" 2>&1
