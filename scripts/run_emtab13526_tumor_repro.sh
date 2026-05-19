#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${REPO_ROOT}/data/raw/nc2024_nsclc_emtab13526"
BASEURL=https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/526/E-MTAB-13526/Files/
mkdir -p "$ROOT"
cd "$ROOT"
curl -L --max-time 60 "$BASEURL" \
  | grep -Eo 'href="[^"]+"' \
  | sed 's/href="//;s/"//' \
  | grep -E '(^[^/]+-matrix\.mtx\.gz$|^[^/]+-barcodes\.tsv\.gz$|^[^/]+-features\.tsv\.gz$|^E-MTAB-13526\.(sdrf|idf)\.txt$)' \
  | sed "s#^#$BASEURL#" > download_urls.txt
wget --continue --tries=10 --waitretry=5 --read-timeout=60 --timeout=60 -i download_urls.txt
/home/zerlinshen/conda/bin/conda run -n sc_gpu python "${REPO_ROOT}/scripts/prepare_emtab13526_input.py"
cd "${REPO_ROOT}"
/home/zerlinshen/conda/bin/conda run -n sc_gpu python -m workflow.modular.cli \
  --project NC2024_NSCLC_TUMOR_REPRO_AUTO \
  --sample-root "${ROOT}/tumor" \
  --optional-modules clustering,batch_correction,differential_expression,annotation,composition,immune_phenotyping,tumor_microenvironment,trajectory,cnv_inference \
  --batch-key sample \
  --checkpoint \
  --gpu-mode auto \
  --scale-mode large
