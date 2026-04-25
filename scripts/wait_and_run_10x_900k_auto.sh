#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/zerlinshen/singlecell_factory
DATA=$ROOT/data/raw/10x_nsclc_900k_flex
TAR=$DATA/16plex_900k_32_NSCLC_multiplex_count_filtered_feature_bc_matrix.tar.gz
DONE=$DATA/.download_complete
LOG=$DATA/auto_run_900k.log
RUNLOG=$ROOT/results/NSCLC_900K_BASELINE_AUTO.launch.log
{
  echo "[$(date)] waiting for completed 900k download"
  until [[ -f "$DONE" ]] || { [[ -s "$TAR" ]] && tar -tzf "$TAR" >/dev/null 2>&1; }; do
    sleep 300
  done
  echo "[$(date)] 900k tarball verified complete"
  bash $ROOT/scripts/prepare_10x_900k_layout.sh
  echo "[$(date)] extracted matrix layout"
  cd $ROOT
  env MPLCONFIGDIR=$ROOT/.mplconfig NUMBA_CACHE_DIR=/tmp/numba_cache /home/zerlinshen/conda/bin/conda run -n sc_gpu \
    python -m workflow.modular.cli \
      --project NSCLC_900K_BASELINE_AUTO \
      --sample-root $DATA \
      --optional-modules clustering,differential_expression,annotation,immune_phenotyping,tumor_microenvironment,metacell \
      --checkpoint \
      --gpu-mode auto \
      --n-top-genes 2000 \
      --n-pcs 30 \
      --n-neighbors 15 \
      --leiden-resolution 0.6 \
      --expected-doublet-rate 0.03 > "$RUNLOG" 2>&1
  echo "[$(date)] 900k baseline run finished"
} >> "$LOG" 2>&1
