#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/zerlinshen/singlecell_factory
DATA900K=$ROOT/data/raw/10x_nsclc_900k_flex
DATA40K=$ROOT/data/raw/10x_nsclc_40k_dtc_7donors
RESULTS=$ROOT/results
STATE=$ROOT/runtime_monitor
RUN_NAME=NSCLC_900K_CSS_PROXY_V3_AUTO
LOG=$STATE/monitor_900k_40k.log
LAUNCH_LOCK=$STATE/${RUN_NAME}.launch_attempted.flag
RUNLOG=$RESULTS/${RUN_NAME}.launch.log
mkdir -p "$STATE" "$RESULTS"

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >> "$LOG"; }

start_900k_if_ready() {
  local prepared="$DATA900K/prepared_input.h5ad"
  if [[ ! -s "$prepared" ]]; then
    return 0
  fi
  if ps -ef | grep -E "python -m workflow.modular.cli .*${RUN_NAME}" | grep -v grep >/dev/null 2>&1; then
    return 0
  fi
  if ls "$RESULTS" | grep -E "^${RUN_NAME}(_|$)" >/dev/null 2>&1; then
    return 0
  fi
  if [[ -f "$LAUNCH_LOCK" ]]; then
    return 0
  fi
  touch "$LAUNCH_LOCK"
  log "prepared_input.h5ad is ready; launching ${RUN_NAME}"
  nohup bash -lc "cd '$ROOT' && env MPLCONFIGDIR='$ROOT/.mplconfig' NUMBA_CACHE_DIR=/tmp/numba_cache /home/zerlinshen/conda/bin/conda run -n sc_gpu python -m workflow.modular.cli --project '${RUN_NAME}' --sample-root '$DATA900K' --optional-modules clustering --checkpoint --gpu-mode auto --n-top-genes 1000 --n-pcs 20 --n-neighbors 10 --leiden-resolution 0.4 --expected-doublet-rate 0.03 --scale-mode massive" > "$RUNLOG" 2>&1 &
  log "launched ${RUN_NAME} pid=$!"
}

summarize_once() {
  local size900k size40k prep_size
  size900k=$(du -sh "$DATA900K" 2>/dev/null | awk '{print $1}')
  size40k=$(du -sh "$DATA40K" 2>/dev/null | awk '{print $1}')
  if [[ -f "$DATA900K/prepared_input.h5ad" ]]; then
    prep_size=$(du -sh "$DATA900K/prepared_input.h5ad" 2>/dev/null | awk '{print $1}')
  else
    prep_size='missing'
  fi
  log "status 900k_dir=$size900k prepared_input=$prep_size 40k_dir=$size40k"

  local proc40k proc900k
  proc40k=$(ps -ef | grep -E '40k_NSCLC_DTC_3p_HT_nextgem_donor_[0-9]+_count_sample_alignments.bam|40k_NSCLC_DTC_3p_HT_nextgem_donor_[0-9]+_count_' | grep -v grep | head -n 1 || true)
  proc900k=$(ps -ef | grep -E "${RUN_NAME}|prepare_900k_from_h5" | grep -v grep | head -n 2 || true)
  [[ -n "$proc40k" ]] && log "40k_active: $proc40k"
  [[ -n "$proc900k" ]] && log "900k_active: $proc900k"

  if [[ -f "$DATA40K/download.log" ]]; then
    log "40k_tail: $(tail -n 1 "$DATA40K/download.log" | tr '\0' ' ' | sed 's/[[:space:]]\+/ /g')"
  fi
  if [[ -f "$DATA900K/prepare_from_h5.log" ]]; then
    log "900k_prepare_tail: $(tail -n 1 "$DATA900K/prepare_from_h5.log" | tr '\0' ' ' | sed 's/[[:space:]]\+/ /g')"
  fi
  if [[ -f "$RUNLOG" ]]; then
    log "900k_run_tail: $(tail -n 1 "$RUNLOG" | tr '\0' ' ' | sed 's/[[:space:]]\+/ /g')"
  fi
}

log 'monitor started'
while true; do
  start_900k_if_ready || log 'launch-check failed'
  summarize_once || log 'summary failed'
  sleep 300
done
