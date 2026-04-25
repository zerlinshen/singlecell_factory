#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/zerlinshen/singlecell_factory
printf '== host ==\n'
hostname
uptime
printf '\n== memory ==\n'
free -h
printf '\n== swap ==\n'
swapon --show || true
printf '\n== sizes ==\n'
du -sh "$ROOT/data/raw/10x_nsclc_900k_flex" "$ROOT/data/raw/10x_nsclc_40k_dtc_7donors" 2>/dev/null || true
printf '\n== active processes ==\n'
ps -ef | grep -E 'NSCLC_900K_MASSIVE_V2_AUTO|40k_NSCLC_DTC_3p_HT_nextgem_donor_|remote_watch_900k_40k' | grep -v grep || true
printf '\n== latest monitor ==\n'
tail -n 20 "$ROOT/runtime_monitor/monitor_900k_40k.log" 2>/dev/null || true
printf '\n== latest 900k run ==\n'
tail -n 20 "$ROOT/results/NSCLC_900K_MASSIVE_V2_AUTO.launch.log" 2>/dev/null || true
