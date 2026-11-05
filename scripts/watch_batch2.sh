#!/bin/bash
# Periodic progress line for Auriga batch 2. Read-only.
#   <search pts>/257  s=merged  R=rescored  E=scheme-B posterior  A=archived
set -u
cd "$(dirname "$0")"; ROOT=$(cd .. && pwd)
U=${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}
HALOS=${HALOS:-"5 12 14 20 22 28"}; EVERY_MIN=${EVERY_MIN:-5}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
fg () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }
last=-1; stall=0
while :; do
  tot=0; fin=0; line=""
  for H in $HALOS; do
    R=test2_au${H}_sf_12p5pct; d=$ROOT/runs/$R
    n=$(find "$d" -name 'fourier_*.npz' 2>/dev/null | wc -l); tot=$((tot+n))
    f=""
    [ -f "$ROOT/results/$R/evals.csv" ] && f="${f}s"
    [ "$(find "$d" -name 'sfratio_*.npz' 2>/dev/null | wc -l)" -ge 257 ] && f="${f}R"
    [ -f "$ROOT/results/$R/emulator_results_logLs_${B}_sfratio.json" ] && { f="${f}E"; }
    [ -d "$U/runs/$R" ] && { f="${f}A"; fin=$((fin+1)); }
    line="$line au${H}:${n}${f}"
  done
  [ "$tot" -eq "$last" ] && stall=$((stall+EVERY_MIN)) || stall=0; last=$tot
  printf '%s  %4d/1542 search pts  %d/6 complete |%s | %sG here, %sG usb' \
     "$(date '+%F %T')" "$tot" "$fin" "$line" "$(fg "$ROOT/runs")" "$(fg "$U")"
  [ "$stall" -ge 45 ] && printf '  ** no search progress for %d min **' "$stall"
  printf '\n'
  [ "$fin" -ge 6 ] && { echo "BATCH 2 COMPLETE"; break; }
  sleep $((EVERY_MIN*60))
done
