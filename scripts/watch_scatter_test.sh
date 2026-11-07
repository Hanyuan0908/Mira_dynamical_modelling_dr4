#!/bin/bash
# Progress line for the per-star distance-error test.
#   <search pts>/257  s=merged  R=rescored(sfratio,K=10)  E=schemeB posterior  A=archived
set -u
cd "$(dirname "$0")"; ROOT=$(cd .. && pwd)
U=${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}
TAGS=${TAGS:-"t3 t5.5 t7 au18 au26"}; EVERY_MIN=${EVERY_MIN:-5}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
fg(){ df -BG --output=avail "$1" 2>/dev/null|tail -1|tr -dc '0-9'; }
last=-1; stall=0
while :; do
  tot=0; fin=0; line=""
  for T in $TAGS; do
    R=test2_${T}_sf_12p5pct_scat; d=$ROOT/runs/$R; f=""
    n=$(find "$d" -name 'fourier_*.npz' 2>/dev/null|wc -l); tot=$((tot+n))
    [ -f "$ROOT/results/$R/evals.csv" ] && f="${f}s"
    [ "$(find "$d" -name 'sfratio_*.npz' 2>/dev/null|wc -l)" -ge 257 ] && f="${f}R"
    [ -f "$ROOT/results/$R/emulator_results_logLs_${B}_sfratio.json" ] && f="${f}E"
    [ -d "$U/runs/$R" ] && { f="${f}A"; fin=$((fin+1)); }
    line="$line ${T}:${n}${f}"
  done
  [ "$tot" -eq "$last" ] && stall=$((stall+EVERY_MIN)) || stall=0; last=$tot
  N=$(printf '%s' "$TAGS" | wc -w)
  printf '%s  %4d/%d pts  %d/%d done |%s | %sG here, %sG usb' "$(date '+%F %T')" \
     "$tot" "$(( 257 * N ))" "$fin" "$N" "$line" "$(fg "$ROOT/runs")" "$(fg "$U")"
  [ "$stall" -ge 45 ] && printf '  ** no progress %d min **' "$stall"; printf '\n'
  [ "$fin" -ge "$N" ] && { echo "SCATTER TEST COMPLETE"; break; }
  sleep $((EVERY_MIN*60))
done
