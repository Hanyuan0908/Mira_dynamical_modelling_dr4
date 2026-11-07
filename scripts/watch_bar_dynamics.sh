#!/bin/bash
# Progress line for the dynamical-bar-length re-integration.
#   per galaxy:  <draws>/N_DRAW  T = truth done
set -u
cd "$(dirname "$0")"; ROOT=$(cd .. && pwd)
OUT=$ROOT/results/bar_dynamics
TAGS=${TAGS:-"t2 t2.375 t2.625 t3 t3.375 t3.625 t4 t4.5 t5 t5.5 t6.25 t7 au5 au9 au10 au12 au14 au17 au18 au20 au22 au24 au26 au27 au28"}
N_DRAW=${N_DRAW:-30}; EVERY_MIN=${EVERY_MIN:-5}
gal(){ sed -n 's/.*"galaxy": "\([^"]*\)".*/\1/p' "$ROOT/results/test2_$1_sf_12p5pct/run_meta.json"; }
last=-1; stall=0
while :; do
  tot=0; fin=0; nt=0; line=""
  for T in $TAGS; do
    G=$(gal "$T"); n=$(ls "$OUT/test2_${T}_sf_12p5pct"/draw_*.npz 2>/dev/null | wc -l)
    tot=$(( tot + n )); t=""
    [ -f "$OUT/truth_${G}.npz" ] && { t="T"; nt=$(( nt + 1 )); }
    [ "$n" -ge "$N_DRAW" ] && [ -n "$t" ] && fin=$(( fin + 1 ))
    line="$line ${T}:${n}${t}"
  done
  [ "$tot" -eq "$last" ] && stall=$(( stall + EVERY_MIN )) || stall=0; last=$tot
  N=$(printf '%s' "$TAGS" | wc -w)
  printf '%s  %4d/%d draws  %d/%d truths  %d/%d complete |%s' "$(date '+%F %T')" \
     "$tot" "$(( N * N_DRAW ))" "$nt" "$N" "$fin" "$N" "$line"
  [ "$stall" -ge 30 ] && printf '  ** no progress %d min **' "$stall"; printf '\n'
  [ "$fin" -ge "$N" ] && { echo "BAR DYNAMICS COMPLETE"; break; }
  sleep $(( EVERY_MIN * 60 ))
done
