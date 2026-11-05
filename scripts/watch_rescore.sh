#!/bin/bash
# Periodic progress line for the sfratio re-score. Read-only.
set -u
cd "$(dirname "$0")"; ROOT=$(cd .. && pwd)
HALOS=${HALOS:-"9 10 17 18 24 26 27"}; EVERY_MIN=${EVERY_MIN:-5}
while :; do
  tot=0; done_g=0; line=""
  for H in $HALOS; do
    R=test2_au${H}_sf_12p5pct
    n=$(find "$ROOT/runs/$R" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
    e=$(ls "$ROOT/results/$R"/emulator_results_*sfratio.json 2>/dev/null | wc -l)
    tot=$(( tot + n )); [ "$e" -ge 1 ] && done_g=$(( done_g + 1 ))
    line="$line au${H}:${n}$([ "$e" -ge 1 ] && echo E)"
  done
  echo "$(date '+%F %T')  $tot/$(( 257 * $(printf '%s' "$HALOS" | wc -w) )) pts  ${done_g}/$(printf '%s' "$HALOS" | wc -w) with posteriors |$line"
  [ "$done_g" -eq "$(printf '%s' "$HALOS" | wc -w)" ] && { echo "ALL DONE"; break; }
  sleep $(( EVERY_MIN * 60 ))
done
