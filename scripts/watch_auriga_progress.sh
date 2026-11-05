#!/bin/bash
# Progress log for the Auriga campaign, for reading the morning after.
#
#   sh watch_auriga_progress.sh                      # log every 5 min until all galaxies are done
#   HALOS="9 10" EVERY_MIN=2 sh watch_auriga_progress.sh
#   sh watch_auriga_progress.sh --once               # one status block, then exit
#
# Reports, per galaxy: design points finished (of 257), whether the search has merged to an
# evals.csv, whether the GP posterior landed, and whether it has been archived to the portable
# disk. Plus free space on both disks, and a stall warning if nothing has moved for a while --
# which is the signal that the cluster job died, or that the archiver has stopped keeping up and
# the compute script is blocked waiting for space.
#
# Read-only: it writes nothing but its own log.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
DEST=${DEST:-${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}}
TAG=${TAG:-12p5pct}
HALOS=${HALOS:-"9 10 17 18 24 26 27"}
EVERY_MIN=${EVERY_MIN:-5}
NPTS=${NPTS:-257}
STALL_MIN=${STALL_MIN:-45}

ONCE=0
for a in "$@"; do case "$a" in --once) ONCE=1 ;; *) echo "unknown option: $a" >&2; exit 2 ;; esac; done

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

last_total=-1
stalled=0
while : ; do
  total=0; done_gals=0; line=""
  for H in $HALOS; do
    RUN="test2_au${H}_sf_${TAG}"
    SRC="$RUNS_DIR/$RUN"
    n=$(find "$SRC" -name 'fourier_*.npz' 2>/dev/null | wc -l)
    total=$(( total + n ))
    ev=$([ -f "$ROOT/results/$RUN/evals.csv" ] && echo m || echo .)
    em=$(ls "$ROOT/results/$RUN"/emulator_results_*.json 2>/dev/null | head -1)
    em=$([ -n "$em" ] && echo E || echo .)
    ar=$([ -d "$DEST/runs/$RUN" ] && echo A || echo .)
    [ "$n" -ge "$NPTS" ] && done_gals=$(( done_gals + 1 ))
    line="$line au${H}:${n}${ev}${em}${ar}"
  done

  if [ "$total" -eq "$last_total" ]; then
    stalled=$(( stalled + EVERY_MIN ))
  else
    stalled=0
  fi
  last_total=$total

  printf '%s  %4d/%d pts  %d/%d galaxies  |%s  | %sG here, %sG usb' \
    "$(date '+%F %T')" "$total" "$(( NPTS * $(printf '%s' "$HALOS" | wc -w) ))" \
    "$done_gals" "$(printf '%s' "$HALOS" | wc -w)" "$line" \
    "$(free_gb "$RUNS_DIR")" "$(free_gb "$DEST")"
  if [ "$stalled" -ge "$STALL_MIN" ]; then
    printf '   ** NO PROGRESS FOR %d MIN **' "$stalled"
  fi
  printf '\n'

  [ "$ONCE" = 1 ] && break
  if [ "$done_gals" -eq "$(printf '%s' "$HALOS" | wc -w)" ]; then
    echo "$(date '+%F %T')  all $done_gals galaxies have finished their 257 points."
    echo "key:  <points><m=merged><E=emulator done><A=archived to usb>"
    break
  fi
  sleep $(( EVERY_MIN * 60 ))
done
