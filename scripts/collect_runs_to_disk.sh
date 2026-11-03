#!/bin/bash
# Move a finished search's output to the archive disk and free the space here.
#
# Per-run, because a 257-point run writes ~12.3 GB and /data/hz420 does not have room for a whole
# six-run campaign at once. Order is strict: copy, verify, only then delete. Nothing is deleted
# unless --delete is given AND every check passed.
#
#   sh collect_runs_to_disk.sh test2_t7_sf_12p5pct                  # copy + verify, delete nothing
#   sh collect_runs_to_disk.sh --delete test2_t7_sf_12p5pct          # ... then free the trajectories
#   sh collect_runs_to_disk.sh --delete --all-npz <run>             # ... and the per-point npz too
#   DEST=/elsewhere sh collect_runs_to_disk.sh <run>
#
# --delete removes only traj_NNN.npy, which is 98% of the bytes (12.1 of 12.3 GB). The small
# fourier_NNN.npz stay so that re-scoring and the emulator still work here without the archive disk
# mounted; --all-npz removes those as well. results/<run>/ is always copied and never deleted: it is
# ~1 MB and holds evals.csv, run_meta.json and the emulator output.
#
# A run counts as finished when it has as many fourier_NNN.npz as the design has points and an
# evals.csv beside it. Unfinished runs are refused, so this cannot eat a search in progress.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
DEST=${DEST:-${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}}
DESIGN=$ROOT/results/test1_nbody_qmc_focused/design.csv
RSYNC="rsync -rt --no-perms --no-owner --no-group --modify-window=2"

DEL=0; ALLNPZ=0; QUEUE=""
for a in "$@"; do
  case "$a" in
    --delete)   DEL=1 ;;
    --all-npz)  ALLNPZ=1 ;;
    -*)         echo "unknown option: $a" >&2; exit 2 ;;
    *)          QUEUE="$QUEUE $a" ;;
  esac
done
[ -z "$QUEUE" ] && { echo "usage: sh collect_runs_to_disk.sh [--delete] [--all-npz] <run> [<run> ...]" >&2; exit 2; }
[ -d "$DEST" ] || { echo "ERROR: archive disk not mounted at $DEST" >&2; exit 1; }
NPTS=$(( $(wc -l < "$DESIGN") - 1 ))
echo "archive  $DEST"
echo "design   $NPTS points"

for RUN in $QUEUE; do
  SRC="$RUNS_DIR/$RUN"
  echo ""
  echo "########## $RUN ##########"
  [ -d "$SRC" ] || { echo "  ERROR: no such run directory: $SRC"; continue; }

  n_npz=$(find "$SRC" -name 'fourier_*.npz' | wc -l)
  if [ "$n_npz" -ne "$NPTS" ] || [ ! -f "$ROOT/results/$RUN/evals.csv" ]; then
    echo "  REFUSING: $n_npz/$NPTS points done, evals.csv $([ -f "$ROOT/results/$RUN/evals.csv" ] && echo present || echo missing)"
    echo "  (the search is still running; run this again when it has finished)"
    continue
  fi
  echo "  finished: $n_npz/$NPTS points, $(du -sh "$SRC" | cut -f1)"

  echo "  copying ..."
  $RSYNC "$SRC/" "$DEST/runs/$RUN/" || { echo "  ERROR: rsync of runs/ failed; nothing deleted"; continue; }
  $RSYNC "$ROOT/results/$RUN/" "$DEST/results/$RUN/" || { echo "  ERROR: rsync of results/ failed; nothing deleted"; continue; }

  # verify: counts, apparent bytes, a re-scan, and content hashes on the biggest and a random sample
  ns=$(find "$SRC" -type f | wc -l); nd=$(find "$DEST/runs/$RUN" -type f | wc -l)
  bs=$(find "$SRC" -type f -printf '%s\n' | awk '{s+=$1} END {print s}')
  bd=$(find "$DEST/runs/$RUN" -type f -printf '%s\n' | awk '{s+=$1} END {print s}')
  echo "  counts   src $ns files / $bs bytes    dst $nd files / $bd bytes"
  [ "$ns" = "$nd" ] && [ "$bs" = "$bd" ] || { echo "  ERROR: count or size mismatch; nothing deleted"; continue; }
  left=$($RSYNC --itemize-changes --dry-run "$SRC/" "$DEST/runs/$RUN/" | wc -l)
  [ "$left" -eq 0 ] || { echo "  ERROR: rsync re-scan still lists $left items; nothing deleted"; continue; }

  bad=0
  for f in $( (find "$SRC" -type f -size +1M -printf '%f\n' | sort | awk 'NR%37==1';
               find "$SRC" -type f -printf '%s %f\n' | sort -rn | head -2 | cut -d' ' -f2) | sort -u); do
    a=$(md5sum "$SRC/$f" | cut -d' ' -f1); b=$(md5sum "$DEST/runs/$RUN/$f" 2>/dev/null | cut -d' ' -f1)
    [ "$a" = "$b" ] || { echo "  HASH MISMATCH on $f"; bad=1; }
  done
  [ "$bad" -eq 0 ] || { echo "  ERROR: hash mismatch; nothing deleted"; continue; }
  echo "  verified: counts, bytes, re-scan and sampled hashes all match"

  if [ "$DEL" -eq 1 ]; then
    before=$(df -h "$SRC" | awk 'END{print $4}')
    find "$SRC" -name 'traj_*.npy' -delete
    [ "$ALLNPZ" -eq 1 ] && find "$SRC" -name 'fourier_*.npz' -delete
    echo "  freed here: $(du -sh "$SRC" | cut -f1) left in $SRC (was $before free, now $(df -h "$SRC" | awk 'END{print $4}'))"
  else
    echo "  kept here (pass --delete to free the trajectories)"
  fi
done
echo ""
echo "=== done; $(df -h "$RUNS_DIR" | awk 'END{print $4}') free on the working disk ==="
