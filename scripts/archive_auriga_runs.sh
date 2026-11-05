#!/bin/bash
# Drain finished Auriga runs to the portable disk while the campaign is still going.
# Run this ON THE WORKSTATION (the only machine with the USB disk mounted):
#
#   sh archive_auriga_runs.sh --delete          # watch, archive and free space until all 13 are done
#   sh archive_auriga_runs.sh                   # copy + verify only, delete nothing (dry-ish run)
#   sh archive_auriga_runs.sh --delete --once   # one pass over whatever is finished, then exit
#
# Why this exists: a 257-point run writes ~12.3 GB and /data/hz420 has ~55 GB free, so the 13-galaxy
# campaign needs ~160 GB it does not have. The compute script blocks when the disk fills; this moves
# finished runs out from under it so it can carry on. Start it before or alongside the campaign.
#
# All the real work is collect_runs_to_disk.sh, which copies, verifies (file counts, byte totals, an
# rsync re-scan and md5 on the largest and a sample of files) and only then deletes the trajectories.
# A run is touched only when it has all 257 fourier_*.npz and an evals.csv, so an in-progress search
# is never eaten. --delete removes traj_*.npy (98% of the bytes) and keeps the small fourier_*.npz
# here, so re-scoring and the emulator still work without the USB disk mounted.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
DEST=${DEST:-${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}}
TAG=${TAG:-12p5pct}
HALOS=${HALOS:-"5 9 10 12 14 17 18 20 22 24 26 27 28"}
WAIT_MIN=${WAIT_MIN:-15}

DEL=""; ONCE=0
for a in "$@"; do
  case "$a" in
    --delete) DEL="--delete" ;;
    --once)   ONCE=1 ;;
    *)        echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

[ -d "$DEST" ] || { echo "ERROR: portable disk not mounted at $DEST" >&2; exit 1; }
[ -n "$DEL" ] || echo "NOTE: no --delete, so nothing will be freed here (copy + verify only)."

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

while : ; do
  pending=0; archived=0
  for H in $HALOS; do
    RUN="test2_au${H}_sf_${TAG}"
    SRC="$RUNS_DIR/$RUN"
    # already drained: results kept, trajectories gone
    if [ -d "$DEST/runs/$RUN" ] && [ -z "$(find "$SRC" -name 'traj_*.npy' 2>/dev/null | head -1)" ] \
       && [ -d "$SRC" ]; then
      continue
    fi
    [ -d "$SRC" ] || { pending=$(( pending + 1 )); continue; }
    n_npz=$(find "$SRC" -name 'fourier_*.npz' 2>/dev/null | wc -l)
    if [ "$n_npz" -ne 257 ] || [ ! -f "$ROOT/results/$RUN/evals.csv" ]; then
      pending=$(( pending + 1 )); continue
    fi
    echo ""
    echo "=== $(date '+%F %T')  archiving $RUN ==="
    sh ./collect_runs_to_disk.sh $DEL "$RUN" && archived=$(( archived + 1 ))
  done

  echo ""
  echo "$(date '+%F %T'): $archived archived this pass, $pending still running or not started; "\
"$(free_gb "$RUNS_DIR") GB free here, $(free_gb "$DEST") GB free on the portable disk"
  [ "$ONCE" = 1 ] && break
  [ "$pending" -eq 0 ] && { echo "all 13 runs archived."; break; }
  sleep $(( WAIT_MIN * 60 ))
done
