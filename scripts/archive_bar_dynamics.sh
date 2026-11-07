#!/bin/bash
# Mirror the dynamical-bar-length output to the portable disk. Copy, verify, then (optionally) delete.
#
#   sh archive_bar_dynamics.sh                 # copy + verify once, delete nothing
#   sh archive_bar_dynamics.sh --watch         # keep mirroring every WAIT_MIN while the job runs
#   sh archive_bar_dynamics.sh --delete        # copy + verify, then remove the local copy
#
# Order is strict, as in collect_runs_to_disk.sh: nothing is deleted until the copy has been
# verified by count AND by content hash on a sample.
#
# NOTE ON SIZE. This campaign writes no trajectories, so it is ~600 MB in total (25 truths at ~3.4 MB
# and 750 draws at ~0.7 MB), not the ~12.6 GB per run the orbit-library campaigns produced. It is
# also the ANALYSIS INPUT, not an intermediate - every dynamical-length plot reads it - so --delete
# makes the analysis depend on the disk being mounted. Mirror by default; delete only on request.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
SRC=$ROOT/results/bar_dynamics
DEST=${DEST:-${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}}
WAIT_MIN=${WAIT_MIN:-10}
N_DRAW=${N_DRAW:-30}
N_GAL=${N_GAL:-25}
RSYNC="rsync -rt --no-perms --no-owner --no-group --modify-window=2"

DEL=0; WATCH=0
for a in "$@"; do case "$a" in
  --delete) DEL=1 ;; --watch) WATCH=1 ;;
  *) echo "unknown option: $a" >&2; exit 2 ;; esac; done
[ -d "$DEST" ] || { echo "ERROR: portable disk not mounted at $DEST" >&2; exit 1; }
mkdir -p "$SRC" "$DEST/results/bar_dynamics"

once () {
  nt=$(ls "$SRC"/truth_*.npz 2>/dev/null | wc -l)
  nd=$(find "$SRC" -name 'draw_*.npz' 2>/dev/null | wc -l)
  [ "$nt" -eq 0 ] && [ "$nd" -eq 0 ] && { echo "$(date '+%F %T')  nothing to archive yet"; return 1; }

  $RSYNC "$SRC/" "$DEST/results/bar_dynamics/" || { echo "  ERROR: rsync failed; nothing deleted"; return 1; }

  # verify: the destination must hold at least what the source does, and a sample must hash equal
  mt=$(ls "$DEST/results/bar_dynamics"/truth_*.npz 2>/dev/null | wc -l)
  md=$(find "$DEST/results/bar_dynamics" -name 'draw_*.npz' 2>/dev/null | wc -l)
  if [ "$mt" -lt "$nt" ] || [ "$md" -lt "$nd" ]; then
    echo "  ERROR: destination has $mt truths / $md draws against $nt / $nd here; nothing deleted"
    return 1
  fi
  bad=0; checked=0
  for f in $(find "$SRC" -name '*.npz' | sort | sed -n '1p;$p' ; find "$SRC" -name '*.npz' | shuf -n 3 2>/dev/null); do
    rel=${f#"$SRC"/}; checked=$(( checked + 1 ))
    cmp -s "$f" "$DEST/results/bar_dynamics/$rel" || { echo "  MISMATCH: $rel"; bad=$(( bad + 1 )); }
  done
  [ "$bad" -eq 0 ] || { echo "  ERROR: $bad/$checked sampled files differ; nothing deleted"; return 1; }
  echo "$(date '+%F %T')  mirrored $nt truths + $nd draws ($(du -sh "$SRC" | cut -f1)); "\
"verified $checked files byte-identical; usb $(df -h "$DEST" | awk 'END{print $4}') free"

  if [ "$DEL" = 1 ]; then
    if [ "$nt" -ge "$N_GAL" ] && [ "$nd" -ge "$(( N_GAL * N_DRAW ))" ]; then
      rm -rf "$SRC"/truth_*.npz "$SRC"/test2_*
      echo "  deleted the local copy; it now lives only on $DEST"
    else
      echo "  --delete ignored: the campaign is not complete ($nt/$N_GAL truths, $nd/$(( N_GAL * N_DRAW )) draws)"
    fi
  fi
  return 0
}

if [ "$WATCH" = 1 ]; then
  while :; do
    once
    nt=$(ls "$SRC"/truth_*.npz 2>/dev/null | wc -l)
    nd=$(find "$SRC" -name 'draw_*.npz' 2>/dev/null | wc -l)
    [ "$nt" -ge "$N_GAL" ] && [ "$nd" -ge "$(( N_GAL * N_DRAW ))" ] && { echo "campaign complete and mirrored."; break; }
    sleep $(( WAIT_MIN * 60 ))
  done
else
  once
fi
