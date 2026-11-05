#!/bin/bash
# Workstation side of Auriga batch 2: watch for each galaxy the cluster finishes, then
#
#   1. re-score it into the adopted configuration   sfratio (rho_obs/rho_lib) + K = 10 + 50 Myr
#   2. build the posterior on chi2 SCHEME B          logLs_..._vR_sigR_sigphi_m24_vphi_m024_sfratio
#   3. copy run + results to the portable disk, verify, and free the trajectories here
#
# Run it here before or alongside the cluster job:
#
#   nohup sh postprocess_auriga_batch2.sh > ../logs/batch2_post.log 2>&1 &
#   HALOS="5 12" sh postprocess_auriga_batch2.sh --once     # one pass over what is already done
#
# Ordering matters and is the one improvement over batch 1: the re-score needs the trajectories, so
# it runs BEFORE the archive step. Batch 1 was archived first and every galaxy then had to be pulled
# back off the USB, 13 GB at a time.
#
# Cadence: the rescore thins by --cadence-step UNCONDITIONALLY -- there is no guard against
# double-thinning. Batch 2 is integrated at 50 Myr (21 of 41 snapshots), so it must be re-scored
# with --cadence-step 1, where batch 1 (saved at 25 Myr) needed 2. This script reads the snapshot
# count out of a trajectory file and picks the right value rather than trusting either.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
DEST=${DEST:-${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}}
HALOS=${HALOS:-"5 12 14 20 22 28"}
TAG=${TAG:-12p5pct}
WORKERS=${WORKERS:-12}
WAIT_MIN=${WAIT_MIN:-10}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
N_FULL=41                       # ssl.N_SNAP: snapshots over 1 Gyr at 25 Myr

ONCE=0
for a in "$@"; do case "$a" in --once) ONCE=1 ;; *) echo "unknown option: $a" >&2; exit 2 ;; esac; done
[ -d "$DEST" ] || { echo "ERROR: portable disk not mounted at $DEST" >&2; exit 1; }

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

done_file () { [ -d "$DEST/runs/$1" ] && [ -f "$ROOT/results/$1/emulator_results_logLs_${B}_sfratio.json" ]; }

while : ; do
  pending=0
  for H in $HALOS; do
    RUN=test2_au${H}_sf_${TAG}
    SRC=$RUNS_DIR/$RUN
    done_file "$RUN" && continue
    n=$(find "$SRC" -name 'fourier_*.npz' 2>/dev/null | wc -l)
    if [ "$n" -lt 257 ] || [ ! -f "$ROOT/results/$RUN/evals.csv" ]; then
      pending=$(( pending + 1 )); continue
    fi

    echo ""
    echo "########## $RUN  $(date '+%F %T') ##########"
    echo "    search finished ($n/257); $(free_gb "$RUNS_DIR") GB free"

    # --- 1. re-score: sfratio + K=10, at whatever cadence step this file still needs -------------
    nsnap=$("$PY" -c "
import numpy as np, glob
f=sorted(glob.glob('$SRC/traj_*.npy'))
print(np.load(f[0], mmap_mode='r').shape[1] if f else 0)")
    if [ "$nsnap" -ge "$N_FULL" ]; then STEP=2; else STEP=1; fi
    echo "    trajectory has $nsnap snapshots -> --cadence-step $STEP"
    if [ "$nsnap" = "0" ]; then
      echo "    ERROR: no trajectories (already archived?); cannot re-score. Skipping."; continue
    fi

    nsf=$(find "$SRC" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
    if [ "$nsf" -lt 257 ]; then
      echo "    re-scoring sfratio + K=10 ($nsf/257 done) ..."
      pids=""; w=0
      while [ "$w" -lt "$WORKERS" ]; do
        OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" \
          --estimator ratio --n-draw 10 --cadence-step "$STEP" --stride "$WORKERS" --offset "$w" \
          > "/tmp/${RUN}_rs${w}.log" 2>&1 &
        pids="$pids $!"; w=$(( w + 1 ))
      done
      for p in $pids; do wait "$p"; done
    fi
    OMP_NUM_THREADS=2 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" \
      --estimator ratio --n-draw 10 --cadence-step "$STEP" --merge > "/tmp/${RUN}_merge.log" 2>&1
    nsf=$(find "$SRC" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
    echo "    scored $nsf/257"
    if [ "$nsf" -lt 257 ]; then
      # a killed worker leaves an O_EXCL lock and that point is dropped to NaN for ever
      echo "    WARNING: $nsf/257 only. Stale locks: $(find "$SRC" -name 'sfratio_*.lock' | wc -l)."
      echo "    Not archiving; clear the locks and re-run this script."
      continue
    fi

    # --- 2. posterior on scheme B -----------------------------------------------------------------
    OMP_NUM_THREADS=4 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" \
      --column "logLs_${B}_sfratio" --loo-top 32 > "/tmp/${RUN}_emuB.log" 2>&1 \
      && echo "    posterior ok (scheme B / sfratio)" \
      || { echo "    posterior FAILED (see /tmp/${RUN}_emuB.log)"; tail -3 "/tmp/${RUN}_emuB.log"; continue; }
    OMP_NUM_THREADS=4 "$PY" -W ignore plot_test1_logl_corner.py --run "$RUN" \
      --column "logLs_${B}_sfratio" --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2 \
      > "/tmp/${RUN}_cornerB.log" 2>&1 || true

    # --- 3. archive and free, after the re-score so the sfratio files go too ---------------------
    sh ./collect_runs_to_disk.sh --delete "$RUN" || echo "    WARNING: archiving failed; kept here"
    echo "    $(free_gb "$RUNS_DIR") GB free after $RUN"
  done

  left=0
  for H in $HALOS; do done_file "test2_au${H}_sf_${TAG}" || left=$(( left + 1 )); done
  echo "$(date '+%F %T')  $(( $(printf '%s' "$HALOS" | wc -w) - left ))/$(printf '%s' "$HALOS" | wc -w) complete; $left to go; $(free_gb "$RUNS_DIR") GB free here, $(free_gb "$DEST") GB on the portable disk"
  [ "$ONCE" = 1 ] && break
  [ "$left" -eq 0 ] && { echo "batch 2 complete."; break; }
  sleep $(( WAIT_MIN * 60 ))
done
