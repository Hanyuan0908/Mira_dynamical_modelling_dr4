#!/bin/bash
# Workstation side of the per-star distance-error test: re-score, build the posteriors, archive.
#
#   nohup sh postprocess_scatter_test.sh > ../logs/scatter_post.log 2>&1 &
#   sh postprocess_scatter_test.sh --once
#
# Identical treatment to the fiducial campaign, which is the point: the data has a per-star
# distance error but everything downstream assumes the nominal uniform 12.5%, so any shift in the
# recovered parameters against the fiducial run is caused by that mismatch alone.
#   re-score   sfratio (rho_obs/rho_lib, nside 32, 1 kpc) + K = 10
#   posterior  scheme B, logLs_..._vR_sigR_sigphi_m24_vphi_m024_sfratio
#   archive    collect_runs_to_disk.sh --delete
#
# Cadence: these runs are integrated at 50 Myr, so the re-score takes --cadence-step 1. The script
# reads the snapshot count out of a trajectory rather than assuming, because the rescore thins
# UNCONDITIONALLY and passing 2 to an already-thinned file would silently halve it again.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
DEST=${DEST:-${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}}
WORKERS=${WORKERS:-12}
WAIT_MIN=${WAIT_MIN:-10}
TAGS=${TAGS:-"t3 t5.5 t7 au18 au26"}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
N_FULL=41

ONCE=0
for a in "$@"; do case "$a" in --once) ONCE=1 ;; *) echo "unknown option: $a" >&2; exit 2 ;; esac; done
[ -d "$DEST" ] || { echo "ERROR: portable disk not mounted at $DEST" >&2; exit 1; }
free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

while : ; do
  pending=0
  for T in $TAGS; do
    RUN=test2_${T}_sf_12p5pct_scat
    SRC=$RUNS_DIR/$RUN
    [ -f "$ROOT/results/$RUN/emulator_results_logLs_${B}_sfratio.json" ] && [ -d "$DEST/runs/$RUN" ] && continue
    n=$(find "$SRC" -name 'fourier_*.npz' 2>/dev/null | wc -l)
    if [ "$n" -lt 257 ] || [ ! -f "$ROOT/results/$RUN/evals.csv" ]; then
      pending=$(( pending + 1 )); continue
    fi
    echo ""
    echo "########## $RUN  $(date '+%F %T') ##########"
    nsnap=$("$PY" -c "
import numpy as np, glob
f=sorted(glob.glob('$SRC/traj_*.npy'))
print(np.load(f[0], mmap_mode='r').shape[1] if f else 0)")
    [ "$nsnap" = "0" ] && { echo "    no trajectories; cannot re-score"; continue; }
    if [ "$nsnap" -ge "$N_FULL" ]; then STEP=2; else STEP=1; fi
    echo "    $nsnap snapshots -> --cadence-step $STEP"

    nsf=$(find "$SRC" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
    if [ "$nsf" -lt 257 ]; then
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
    [ "$nsf" -lt 257 ] && { echo "    incomplete ($(find "$SRC" -name 'sfratio_*.lock'|wc -l) stale locks); not archiving"; continue; }

    OMP_NUM_THREADS=4 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" \
      --column "logLs_${B}_sfratio" --loo-top 32 > "/tmp/${RUN}_emuB.log" 2>&1 \
      && echo "    posterior ok" || { echo "    posterior FAILED"; tail -3 "/tmp/${RUN}_emuB.log"; continue; }
    OMP_NUM_THREADS=4 "$PY" -W ignore plot_test1_logl_corner.py --run "$RUN" \
      --column "logLs_${B}_sfratio" --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2 \
      > "/tmp/${RUN}_cornerB.log" 2>&1 || true
    sh ./collect_runs_to_disk.sh --delete "$RUN" || echo "    WARNING: archiving failed"
    echo "    $(free_gb "$RUNS_DIR") GB free"
  done
  left=0
  for T in $TAGS; do
    [ -f "$ROOT/results/test2_${T}_sf_12p5pct_scat/emulator_results_logLs_${B}_sfratio.json" ] || left=$(( left + 1 ))
  done
  N=$(printf '%s' "$TAGS" | wc -w)
  echo "$(date '+%F %T')  $(( N - left ))/$N complete; $(free_gb "$RUNS_DIR") GB free"
  [ "$ONCE" = 1 ] && break
  [ "$left" -eq 0 ] && { echo "scatter test complete."; break; }
  sleep $(( WAIT_MIN * 60 ))
done
