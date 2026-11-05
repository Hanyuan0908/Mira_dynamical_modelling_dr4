#!/bin/bash
# Re-score the Auriga campaign with the configuration the record actually calls for.
#
#   sh rescore_auriga_sfratio.sh            # all 7 of batch 1
#   HALOS="9 10" WORKERS=8 sh rescore_auriga_sfratio.sh
#
# The campaign was scored with `sf` (the true S tables from the mock), K = 1, and a 25 Myr
# cadence. All three are wrong for this project:
#
#   sfratio   S estimated as rho_obs/rho_lib per (HEALPix nside 32 x 1 kpc) cell, which is what a
#             real survey can do. PROGRESS 2026-11-02/03: at 12.5% it gives 3.1% mean |Omega_b
#             bias| against 4.1% for the true S and 9.3% for no weighting -- and the true S carries
#             a POSITIVE Omega_b bias (t=3 +6.0%, t=4 +9.5%, t=7 +2.6%) that the ratio removes
#             (+2.5%, +2.3%, -2.0%). The uniformly positive +7.7% seen across the 7 Auriga galaxies
#             on the `sf` column is that artefact.
#   K = 10    ten pooled distance-error draws per library sample. A one-sample library biases chi2
#             high by a parameter-dependent amount that averaging separate chi2 values cannot fix;
#             K = 10 cuts the emulator leave-one-out error 1.47 -> 0.65 (selection.py).
#   50 Myr    --cadence-step 2 thins the saved 25 Myr trajectory while scoring, putting these runs
#             on the same grid as test2_t7_sf_12p5pct and the rest of the recent N-body set.
#
# All of this is a re-score: no orbit integration. It reads the saved trajectories, which live on
# the portable disk, so each galaxy is restored from there, scored, and its trajectories dropped
# again before the next one -- /data has room for about three at 13 GB each.
#
# Outputs per galaxy: runs/<run>/sfratio_NNN.npz, results/<run>/evals_sfratio.csv, the sfratio
# columns folded into results/<run>/evals.csv, and the GP posterior + corner plot for
# logLs_<scheme>_sfratio.
set -u
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
U=${MIRA_DR4_ARCHIVE:-/run/media/hz420/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4}
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
WORKERS=${WORKERS:-12}
HALOS=${HALOS:-"9 10 17 18 24 26 27"}
TAG=${TAG:-12p5pct}
SCHEME=fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024
NEED_GB=${NEED_GB:-16}
RSYNC="rsync -rt --no-perms --no-owner --no-group --modify-window=2"

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

[ -d "$U" ] || { echo "ERROR: portable disk not mounted at $U" >&2; exit 1; }

# Never run two sets of workers over the same run: a killed worker leaves an O_EXCL lock behind
# and that point is then skipped for good. Wait for anything already scoring to finish.
while pgrep -f 'rescore_empirical_selection.py' > /dev/null 2>&1; do
  echo "$(date '+%H:%M:%S')  waiting for a running rescore to finish ..."
  sleep 60
done

echo "=== sfratio + K=10 + 50 Myr re-score: halos $HALOS, $WORKERS workers ==="
for H in $HALOS; do
  RUN=test2_au${H}_sf_${TAG}
  SRC=$RUNS_DIR/$RUN
  echo ""
  echo "########## $RUN  $(date '+%F %T') ##########"

  ndone=$(find "$SRC" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
  if [ "$ndone" -lt 257 ]; then
    ntraj=$(find "$SRC" -name 'traj_*.npy' 2>/dev/null | wc -l)
    if [ "$ntraj" -lt 257 ]; then
      while [ "$(free_gb "$RUNS_DIR")" -lt "$NEED_GB" ]; do
        echo "    only $(free_gb "$RUNS_DIR") GB free, need $NEED_GB; waiting 5 min"
        sleep 300
      done
      echo "    restoring trajectories from the portable disk ..."
      $RSYNC "$U/runs/$RUN/" "$SRC/" || { echo "    ERROR: restore failed; skipping"; continue; }
    fi
    echo "    scoring ($ndone/257 already done) ..."
    pids=""; w=0
    while [ "$w" -lt "$WORKERS" ]; do
      OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" \
        --estimator ratio --n-draw 10 --cadence-step 2 --stride "$WORKERS" --offset "$w" \
        > "/tmp/${RUN}_rs${w}.log" 2>&1 &
      pids="$pids $!"; w=$(( w + 1 ))
    done
    for p in $pids; do wait "$p"; done
  else
    echo "    already scored (257/257)"
  fi

  OMP_NUM_THREADS=2 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" \
    --estimator ratio --n-draw 10 --cadence-step 2 --merge >> "/tmp/${RUN}_merge.log" 2>&1
  n=$(find "$SRC" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
  echo "    scored $n/257  $(date '+%H:%M:%S')"

  OMP_NUM_THREADS=4 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" \
    --column "logLs_${SCHEME}_sfratio" --loo-top 32 > "/tmp/${RUN}_emu_sfratio.log" 2>&1 \
    && echo "    posterior ok" || echo "    posterior FAILED (see /tmp/${RUN}_emu_sfratio.log)"
  OMP_NUM_THREADS=4 "$PY" -W ignore plot_test1_logl_corner.py --run "$RUN" \
    --column "logLs_${SCHEME}_sfratio" \
    --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2 \
    > "/tmp/${RUN}_corner_sfratio.log" 2>&1 || true

  # trajectories are safe on the portable disk; drop them here to make room for the next galaxy
  find "$SRC" -name 'traj_*.npy' -delete
  echo "    trajectories dropped; $(free_gb "$RUNS_DIR") GB free"
done
echo ""
echo "=== re-score done $(date '+%F %T') ==="
