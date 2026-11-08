#!/bin/bash
# Stages 2-4 of run_t7_freq_local.sh, which died at the handover because the script was EDITED
# WHILE RUNNING: sh reads a script incrementally by byte offset, so changing its length made the
# shell resume mid-token. The search itself completed (257/257, merged). NEVER edit a running
# shell script; copy it, edit the copy, and launch that.
# Trajectories are already 21 snapshots, so the re-score takes --cadence-step 1.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
ROOT=$(cd .. && pwd)
RUN=test2_t7_sf_12p5pct_freq
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
echo "=== re-score (sfratio, K=10, 50 Myr) $(date '+%F %T') ==="
pids=""; w=0
while [ "$w" -lt 12 ]; do
  OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" \
    --estimator ratio --n-draw 10 --cadence-step 1 --stride 12 --offset "$w" \
    > "/tmp/${RUN}_rs${w}.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
OMP_NUM_THREADS=2 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" \
  --estimator ratio --n-draw 10 --cadence-step 1 --merge
echo "=== scored $(ls "$ROOT/runs/$RUN"/sfratio_*.npz 2>/dev/null | wc -l)/257  $(date '+%F %T') ==="
OMP_NUM_THREADS=8 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" \
  --column "logLs_${B}_sfratio" --loo-top 32
OMP_NUM_THREADS=8 "$PY" -W ignore plot_corner_with_kinematic_prior.py --run "$RUN" \
  --column "logLs_${B}_sfratio" --measure-from-mock || true
echo "=== ALL DONE $(date '+%F %T') ==="
