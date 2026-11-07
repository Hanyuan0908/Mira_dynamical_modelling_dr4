#!/bin/bash
# The whole thing in one pass, locally: t=7 mock Miras at 12.5% distance error, 257 Sobol points,
# trajectories at 50 Myr AND the five bar-frame orbital frequencies per star, then the usual
# re-score + posterior. Nothing new here - it is the existing search with --freq-snap 1001 added.
#
#   nohup sh run_t7_freq_local.sh > ../logs/t7_freq.log 2>&1 &
#
# The dense sampling is bitwise identical to trajsize 41 before thinning (verified: control with
# --freq-snap off agrees to 8e-15), so --freq-snap moves nothing. It does NOT follow that this run
# reproduces test2_t7_sf_12p5pct column for column: that one was SEARCHED at 25 Myr (its saved
# trajectory is (50218, 41, 6)) and only thinned to 50 Myr at re-score, so its search-time scheme-A
# columns use 30 library snapshots against 15 here - a real 7.8% difference that is the cadence and
# nothing else. The column to compare is the re-scored sfratio one, 21 snapshots at the same times
# in both. NOTE: run_meta.json for that run claims cadence_step 2 / 50 Myr and is WRONG - it is
# rewritten by later invocations and does not describe the integration (PROGRESS.md, 2026-11-05).
# Measured: 148 s integrate + 2 s score, so 257 points on 2 workers is ~5.5 h and ~7.2 GB.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
ROOT=$(cd .. && pwd)
RUN=${RUN:-test2_t7_sf_12p5pct_freq}
MOCK=nbody_mock_selection_sigmaDM0.271_N50k.hdf5
THREADS=${THREADS:-8}; WORKERS=${WORKERS:-2}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024

echo "=== $RUN: 257 points, 50 Myr trajectories + frequencies at 1 Myr  $(date '+%F %T') ==="
pids=""; w=0
while [ "$w" -lt "$WORKERS" ]; do
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
    --mock "$MOCK" --run "$RUN" --galaxy nbody --cadence-step 2 --freq-snap 1001 \
    --stride "$WORKERS" --offset "$w" > "/tmp/${RUN}_w${w}.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
  --mock "$MOCK" --run "$RUN" --galaxy nbody --cadence-step 2 --merge
echo "=== search done $(date '+%F %T'); $(find "$ROOT/runs/$RUN" -name 'freq_*.npz' | wc -l) frequency tables ==="

# re-score into the adopted configuration. The trajectories are ALREADY 50 Myr, so step 1.
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
echo "=== rescored $(find "$ROOT/runs/$RUN" -name 'sfratio_*.npz' | wc -l)/257  $(date '+%F %T') ==="

OMP_NUM_THREADS=4 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" \
  --column "logLs_${B}_sfratio" --loo-top 32
OMP_NUM_THREADS=4 "$PY" -W ignore plot_corner_with_kinematic_prior.py --run "$RUN" \
  --column "logLs_${B}_sfratio" --measure-from-mock || true
echo "=== ALL DONE $(date '+%F %T') ==="
