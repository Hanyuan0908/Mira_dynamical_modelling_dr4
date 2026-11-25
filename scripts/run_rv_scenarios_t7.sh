#!/bin/bash
# Gaia-RV sub-samples of the t = 7 mock (12.5%), re-scored from the saved trajectories (no integration):
# masks from make_rv_selection_masks.py; adopted configuration (sfratio, K = 10, 50 Myr, scheme B).
#   caffeinate -i -s nohup bash run_rv_scenarios_t7.sh > ~/mira_dr4_runs/rv_scenarios_t7.log 2>&1 &
# Trajectories live on the USB copy of the cluster runs (MIRA_DR4_RUNS). Control (2026-11-21): the all-stars
# mask reproduces the run's plain sfratio column exactly at three design points.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-$HOME/miniforge3/envs/astro312/bin/python}
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs}
RUN=test2_t7_sf_12p5pct
STEM=nbody_mock_selection_sigmaDM0.271_N50k
WORKERS=${WORKERS:-5}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
for TAG in random_osargN osarg_cover g16 g15 g14 grvs14; do
  MASK=results/rv_selection_masks/${STEM}__${TAG}.npy
  LOGS=$MIRA_DR4_RUNS/$RUN/logs_rv; mkdir -p "$LOGS"
  echo "=== $TAG: re-score $(date '+%T') ==="
  pids=""; w=0
  while [ "$w" -lt "$WORKERS" ]; do
    OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 \
      --cadence-step 2 --keep "$MASK" --keep-tag "$TAG" --stride "$WORKERS" --offset "$w" \
      > "$LOGS/${TAG}_w$w.log" 2>&1 &
    pids="$pids $!"; w=$(( w + 1 ))
  done
  for p in $pids; do wait "$p"; done
  n=$(ls "$MIRA_DR4_RUNS/$RUN" 2>/dev/null | grep -E "^sfratio_${TAG}_[0-9]+\.npz$" | wc -l | tr -d ' ')
  echo "$TAG: $n/257 re-scored"
  if [ "$n" -ne 257 ]; then echo "FAILED $TAG"; grep -h -E "Error|Trace" "$LOGS/${TAG}"_w*.log | head -5; exit 1; fi
  "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 2 \
    --keep "$MASK" --keep-tag "$TAG" --merge | tail -1
  OMP_NUM_THREADS=8 "$PY" -W ignore test1_gp_emulator.py --run "${RUN}_${TAG}" --column "logLs_${B}_sfratio_${TAG}" \
    --loo-top 32 > "$LOGS/${TAG}_emulator.log" 2>&1 || { echo "EMULATOR FAILED $TAG"; tail -3 "$LOGS/${TAG}_emulator.log"; exit 1; }
  echo "=== $TAG done $(date '+%T') ==="
done
echo "ALL OK"
