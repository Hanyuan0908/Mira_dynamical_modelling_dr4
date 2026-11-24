#!/bin/bash
# Fourier-grid strategies for signal/noise (user, 2026-11-24): re-score a run's saved trajectories on coarser
# (ring, |z| layer) grids and refit the emulator on each. No new integration.
#   caffeinate -i -s nohup bash run_grid_tests.sh <run> > ~/mira_dr4_runs/grid_tests_<run>.log 2>&1 &
# Grids (ring width inside 6 kpc / z max / layer width):  G1 0.5/2/1   G2 0.5/2/2   G3 1/2/2   G4 0.5/3/3   G5 0.3/2/2
set -u
cd "$(dirname "$0")"
RUN=$1
PY=$HOME/miniforge3/envs/astro312/bin/python
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-$HOME/mira_dr4_runs}
RS_WORKERS=${RS_WORKERS:-6}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
LOGS=$MIRA_DR4_RUNS/$RUN/logs; mkdir -p "$LOGS"
for g in "0.5 2 1" "0.5 2 2" "1 2 2" "0.5 3 3" "0.3 2 2"; do
  set -- $g
  RW=$1; ZM=$2; LW=$3
  echo "=== grid ring $RW, z max $ZM, layer $LW: re-score $(date '+%T') ==="
  pids=""; w=0
  while [ "$w" -lt "$RS_WORKERS" ]; do
    OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 \
      --cadence-step 1 --ring-width "$RW" --z-max "$ZM" --layer-width "$LW" --stride "$RS_WORKERS" --offset "$w" \
      > "$LOGS/rescore_grid_${RW}_${ZM}_${LW}_w$w.log" 2>&1 &
    pids="$pids $!"; w=$(( w + 1 ))
  done
  for p in $pids; do wait "$p"; done
  "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 1 \
    --ring-width "$RW" --z-max "$ZM" --layer-width "$LW" --merge | tail -1
  TAG=$("$PY" -c "import test1_nbody_qmc_focused_fourier_stars as fs; print(fs.set_grid($RW, $ZM, $LW))")
  ( OMP_NUM_THREADS=2 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" --column "logLs_${B}_sfratio${TAG}" --loo-top 32 \
      > "$LOGS/emulator_grid${TAG}.log" 2>&1 && echo "emulator OK${TAG}" || echo "emulator FAILED${TAG}" ) &
done
wait
echo "ALL DONE $(date '+%T')"
