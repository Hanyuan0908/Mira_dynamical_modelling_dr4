#!/bin/bash
# Radial ring schemes refined around the Milky Way's ILR (user, 2026-11-25), all with one |z| < 2 layer: re-score a run's
# saved trajectories, merge (once, after the workers), refit the emulator. No integration.
#   caffeinate -i -s nohup bash run_grid_tests_radial.sh <run> > ~/mira_dr4_runs/grid_tests_radial_<run>.log 2>&1 &
set -u
cd "$(dirname "$0")"
RUN=$1
PY=$HOME/miniforge3/envs/astro312/bin/python
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-$HOME/mira_dr4_runs}
RS_WORKERS=${RS_WORKERS:-6}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
LOGS=$MIRA_DR4_RUNS/$RUN/logs; mkdir -p "$LOGS"
for SCH in geo20 geo32 ilr01; do
  echo "=== scheme $SCH, |z| < 2: re-score $(date '+%T') ==="
  pids=""; w=0
  while [ "$w" -lt "$RS_WORKERS" ]; do
    OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 \
      --cadence-step 1 --ring-scheme "$SCH" --z-max 2 --layer-width 2 --stride "$RS_WORKERS" --offset "$w" \
      > "$LOGS/rescore_grid_${SCH}_w$w.log" 2>&1 &
    pids="$pids $!"; w=$(( w + 1 ))
  done
  for p in $pids; do wait "$p"; done
  "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 1 \
    --ring-scheme "$SCH" --z-max 2 --layer-width 2 --merge | tail -1
  OMP_NUM_THREADS=3 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" --column "logLs_${B}_sfratio_${SCH}_z2_dz2" \
    --loo-top 32 > "$LOGS/emulator_grid_${SCH}.log" 2>&1 && echo "emulator OK $SCH" || echo "emulator FAILED $SCH"
done
echo "ALL DONE $(date '+%T')"
