#!/bin/bash
# One Fourier grid on one run: re-score the saved trajectories with a named radial scheme and given |z| layers, merge
# once after the workers, refit the pipeline emulator (diagnostics). No integration.
#   bash run_grid_test_scheme.sh <run> <ring_scheme> <z_max> <layer_width>
set -u
cd "$(dirname "$0")"
RUN=$1; SCH=$2; ZMAX=$3; LW=$4
PY=$HOME/miniforge3/envs/astro312/bin/python
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-$HOME/mira_dr4_runs}
RS_WORKERS=${RS_WORKERS:-6}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
TAG=$("$PY" -c "import test1_nbody_qmc_focused_fourier_stars as fs; print(fs.set_grid(0.3, $ZMAX, $LW, '$SCH'))")
LOGS=$MIRA_DR4_RUNS/$RUN/logs; mkdir -p "$LOGS"
echo "=== $RUN: scheme $SCH, z_max $ZMAX, layer width $LW (tag $TAG): re-score $(date '+%T') ==="
pids=""; w=0
while [ "$w" -lt "$RS_WORKERS" ]; do
  OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 \
    --cadence-step 1 --ring-scheme "$SCH" --z-max "$ZMAX" --layer-width "$LW" --stride "$RS_WORKERS" --offset "$w" \
    > "$LOGS/rescore_grid${TAG}_w$w.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
"$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 1 \
  --ring-scheme "$SCH" --z-max "$ZMAX" --layer-width "$LW" --merge | tail -1
OMP_NUM_THREADS=3 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" --column "logLs_${B}_sfratio${TAG}" \
  --loo-top 32 > "$LOGS/emulator_grid${TAG}.log" 2>&1 && echo "emulator OK $TAG" || echo "emulator FAILED $TAG"
echo "DONE $RUN $TAG $(date '+%T')"
