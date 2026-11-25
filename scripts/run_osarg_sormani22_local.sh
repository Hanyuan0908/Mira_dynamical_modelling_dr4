#!/bin/bash
# The OSARG data (Zhang et al. 2024) in the Sormani+22 Milky Way family instead of Hunter+24 (user, 2026-11-24:
# "test Sormani et al. 2022 anyway"). Same as run_osarg_mw_local.sh in every other respect: the OSARG 384-point
# design (results/osarg_zhang2024_mw_hunter24/design.csv), 50 Myr, sfratio with K = 10, scheme B emulator.
# Family: build_sormani22_potential.py (Portail17.ini with its halo made spherical).
#   caffeinate -i -s nohup bash run_osarg_sormani22_local.sh > ~/mira_dr4_runs/osarg_sormani22.log 2>&1 &
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-$HOME/miniforge3/envs/astro312/bin/python}
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-$HOME/mira_dr4_runs}
RUN=${RUN:-osarg_zhang2024_sormani22}   # RUN=osarg_zhang2024_sormani22_widebar for the wider bar box (2026-11-25)
MOCK=osarg_zhang2024_observed_R8.hdf5
DESIGN=${DESIGN:-results/osarg_zhang2024_mw_hunter24/design.csv}   # the shared OSARG design
THREADS=${THREADS:-7}; WORKERS=${WORKERS:-2}; RS_WORKERS=${RS_WORKERS:-6}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
N=$(( $(wc -l < "../$DESIGN") - 1 ))
LOGS=$MIRA_DR4_RUNS/$RUN/logs; mkdir -p "$LOGS"

echo "=== $RUN: search, $N points  $(date '+%F %T') ==="
pids=""; w=0
while [ "$w" -lt "$WORKERS" ]; do
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py --mock "$MOCK" --run "$RUN" \
    --galaxy sormani22 --design "$DESIGN" --cadence-step 2 --stride "$WORKERS" --offset "$w" > "$LOGS/search_w$w.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
n=$(ls "$MIRA_DR4_RUNS/$RUN" 2>/dev/null | grep -E "^fourier_[0-9]+\.npz$" | wc -l | tr -d ' ')
echo "search: $n/$N"
if [ "$n" -ne "$N" ]; then echo "FAILED search"; grep -h -E "Error|Trace" "$LOGS"/search_w*.log | head -5; exit 1; fi
"$PY" -W ignore test2b_qmc_distance_error_run.py --mock "$MOCK" --run "$RUN" --galaxy sormani22 --design "$DESIGN" \
  --cadence-step 2 --merge | tail -1

pids=""; w=0
while [ "$w" -lt "$RS_WORKERS" ]; do
  OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 \
    --cadence-step 1 --stride "$RS_WORKERS" --offset "$w" > "$LOGS/rescore_w$w.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
n=$(ls "$MIRA_DR4_RUNS/$RUN" 2>/dev/null | grep -E "^sfratio_[0-9]+\.npz$" | wc -l | tr -d ' ')
echo "re-score: $n/$N"
if [ "$n" -ne "$N" ]; then echo "FAILED re-score"; grep -h -E "Error|Trace" "$LOGS"/rescore_w*.log | head -5; exit 1; fi
"$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 1 --merge | tail -1
OMP_NUM_THREADS=8 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" --column "logLs_${B}_sfratio" --loo-top 32 \
  > "$LOGS/emulator.log" 2>&1 || { echo "FAILED emulator"; tail -3 "$LOGS/emulator.log"; exit 1; }
echo "ALL OK $(date '+%F %T')"
