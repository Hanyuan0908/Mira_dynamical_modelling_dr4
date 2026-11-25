#!/bin/bash
# Is OSARG's flat likelihood a potential-family mismatch? The t = 7 mock (12.5%, Mira-like selection; truth
# |Omega_b| 25.20, phi_b 27) integrated in the Milky Way family (Hunter+24, galaxy "mw") on EXACTLY the OSARG
# design (results/osarg_zhang2024_mw_hunter24/design.csv: 384 points, |Omega_b| 15-45, phi_b 20-30, S and
# bar_frac 0.5-1.5), then the adopted re-score (sfratio, K = 10, 50 Myr) and the scheme B emulator. Only the data
# differ from the OSARG run. The family mismatch here is large: Hunter+24's v_c is 25-35% below t = 7's inside 2 kpc.
#   caffeinate -i -s nohup bash run_t7_in_mw_hunter24_family.sh > ~/mira_dr4_runs/t7_mwH24.log 2>&1 &
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-$HOME/miniforge3/envs/astro312/bin/python}
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs}
RUN=test2_t7_sf_12p5pct_mwH24_osargdesign
MOCK=nbody_mock_selection_sigmaDM0.271_N50k.hdf5
DESIGN=results/osarg_zhang2024_mw_hunter24/design.csv
THREADS=${THREADS:-7}; WORKERS=${WORKERS:-2}; RS_WORKERS=${RS_WORKERS:-5}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
N=$(( $(wc -l < "../$DESIGN") - 1 ))
LOGS=$MIRA_DR4_RUNS/$RUN/logs; mkdir -p "$LOGS"

echo "=== $RUN: search, $N points  $(date '+%F %T') ==="
pids=""; w=0
while [ "$w" -lt "$WORKERS" ]; do
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py --mock "$MOCK" --run "$RUN" \
    --galaxy mw --design "$DESIGN" --cadence-step 2 --stride "$WORKERS" --offset "$w" > "$LOGS/search_w$w.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
n=$(ls "$MIRA_DR4_RUNS/$RUN" 2>/dev/null | grep -E "^fourier_[0-9]+\.npz$" | wc -l | tr -d ' ')
echo "search: $n/$N"
if [ "$n" -ne "$N" ]; then echo "FAILED search"; grep -h -E "Error|Trace" "$LOGS"/search_w*.log | head -5; exit 1; fi
"$PY" -W ignore test2b_qmc_distance_error_run.py --mock "$MOCK" --run "$RUN" --galaxy mw --design "$DESIGN" \
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
