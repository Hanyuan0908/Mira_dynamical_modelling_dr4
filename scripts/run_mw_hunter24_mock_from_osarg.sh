#!/bin/bash
# Apple-to-apple test of the OSARG result (user, 2026-11-23): "if the OSARGs were sampled from a simulation of the
# Milky Way, would they constrain the pattern speed?" The mock (make_mw_hunter24_mock_from_osarg.py) is the 26,630
# OSARG stars phase-mixed for 6 Gyr in the Hunter+24 family at |Omega_b| 32.5, phi_b 25, S = bar_frac = 1, observed
# with OSARG's sigma_DM 0.2084 and selected to OSARG's observed spatial distribution and N. It is run here EXACTLY as
# the OSARG data were: galaxy mw, the OSARG 384-point design (results/osarg_zhang2024_mw_hunter24/design.csv), 50 Myr,
# sfratio with K = 10, scheme B emulator. Only the data's origin differs.
#   caffeinate -i -s nohup bash run_mw_hunter24_mock_from_osarg.sh > ~/mira_dr4_runs/mw_hunter24_mock_from_osarg.log 2>&1 &
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-$HOME/miniforge3/envs/astro312/bin/python}
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-$HOME/mira_dr4_runs}
RUN=mw_hunter24_mock_from_osarg_ics
MOCK=mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5
DESIGN=results/osarg_zhang2024_mw_hunter24/design.csv
THREADS=${THREADS:-7}; WORKERS=${WORKERS:-2}; RS_WORKERS=${RS_WORKERS:-6}
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
n=$(ls "$MIRA_DR4_RUNS/$RUN"/fourier_[0-9][0-9][0-9].npz 2>/dev/null | wc -l | tr -d ' ')
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
n=$(ls "$MIRA_DR4_RUNS/$RUN"/sfratio_[0-9][0-9][0-9].npz 2>/dev/null | wc -l | tr -d ' ')
echo "re-score: $n/$N"
if [ "$n" -ne "$N" ]; then echo "FAILED re-score"; grep -h -E "Error|Trace" "$LOGS"/rescore_w*.log | head -5; exit 1; fi
"$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 1 --merge | tail -1
OMP_NUM_THREADS=8 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" --column "logLs_${B}_sfratio" --loo-top 32 \
  > "$LOGS/emulator.log" 2>&1 || { echo "FAILED emulator"; tail -3 "$LOGS/emulator.log"; exit 1; }
echo "ALL OK $(date '+%F %T')"
