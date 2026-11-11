#!/bin/bash
# The real data, first pass (user, 2026-11-11: a preliminary result for a public talk): the Zhang et al. (2024)
# OSARG sample with the adopted pipeline, in the Milky Way trial-potential family (Hunter+24).
#
#   caffeinate -i -s nohup sh run_osarg_mw_local.sh > ~/mira_dr4_runs/osarg_zhang2024_mw_hunter24.log 2>&1 &
#
# Inputs, made once:
#   python build_mw_potential.py          # mock/mw_potential/ (Hunter+24 family, bar along x)
#   python make_osarg_data.py             # mock/osarg_zhang2024_observed_R8.hdf5 (26,630 stars, sigma_DM 0.2084)
#   python make_sobol_design.py --run osarg_zhang2024_mw_hunter24 --n 256 --seed 2041 \
#       --omega 25 45 --phi 20 30 --S 0.5 1.5 --bar-frac 0.5 1.5
# Configuration as PROGRESS 2026-11-05 ("THE PIPELINE"): 50 Myr trajectories (--cadence-step 2), sfratio
# (rho_obs/rho_lib per HEALPix nside 32 x 1 kpc cell) with K = 10 distance draws, scheme B, tempered logLs.
# Measured on the M5 Pro: 4 s integrate per point (26.6k stars), 6 s and 1.8 GB per re-scored point.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-$HOME/miniforge3/envs/astro312/bin/python}
export MIRA_DR4_RUNS=${MIRA_DR4_RUNS:-$HOME/mira_dr4_runs}
RUN=osarg_zhang2024_mw_hunter24
MOCK=osarg_zhang2024_observed_R8.hdf5
DESIGN=results/$RUN/design.csv
THREADS=${THREADS:-7}; WORKERS=${WORKERS:-2}; RS_WORKERS=${RS_WORKERS:-6}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
LOGS=$MIRA_DR4_RUNS/$RUN/logs; mkdir -p "$LOGS"

echo "=== $RUN: search, $(wc -l < ../$DESIGN) design rows  $(date '+%F %T') ==="
pids=""; w=0
while [ "$w" -lt "$WORKERS" ]; do
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py --mock "$MOCK" --run "$RUN" \
    --galaxy mw --design "$DESIGN" --cadence-step 2 --stride "$WORKERS" --offset "$w" > "$LOGS/search_w$w.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
"$PY" -W ignore test2b_qmc_distance_error_run.py --mock "$MOCK" --run "$RUN" --galaxy mw --design "$DESIGN" \
  --cadence-step 2 --merge
echo "=== search done $(date '+%F %T'): $(ls "$MIRA_DR4_RUNS/$RUN"/fourier_*.npz | wc -l) points ==="

# the trajectories are already 50 Myr, so the re-score uses step 1
pids=""; w=0
while [ "$w" -lt "$RS_WORKERS" ]; do
  OMP_NUM_THREADS=1 "$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 \
    --cadence-step 1 --stride "$RS_WORKERS" --offset "$w" > "$LOGS/rescore_w$w.log" 2>&1 &
  pids="$pids $!"; w=$(( w + 1 ))
done
for p in $pids; do wait "$p"; done
"$PY" -W ignore rescore_empirical_selection.py --run "$RUN" --estimator ratio --n-draw 10 --cadence-step 1 --merge
echo "=== re-scored $(ls "$MIRA_DR4_RUNS/$RUN"/sfratio_*.npz | wc -l) points  $(date '+%F %T') ==="

OMP_NUM_THREADS=8 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" --column "logLs_${B}_sfratio" --loo-top 32
echo "=== ALL DONE $(date '+%F %T') ==="
