#!/bin/bash
# Re-score every run with saved libraries without the <v_phi> m = 0 term (rescore_without_vphi_m0.py), then refit
# its emulator on that column. No integration. User, 2026-11-24.
#   caffeinate -i -s nohup bash run_without_vphi_m0_all.sh > ~/mira_dr4_runs/without_vphi_m0.log 2>&1 &
set -u
cd "$(dirname "$0")"
PY=$HOME/miniforge3/envs/astro312/bin/python
USB=/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs
LOCAL=$HOME/mira_dr4_runs
COL=logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m24
LOGS=$LOCAL/without_vphi_m0_logs; mkdir -p "$LOGS"
for g in t2 t2.375 t2.625 t3 t3.375 t3.625 t4 t4.5 t5 t5.5 t6.25 t7 au5 au12 au14 au20 au22 au28; do
  OMP_NUM_THREADS=2 "$PY" -W ignore rescore_without_vphi_m0.py --run test2_${g}_sf_12p5pct --score-dir "$USB/test2_${g}_sf_12p5pct" | head -1
done
OMP_NUM_THREADS=2 "$PY" -W ignore rescore_without_vphi_m0.py --run test2_t7_sf_12p5pct --score-dir "$USB/test2_t7_sf_12p5pct" \
  --variant sfratio_osarg_cover --results-run test2_t7_sf_12p5pct_osarg_cover --mock nbody_mock_selection_sigmaDM0.271_N50k.hdf5 \
  --keep results/rv_selection_masks/nbody_mock_selection_sigmaDM0.271_N50k__osarg_cover.npy | head -1
OMP_NUM_THREADS=2 "$PY" -W ignore rescore_without_vphi_m0.py --run test2_t7_sf_12p5pct_mwH24_osargdesign \
  --score-dir "$USB/test2_t7_sf_12p5pct_mwH24_osargdesign" | head -1
OMP_NUM_THREADS=2 "$PY" -W ignore rescore_without_vphi_m0.py --run mw_hunter24_mock_from_osarg_ics --score-dir "$LOCAL/mw_hunter24_mock_from_osarg_ics" | head -1
OMP_NUM_THREADS=2 "$PY" -W ignore rescore_without_vphi_m0.py --run osarg_zhang2024_mw_hunter24 --score-dir "$LOCAL/osarg_zhang2024_mw_hunter24" | head -1
echo "=== re-scored $(date '+%T'); emulators ==="
jobs_list() {
  for g in t2 t2.375 t2.625 t3 t3.375 t3.625 t4 t4.5 t5 t5.5 t6.25 t7 au5 au12 au14 au20 au22 au28; do
    echo "test2_${g}_sf_12p5pct ${COL}_sfratio"
  done
  echo "test2_t7_sf_12p5pct_osarg_cover ${COL}_sfratio_osarg_cover"
  echo "test2_t7_sf_12p5pct_mwH24_osargdesign ${COL}_sfratio"
  echo "mw_hunter24_mock_from_osarg_ics ${COL}_sfratio"
  echo "osarg_zhang2024_mw_hunter24 ${COL}_sfratio"
}
jobs_list | xargs -P 4 -L 1 sh -c 'OMP_NUM_THREADS=2 '"$PY"' -W ignore test1_gp_emulator.py --run "$0" --column "$1" --loo-top 32 > '"$LOGS"'/emulator_"$0".log 2>&1 && echo "emulator OK $0" || echo "emulator FAILED $0"'
echo "ALL DONE $(date '+%T')"
