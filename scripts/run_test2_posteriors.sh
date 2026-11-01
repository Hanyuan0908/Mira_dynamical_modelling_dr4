#!/bin/bash
# Test 2 posteriors: Sobol search + GP emulator for the six 10%-distance-error mocks.
#
# Runs the searches in the agreed order, then the GP posterior and corner plot for each.
# Safe to run on two machines at once against the same output directory: every design point is
# claimed with an atomic lock file, so workers never duplicate work, and anything already finished
# is skipped. Just run it; it picks up wherever the other machine got to.
#
#   sh run_test2_posteriors.sh              # workers = nproc/8
#   WORKERS=12 sh run_test2_posteriors.sh   # force a worker count
#   PYTHON=/path/to/python sh run_test2_posteriors.sh
#
# Outputs: $MIRA_DR4_RUNS/<run>/{traj_NNN.npy, fourier_NNN.npz}
#          results/<run>/{evals.csv, run_meta.json, emulator_results_*.json}
#          figures/<run>_emulator_{corner,slices,loo}_*.png
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
CORES=$(nproc)
THREADS=${THREADS:-8}                       # OpenMP threads per worker
WORKERS=${WORKERS:-$(( CORES / THREADS ))}
[ "$WORKERS" -lt 1 ] && WORKERS=1
SCHEME=fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024
echo "=== $CORES cores -> $WORKERS workers x $THREADS threads ==="

# run name | mock file | galaxy | emulator columns (library variants to build posteriors for)
RUNS="
test2_t7_nosf_10pct|nbody_mock_distance_error_sigmaDM0.217_R8_N50k.hdf5|nbody|observed asis
test2_t7_sf_10pct|nbody_mock_selection_sigmaDM0.217_N50k.hdf5|nbody|sf asis
test2_t2_sf_10pct|nbody_t2_mock_selection_sigmaDM0.217_N50k.hdf5|nbody_t2|sf
test2_t4_sf_10pct|nbody_t4_mock_selection_sigmaDM0.217_N50k.hdf5|nbody_t4|sf
test2_t5.5_sf_10pct|nbody_t5.5_mock_selection_sigmaDM0.217_N50k.hdf5|nbody_t5.5|sf
test2_t3_sf_10pct|nbody_t3_mock_selection_sigmaDM0.217_N50k.hdf5|nbody_t3|sf
"

# ONLY="name1 name2" restricts the queue to those runs (substring match on the run name)
ONLY=${ONLY:-}

echo "$RUNS" | while IFS='|' read -r RUN MOCK GAL COLS; do
  [ -z "${RUN:-}" ] && continue
  if [ -n "$ONLY" ]; then
    skip=1
    for want in $ONLY; do case "$RUN" in *"$want"*) skip=0 ;; esac; done
    [ "$skip" -eq 1 ] && continue
  fi
  echo ""
  echo "########## $RUN  ($MOCK, galaxy $GAL)  $(date '+%F %T') ##########"

  # The trial potentials are built from the snapshot and are not in git, so make sure this galaxy has
  # them before launching any worker - otherwise every worker dies instantly on a missing .ini.
  # Idempotent: skips in seconds when the files are already there.
  if ! OMP_NUM_THREADS=$CORES $PY -c "
import sys, steady_state_likelihood as ssl
ssl.bar_aligned_potential(1.0, 1.0, '$GAL')" >/dev/null 2>&1; then
    echo "    building the potential for $GAL (one-off, a few minutes) ..."
    OMP_NUM_THREADS=$CORES $PY -W ignore build_nbody_potential.py --galaxy "$GAL" > "/tmp/${RUN}_pot.log" 2>&1       && OMP_NUM_THREADS=$CORES $PY -W ignore potential_family.py --galaxy "$GAL" >> "/tmp/${RUN}_pot.log" 2>&1       || { echo "    ERROR: could not build the potential for $GAL (see /tmp/${RUN}_pot.log); skipping $RUN"; continue; }
    OMP_NUM_THREADS=$CORES $PY -c "
import steady_state_likelihood as ssl
ssl.bar_aligned_potential(1.0, 1.0, '$GAL')" >/dev/null 2>&1       || { echo "    ERROR: potential for $GAL still unusable; skipping $RUN"; continue; }
    echo "    potential for $GAL ready"
  fi

  pids=""
  w=0
  while [ "$w" -lt "$WORKERS" ]; do
    OMP_NUM_THREADS=$THREADS $PY -W ignore test2b_qmc_distance_error_run.py \
      --mock "$MOCK" --run "$RUN" --galaxy "$GAL" --stride "$WORKERS" --offset "$w" \
      > "/tmp/${RUN}_w${w}.log" 2>&1 &
    pids="$pids $!"
    w=$(( w + 1 ))
  done
  for pid in $pids; do wait "$pid"; done
  OMP_NUM_THREADS=$THREADS $PY -W ignore test2b_qmc_distance_error_run.py \
    --mock "$MOCK" --run "$RUN" --galaxy "$GAL" --merge
  echo "--- $RUN search done $(date '+%F %T'); posteriors ---"
  for V in $COLS; do
    OMP_NUM_THREADS=$THREADS $PY -W ignore test1_gp_emulator.py \
      --run "$RUN" --column "logLs_${SCHEME}_${V}" --loo-top 32 \
      > "/tmp/${RUN}_emu_${V}.log" 2>&1 \
      && echo "    posterior ok: $RUN / $V" || echo "    posterior FAILED: $RUN / $V (see /tmp/${RUN}_emu_${V}.log)"
    OMP_NUM_THREADS=$THREADS $PY -W ignore plot_test1_logl_corner.py \
      --run "$RUN" --column "logLs_${SCHEME}_${V}" \
      --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2 \
      > "/tmp/${RUN}_corner_${V}.log" 2>&1 || true
  done
done
echo ""
echo "=== all done $(date '+%F %T') ==="
