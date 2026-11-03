#!/bin/bash
# ============================================================================================
#  Test 2 at a 12.5% distance error - the SEVEN FILL-IN SNAPSHOTS.   Run:  sh <this script>
# ============================================================================================
#
# Turns the five points on the recovery plot into twelve, so the systematic error in Omega_b can
# be separated from snapshot-to-snapshot scatter. Nothing about the method changes - same 257-point
# focused Sobol design, same DR3-Mira-like selection function, same chi2 schemes, same sigma_DM
# (0.2714 mag = 12.5%). Only the galaxy changes.
#
#     t       Omega_b   A2      placed in
#     2.375     40.87  0.416    t = 2 .. 3
#     2.625     40.21  0.442    t = 2 .. 3
#     3.375     36.83  0.459    t = 3 .. 4
#     3.625     34.79  0.450    t = 3 .. 4
#     4.5       32.07  0.470    t = 4 .. 5.5
#     5.0       30.60  0.461    t = 4 .. 5.5
#     6.25      26.58  0.433    t = 5.5 .. 7
#
# With t = 2, 3, 4, 5.5, 7 already done: twelve points from 25.2 to 43.9, largest gap 3.0.
#
# --------------------------------------------------------------------------------------------
# WHAT IT NEEDS
#   - the shared filesystem (/data/hz420/...) and the astro312 environment. Nothing else.
#   - the seven mocks and their trial-potential families, which are already built in the shared
#     tree. The script checks and refuses to start a run whose mock is missing, rather than
#     letting every worker die on it.
#
# WHAT IT WRITES  (shared tree, never a removable disk)
#   runs/<run>/{traj_NNN.npy, fourier_NNN.npz}        ~6.4 GB per run (see CADENCE below)
#   results/<run>/{evals.csv, run_meta.json, emulator_results_*.json}
#   figures/<run>_emulator_*.png
# The workstation collects each finished run to the portable disk and frees the space here, so
# the shared disk should not fill. If it does, the runner stops caching trajectories below 25 GB
# free and keeps going - the results stay correct, only the re-scoring cache is lost.
#
# TIME
#   ~400 s of single-core integration per point on this cluster, 257 points per run, 7 runs.
#   At WORKERS = 11 that was 2.7 h per run on the 10% campaign, so expect ~19 h for all seven.
#   The cadence does not change the integration cost - agama integrates the same orbits and only
#   the output sampling is thinned - so it buys disk, not time.
#   Resumable and restart-safe: every point is claimed with an atomic O_EXCL lock and anything
#   already finished anywhere is skipped. Killing this and re-running it loses nothing.
#
# KNOBS
#   PYTHON=/path/to/python  THREADS=8  WORKERS=12  sh run_test2_12p5pct_fill_cluster.sh
#   ONLY="test2_t4.5_sf_12p5pct"      sh ...   # just one run
#   TIME_BUDGET_H=10                  sh ...   # stop cleanly after N hours
# --------------------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")"

export PYTHON=${PYTHON:-/data/hz420-2/astro312/bin/python}
export THREADS=${THREADS:-8}
export WORKERS=${WORKERS:-$(( $(nproc) / ${THREADS} ))}
[ "$WORKERS" -lt 1 ] && export WORKERS=1
export SIGMA=0.271
export TAG=12p5pct
export CADENCE=${CADENCE:-2}     # keep every 2nd snapshot: 50 Myr instead of 25, half the bytes.
                                 # The library is scored on the thinned trajectory and only the
                                 # thinned trajectory is stored, so the saved file and any later
                                 # rescoring are on the same grid. Measured cost on the posterior:
                                 # 0.06-0.17 km/s/kpc, a tenth of the error bar; chi2/N at the truth
                                 # rises 4-12% (scripts/test_time_cadence.py). CADENCE=1 restores
                                 # the 25 Myr behaviour of the earlier campaigns.
export ONLY=${ONLY:-"test2_t2.375_sf_${TAG} test2_t2.625_sf_${TAG} test2_t3.375_sf_${TAG} \
test2_t3.625_sf_${TAG} test2_t4.5_sf_${TAG} test2_t5_sf_${TAG} test2_t6.25_sf_${TAG}"}

echo "=============================================================="
echo " Test 2 fill-in, 12.5% distance error"
echo " host    $(hostname)   $(nproc) cores -> $WORKERS workers x $THREADS threads"
echo " runs    $(echo $ONLY | wc -w)"
echo " cadence every ${CADENCE} snapshots -> $(( (41 - 1) / CADENCE + 1 )) of 41, $(( 25 * CADENCE )) Myr"
echo " started $(date '+%F %T')"
echo "=============================================================="

"$PYTHON" -c "import agama, sklearn, healpy, h5py, numpy; print('environment ok, numpy', numpy.__version__)" || {
  echo "ERROR: the astro312 environment is not usable here; set PYTHON=/path/to/python" >&2; exit 1; }

missing=0
for R in $ONLY; do
  T=${R#test2_t}; T=${T%%_sf*}
  M="../mock/nbody_t${T}_mock_selection_sigmaDM${SIGMA}_N50k.hdf5"
  [ -f "$M" ] || { echo "ERROR: missing mock for $R: $M" >&2; missing=1; }
done
[ "$missing" -eq 0 ] || { echo "refusing to start with missing mocks" >&2; exit 1; }
echo "all seven mocks present"
echo ""

exec bash ./run_test2_posteriors.sh
