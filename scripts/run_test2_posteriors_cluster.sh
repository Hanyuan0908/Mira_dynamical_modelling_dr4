#!/bin/bash
# Run on the big cluster:   sh run_test2_posteriors_cluster.sh
#
# Does the t = 3 and t = 5.5 selection-function posteriors and leaves the rest
# to the workstation, which is working the same shared directories. Nothing is duplicated: every
# design point is claimed with an atomic lock file, and a point already finished anywhere is skipped,
# so it is also safe to run this twice, or to run the full queue here instead (drop ONLY= below).
#
# Needs the same shared filesystem (/data/hz420/...), nothing else. Outputs land beside the
# workstation's, in $MIRA_DR4_RUNS (default <repo>/runs) and <repo>/results.
#
# Knobs, if the defaults are wrong for the cluster:
#   PYTHON=/path/to/python   THREADS=8   WORKERS=12   sh run_test2_posteriors_cluster.sh
#
# The distance-error level is a parameter (defaults = the 10% campaign):
#   SIGMA=0.271 TAG=12p5pct sh run_test2_posteriors_cluster.sh
set -u
cd "$(dirname "$0")"
export PYTHON=${PYTHON:-/data/hz420-2/astro312/bin/python}
export THREADS=${THREADS:-8}
export WORKERS=${WORKERS:-$(( $(nproc) / ${THREADS} ))}
export SIGMA=${SIGMA:-0.217}
export TAG=${TAG:-10pct}
export ONLY=${ONLY:-"test2_t3_sf_${TAG} test2_t5.5_sf_${TAG}"}

echo "cluster run: $(nproc) cores -> $WORKERS workers x $THREADS threads"
echo "runs: $ONLY"
"$PYTHON" -c "import agama, sklearn, healpy, h5py; print('environment ok')" || {
  echo "ERROR: the astro312 environment is not usable here; set PYTHON=..." >&2; exit 1; }
exec bash ./run_test2_posteriors.sh
