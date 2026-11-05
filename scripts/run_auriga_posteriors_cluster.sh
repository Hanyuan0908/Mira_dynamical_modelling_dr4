#!/bin/bash
# Auriga dynamical-model campaign, BATCH 1 -- the cluster half.
#
#   sh run_auriga_posteriors_cluster.sh --check     # verify prerequisites and exit (seconds)
#   sh run_auriga_posteriors_cluster.sh             # run batch 1
#
# Batch 1 is the seven galaxies agreed for the first pass:
#     Au9  Au10  Au17  Au18  Au24  Au26  Au27
# Batch 2, when you want it:
#     HALOS="5 12 14 20 22 28" sh run_auriga_posteriors_cluster.sh
#
# Identical method and identical 257-point focused Sobol design to the Test 2 N-body campaign --
# results/test1_nbody_qmc_focused/design.csv, the same file the N-body runs read. The only
# per-galaxy difference is the one the pipeline already applies to every N-body snapshot: the
# Omega axis of the design is rescaled by the galaxy's own pattern speed / 26.02
# (test2b_qmc_distance_error_run.scaled_design). Nothing else about the method changes.
#
# PREREQUISITES on whatever machine this runs on:
#   - the shared tree /data/hz420/mira_dynamics_dr4/Mira_dynamical_modelling_dr4 (mocks, potentials,
#     ground truth and the design all live there; this script writes back into it)
#   - the astro312 environment, /data/hz420-2/astro312/bin/python, or set PYTHON=...
# Run --check first: it verifies both, plus the mocks, potentials, ground truth and the frame
# convention, in a few seconds, and launches nothing.
#
# This machine does NOT archive -- the USB disk is on the workstation. Run
#   HALOS="9 10 17 18 24 26 27" sh archive_auriga_runs.sh --delete
# on the workstation at the same time; it watches for finished runs, copies them to the portable
# disk, verifies, and only then frees the trajectories here. Without it the campaign fills
# /data/hz420 after ~4 galaxies and blocks (it waits, it does not crash or lose work).
#
# Knobs:  PYTHON=...  THREADS=8  WORKERS=11  FREE_MIN_GB=20  HALOS="9 10"  ONLY="au9"
#
# On worker geometry: the per-point cost is ~235 s of single-core integration and the OpenMP
# threading inside a point scales poorly, so wall-clock divides by WORKERS, not by
# WORKERS x THREADS. If the cluster has cores to spare, THREADS=2 with a correspondingly larger
# WORKERS finishes a galaxy several times faster than THREADS=8. It changes no numbers.
set -u
cd "$(dirname "$0")"
export PYTHON=${PYTHON:-/data/hz420-2/astro312/bin/python}
export THREADS=${THREADS:-8}
export WORKERS=${WORKERS:-$(( $(nproc) / ${THREADS} ))}
export SIGMA=${SIGMA:-0.271}
export TAG=${TAG:-12p5pct}
export ARCHIVE=0                                    # no portable disk here, ever
export HALOS=${HALOS:-"9 10 17 18 24 26 27"}        # batch 1
export ONLY=${ONLY:-}

echo "cluster run: $(nproc) cores -> $WORKERS workers x $THREADS threads"
echo "batch: $HALOS"
"$PYTHON" -c "import agama, sklearn, healpy, h5py; print('environment ok')" || {
  echo "ERROR: the astro312 environment is not usable here; set PYTHON=..." >&2; exit 1; }
exec bash ./run_auriga_posteriors.sh "$@"
