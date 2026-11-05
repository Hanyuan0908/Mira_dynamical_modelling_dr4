#!/bin/bash
# Auriga BATCH 2 -- the cluster half.   Au5, Au12, Au14, Au20, Au22, Au28.
#
#   sh run_auriga_batch2_cluster.sh --check     # verify prerequisites and exit (seconds)
#   sh run_auriga_batch2_cluster.sh             # run batch 2
#
# This does the SEARCH ONLY: 257-point focused Sobol design per galaxy, then the merge to
# results/<run>/evals.csv. It deliberately does NOT build any posterior.
#
# Why no emulator here: the search produces `observed`, `asis` and `sf` on chi2 scheme A, and none
# of those is the answer. The result column is scheme B + sfratio + K = 10, which only exists after
# rescore_empirical_selection.py (see the 2026-11-05 PIPELINE entry in PROGRESS.md). Building GP
# posteriors on the search columns is what produced the misleading +7.7% Omega_b bias on batch 1.
# The workstation does the re-score, the posteriors and the archiving -- postprocess_auriga_batch2.sh.
#
# Differences from batch 1, both deliberate:
#   * integrated at CADENCE=2 (50 Myr, 21 of 41 snapshots) rather than integrated at 25 Myr and
#     thinned at scoring. The scored time grid is identical either way; this just halves the bytes,
#     ~6.3 GB a run instead of ~12.3.
#   * no GP step, as above.
#
# PREREQUISITES on whatever machine this runs on:
#   - the shared tree /data/hz420/mira_dynamics_dr4/Mira_dynamical_modelling_dr4
#   - the astro312 environment, /data/hz420-2/astro312/bin/python, or set PYTHON=...
# --check verifies both, plus the mocks, the ground-truth potentials and the frame convention
# (clockwise, bar at 27 deg, Omega_b negative), and launches nothing.
#
# The workstation drains finished runs to the portable disk as they complete, so /data does not
# fill. If that is not running, this blocks below FREE_MIN_GB rather than crashing.
#
# Knobs:  PYTHON=...  THREADS=8  WORKERS=11  FREE_MIN_GB=20  HALOS="5 12"  ONLY="au5"
#
# Wall-clock note: the per-point cost is ~235 s of single-core integration and the OpenMP threading
# inside a point scales poorly, so time divides by WORKERS, not by WORKERS x THREADS. If the cluster
# has cores spare, THREADS=2 with a larger WORKERS is several times faster and changes no numbers.
set -u
cd "$(dirname "$0")"
export PYTHON=${PYTHON:-/data/hz420-2/astro312/bin/python}
export THREADS=${THREADS:-8}
export WORKERS=${WORKERS:-$(( $(nproc) / ${THREADS} ))}
export SIGMA=${SIGMA:-0.271}
export TAG=${TAG:-12p5pct}
export CADENCE=${CADENCE:-2}                  # integrate and save at 50 Myr
export ARCHIVE=0                              # no portable disk here
export EMULATE=0                              # search only; posteriors are the workstation's job
export HALOS=${HALOS:-"5 12 14 20 22 28"}     # batch 2
export ONLY=${ONLY:-}

echo "cluster run: $(nproc) cores -> $WORKERS workers x $THREADS threads"
echo "batch 2: $HALOS   (cadence step $CADENCE = $(( 25 * CADENCE )) Myr, search only, no emulator)"
"$PYTHON" -c "import agama, sklearn, healpy, h5py; print('environment ok')" || {
  echo "ERROR: the astro312 environment is not usable here; set PYTHON=..." >&2; exit 1; }
exec bash ./run_auriga_posteriors.sh "$@"
