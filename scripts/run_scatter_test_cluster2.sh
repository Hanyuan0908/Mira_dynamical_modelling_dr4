#!/bin/bash
# Per-star distance-error test, SECOND BATCH: 12 galaxies, for a second cluster.
#
#   sh run_scatter_test_cluster2.sh --check     # verify prerequisites and exit
#   nohup sh run_scatter_test_cluster2.sh > ../logs/scatter2.log 2>&1 &
#
# Same test as run_scatter_test_cluster.sh, different galaxies. The data has a per-star distance
# error, sigma_i ~ N(0.2714, 0.0543) mag = 12.5% +- 2.5%; the mock's sigma_DM attribute still
# records the nominal 0.2714, so everything downstream assumes a uniform 12.5%. That mismatch is
# the whole experiment: nothing else differs from the fiducial runs of the same galaxies.
#
# The 12 are chosen to span pattern speed (16.8-56.7), corotation radius (4.8-11.3 kpc) and bar
# strength (A2_max 0.30-0.57), and include both bar-angle outliers: Au5 (weakest bar in the set)
# and Au20 (slowest bar, the one almost immune to distance smearing).
#
# SEARCH ONLY - no emulator. The result column is scheme B + sfratio + K = 10, which only exists
# after rescore_empirical_selection.py; the workstation does that, the posteriors and the
# archiving. Copy results/<run>/ and runs/<run>/ back when the searches finish.
#
# Runs off the shared tree /data/hz420/mira_dynamics_dr4/Mira_dynamical_modelling_dr4, like the
# other cluster - the mocks, trial potentials, truths and design are all already there, and the
# outputs land beside them. The only thing neither cluster can reach is the portable USB disk, so
# archiving stays on the workstation. Nothing needs copying anywhere.
# Needs the astro312 environment; set PYTHON= if it is not at the default path.
#
# Knobs: PYTHON=...  THREADS=8  WORKERS=11  ONLY="au20 au5"
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
ROOT=$(cd .. && pwd)
THREADS=${THREADS:-8}
WORKERS=${WORKERS:-$(( $(nproc) / THREADS ))}
[ "$WORKERS" -lt 1 ] && WORKERS=1
CADENCE=${CADENCE:-2}
SUF=_scat0.054
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
FREE_MIN_GB=${FREE_MIN_GB:-15}
GALS=${ONLY:-"au20 au24 au9 au27 au12 au5 au17 au10 nbody_t6.25 nbody_t4.5 nbody_t3.625 nbody_t2"}

CHECK=0; CLEAR=0
for a in "$@"; do case "$a" in
  --check) CHECK=1 ;;
  --clear-locks) CLEAR=1 ;;
  *) echo "unknown option: $a" >&2; exit 2 ;;
esac; done

# claim() uses an O_EXCL lock with no expiry, so a worker killed mid-point leaves one behind and
# that design point is skipped for ever (it merges as NaN and the emulator dies on it). Only use
# this when NOTHING is running against these runs - it removes every lock with no .npz beside it.
if [ "$CLEAR" = 1 ]; then
  nrm=0
  for G in $GALS; do
    T=$(case "$G" in nbody) echo t7 ;; nbody_t*) echo "t${G#nbody_t}" ;; *) echo "$G" ;; esac)
    D=$RUNS_DIR/test2_${T}_sf_12p5pct_scat
    [ -d "$D" ] || continue
    for L in "$D"/fourier_*.lock; do
      [ -e "$L" ] || continue
      [ -f "${L%.lock}.npz" ] || { rm -f "$L"; nrm=$(( nrm + 1 )); }
    done
  done
  echo "cleared $nrm stale lock(s) with no completed point beside them"
fi
tag_of () { case "$1" in nbody) echo t7 ;; nbody_t*) echo "t${1#nbody_t}" ;; *) echo "$1" ;; esac; }

echo "per-star distance-error test, batch 2: $(printf '%s' "$GALS" | wc -w) galaxies, "\
"cadence $CADENCE ($(( 25 * CADENCE )) Myr), $WORKERS x $THREADS"
"$PY" -c "import agama, sklearn, healpy, h5py; print('environment ok')" || {
  echo "ERROR: astro312 not usable here; set PYTHON=..." >&2; exit 1; }

GALS="$GALS" "$PY" - <<'PRE' || exit 1
import os, sys, h5py, json
from common import MOCK_DIR, RESULTS_DIR
bad = []
for g in os.environ["GALS"].split():
    m = MOCK_DIR / f"{g}_mock_selection_sigmaDM0.271_N50k_scat0.054.hdf5"
    if not m.exists():
        bad.append(f"{g}: mock missing"); continue
    with h5py.File(m, "r") as h:
        if h["DM_err"][:].std() < 1e-6:
            bad.append(f"{g}: DM_err constant - not a per-star mock")
        if abs(float(h.attrs["sigma_DM"]) - 0.2714) > 1e-3:
            bad.append(f"{g}: sigma_DM attr not the nominal 0.2714")
    # Ask for the truth the way the pipeline does: an Auriga halo has its own ground-truth
    # json, an N-body snapshot falls back to results/nbody_pattern_speed_evolution.csv.
    # Checking for the json directly wrongly fails every N-body snapshot.
    try:
        from test2b_qmc_distance_error_run import truth_omega
        truth_omega(g)
    except Exception as e:
        bad.append(f"{g}: no truth pattern speed - {str(e)[:60]}")
    import steady_state_likelihood as ssl
    try:
        ssl.bar_aligned_potential(1.0, 1.0, g)
    except Exception as e:
        bad.append(f"{g}: trial potential unusable - {str(e)[:60]}")
if not (RESULTS_DIR / "test1_nbody_qmc_focused" / "design.csv").exists():
    bad.append("design.csv missing")
if bad:
    sys.stderr.write("ERROR:\n  " + "\n  ".join(bad) + "\n"); sys.exit(1)
print("    mocks per-star, sigma_DM nominal, ground truths and trial potentials all present")
PRE
[ "$CHECK" = 1 ] && { echo "=== preflight passed; nothing launched (--check) ==="; exit 0; }

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }
for G in $GALS; do
  T=$(tag_of "$G")
  RUN=test2_${T}_sf_12p5pct_scat
  MOCK=${G}_mock_selection_sigmaDM0.271_N50k${SUF}.hdf5
  echo ""
  echo "########## $RUN  ($G)  $(date '+%F %T') ##########"
  n=$(find "$RUNS_DIR/$RUN" -name 'fourier_*.npz' 2>/dev/null | wc -l)
  [ "$n" -ge 257 ] && { echo "    already complete"; continue; }
  while [ "$(free_gb "$RUNS_DIR")" -lt "$FREE_MIN_GB" ]; do
    echo "    only $(free_gb "$RUNS_DIR") GB free; waiting 5 min"; sleep 300
  done
  pids=""; w=0
  while [ "$w" -lt "$WORKERS" ]; do
    OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
      --mock "$MOCK" --run "$RUN" --galaxy "$G" --stride "$WORKERS" --offset "$w" \
      --cadence-step "$CADENCE" > "/tmp/${RUN}_w${w}.log" 2>&1 &
    pids="$pids $!"; w=$(( w + 1 ))
  done
  for p in $pids; do wait "$p"; done
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
    --mock "$MOCK" --run "$RUN" --galaxy "$G" --cadence-step "$CADENCE" --merge
  echo "    done $(date '+%F %T'); $(free_gb "$RUNS_DIR") GB free"
done
echo ""
echo "=== all searches done $(date '+%F %T') ==="
