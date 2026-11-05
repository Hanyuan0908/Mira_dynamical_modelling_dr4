#!/bin/bash
# Full dynamical-model pipeline for the 13 non-interacting barred Auriga galaxies, at a 12.5%
# distance error.   Run it with:  sh run_auriga_posteriors.sh
#
# Same method and same 257-point focused Sobol design as the Test 2 N-body campaign; the only thing
# that changes is the galaxy. Per galaxy it does:
#
#   1. trial-potential family       potential_family.py --galaxy auN   (idempotent, seconds-minutes)
#   2. Sobol search, WORKERS wide   test2b_qmc_distance_error_run.py   (the expensive part)
#   3. merge                        ... --merge                        -> results/<run>/evals.csv
#   4. GP posterior + corner        test1_gp_emulator.py, plot_test1_logl_corner.py
#
# The ground-truth potentials (mock/auN_potential/{dm,stars,gas,bh}.ini) and the mocks
# (mock/auN_mock_selection_sigmaDM0.271_N50k.hdf5) must already exist. They do, for all 13.
#
# Every design point is claimed with an atomic O_EXCL lock file and anything already finished is
# skipped, so this is safe to run on two machines at once against the same shared tree, safe to
# re-run, and safe to interrupt.
#
#   sh run_auriga_posteriors.sh --check                   # preflight only: verify and exit
#   sh run_auriga_posteriors.sh                           # all 13, in the table order
#   ONLY="au5 au20" sh run_auriga_posteriors.sh           # just these, in this order
#   WORKERS=11 THREADS=8 sh run_auriga_posteriors.sh      # force the worker geometry
#   PYTHON=/path/to/python sh run_auriga_posteriors.sh
#
# ---------------------------------------------------------------------------------------------
# DISK.  Each run writes ~6.3 GB of orbit library to $MIRA_DR4_RUNS (default <repo>/runs) at the
# default 50 Myr cadence, and /data/hz420 has ~55 GB free. Batch 1 (7 galaxies, ~44 GB) fits; all
# 13 (~82 GB) does not. So this script REFUSES to start
# a new galaxy unless FREE_MIN_GB is free, and waits for the archiver to drain finished runs to the
# portable disk. Two ways to drain:
#
#   ARCHIVE=1   the script archives each run itself as soon as it finishes. Only works where the
#               portable disk is mounted, i.e. on the workstation. NOT on the cluster.
#   otherwise   run scripts/archive_auriga_runs.sh on the workstation in parallel; it watches for
#               finished runs and moves them off. This is the cluster case.
#
# With neither, the campaign stops after ~4 galaxies with the disk full and waits. It will not
# crash and it will not lose work -- it blocks, prints why, and resumes when space appears.
# ---------------------------------------------------------------------------------------------
#
# Outputs: $MIRA_DR4_RUNS/<run>/{traj_NNN.npy, fourier_NNN.npz}
#          results/<run>/{evals.csv, run_meta.json, emulator_results_*.json}
#          figures/<run>_emulator_{corner,slices,loo}_*.png
set -u
CHECK=0
for a in "$@"; do
  case "$a" in
    --check) CHECK=1 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done
cd "$(dirname "$0")"
ROOT=$(cd .. && pwd)
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
CORES=$(nproc)
THREADS=${THREADS:-8}
WORKERS=${WORKERS:-$(( CORES / THREADS ))}
[ "$WORKERS" -lt 1 ] && WORKERS=1
SCHEME=fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024
SIGMA=${SIGMA:-0.271}                   # sigma_DM in the mock file names: 0.271 = 12.5%
TAG=${TAG:-12p5pct}
# 50 Myr orbit sampling (keep every 2nd of the 41 snapshots over 1 Gyr), which is what every
# N-body run added since the option existed uses -- test2_t7_sf_12p5pct and the t2.375/t2.625/
# t3.375/t3.625/t4.5/t5/t6.25 set all record cadence_step=2, snapshot_dt_myr=50. Only the older
# runs predate it and sit at 25 Myr. Halves the trajectory bytes (~6.3 GB a run, not ~12.3) for a
# measured posterior cost below 0.2 km/s/kpc (scripts/test_time_cadence.py). run_test2_posteriors.sh
# still defaults to 1 because changing it would alter how its older runs resume.
CADENCE=${CADENCE:-2}
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}
FREE_MIN_GB=${FREE_MIN_GB:-20}          # refuse to start a galaxy below this
ARCHIVE=${ARCHIVE:-0}                   # 1 = archive each run here as soon as it finishes
# EMULATE=0 skips the GP step. The search only ever produces `observed`/`asis`/`sf` on scheme A,
# which are controls, not the answer: the result column is scheme B + sfratio + K=10, and that
# only exists after rescore_empirical_selection.py. On the cluster, leave the emulator to the
# post-processing step rather than making 3 posteriors nobody should quote.
EMULATE=${EMULATE:-1}
WAIT_MIN=${WAIT_MIN:-10}                # minutes between disk re-checks while blocked
MAX_WAIT_H=${MAX_WAIT_H:-24}            # give up waiting after this long

# the 13 non-interacting barred Auriga haloes (scripts/galaxies.py:AURIGA)
HALOS=${HALOS:-"5 9 10 12 14 17 18 20 22 24 26 27 28"}
ONLY=${ONLY:-}

echo "=== Auriga dynamical-model campaign, sigma_DM $SIGMA (tag $TAG) ==="
echo "    $CORES cores -> $WORKERS workers x $THREADS threads"
echo "    runs dir  $RUNS_DIR"
echo "    archive   $([ "$ARCHIVE" = 1 ] && echo 'inline (portable disk must be mounted)' \
                      || echo 'external (run archive_auriga_runs.sh on the workstation)')"

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

# ---- preflight ------------------------------------------------------------------------------
"$PY" -c "import agama, sklearn, healpy, h5py; print('    environment ok')" || {
  echo "ERROR: the astro312 environment is not usable here; set PYTHON=..." >&2; exit 1; }

# The pipeline modules special-case the galaxy name "au18"; every other Auriga halo silently falls
# back to the N-body recipe and asks for dm_spherical.ini / bulge / disc, which Auriga does not
# have. Fail loudly here rather than after launching WORKERS processes that all die at once.
"$PY" - <<'PRE' || exit 1
import sys
import potential_family as pf
bad = [g for g in ("au5", "au9") if pf._family(g)["file_prefix"] != f"Au{g[2:]}"]
if bad:
    sys.stderr.write(
        "ERROR: potential_family._family() does not know the Auriga halos yet -- "
        f"{bad[0]} resolves to the '{pf._family(bad[0])['file_prefix']}' recipe.\n"
        "       potential_family.py (_family, _expansion, build_files) and\n"
        "       steady_state_likelihood.bar_aligned_potential() hardcode galaxy == 'au18'.\n"
        "       Generalise them to any au<N> before running this campaign.\n")
    sys.exit(1)
print("    pipeline modules know the Auriga halos")
PRE

missing=""
for H in $HALOS; do
  [ -f "$ROOT/mock/au${H}_mock_selection_sigmaDM${SIGMA}_N50k.hdf5" ] || missing="$missing au${H}:mock"
  [ -f "$ROOT/mock/au${H}_potential/summary.json" ]                   || missing="$missing au${H}:potential"
  [ -f "$ROOT/results/au${H}_ground_truth.json" ]                     || missing="$missing au${H}:truth"
done
[ -z "$missing" ] || { echo "ERROR: missing inputs:$missing" >&2; exit 1; }
echo "    mocks, potentials and ground-truth files present for: $HALOS"

# The frame convention is the one thing a silently-wrong galaxy would not crash on, so assert it
# rather than trust it: every mock must rotate clockwise seen from +z (median L_z < 0) with the bar
# at BAR_ANGLE_DEG, exactly as the N-body mocks do, and the pattern speed handed to agama.orbit
# must come out negative.
"$PY" - <<PRE || exit 1
import sys, numpy as np, h5py
from galaxies import galaxy
import potential_family as pf
bad = []
for h in "$HALOS".split():
    g = galaxy(f"au{h}")
    with h5py.File(g["stars"], "r") as f:
        n = len(f["x"])
        s = np.sort(np.random.default_rng(1).choice(n, min(n, 50_000), replace=False))
        x, y, vx, vy = (f[k][:][s] for k in ("x", "y", "vx", "vy"))
        ba = float(f.attrs["bar_angle_deg"])
    R = np.hypot(x, y)
    m = (R > 2) & (R < 10)
    Lz = float(np.median((x * vy - y * vx)[m]))
    om = pf.omega_true(f"au{h}")
    if Lz >= 0 or abs(ba - 27.0) > 1e-6 or om >= 0:
        bad.append(f"au{h}(L_z={Lz:.0f}, bar={ba}, Omega={om:.2f})")
if bad:
    sys.stderr.write("ERROR: frame convention wrong for: " + ", ".join(bad) + "\n")
    sys.exit(1)
print("    frame convention ok: clockwise (L_z<0), bar at 27 deg, Omega_b negative, as the N-body")
PRE

if [ "$CHECK" = 1 ]; then
  echo ""
  echo "    $(free_gb "$RUNS_DIR") GB free on $RUNS_DIR; 257 design points per galaxy;"
  echo "    cadence step $CADENCE ($(( 25 * CADENCE )) Myr sampling)"
  echo "    $(printf '%s' "$HALOS" | wc -w) galaxies x ~$(( 63 * CADENCE / 2 ))e-1 GB = ~$(( $(printf '%s' "$HALOS" | wc -w) * 63 * CADENCE / 20 )) GB of orbit library"
  echo "=== preflight passed; nothing launched (--check) ==="
  exit 0
fi

wait_for_space () {
  waited=0
  while [ "$(free_gb "$RUNS_DIR")" -lt "$FREE_MIN_GB" ]; do
    if [ "$waited" -ge $(( MAX_WAIT_H * 60 )) ]; then
      echo "    GIVING UP: still under ${FREE_MIN_GB} GB free after ${MAX_WAIT_H} h." >&2
      echo "    Archive some finished runs and re-run this script; nothing is lost." >&2
      exit 2
    fi
    echo "    WAITING: $(free_gb "$RUNS_DIR") GB free, need $FREE_MIN_GB."
    echo "             Drain finished runs with:  sh archive_auriga_runs.sh --delete"
    echo "             (re-checking every ${WAIT_MIN} min; waited ${waited} min so far)"
    sleep $(( WAIT_MIN * 60 ))
    waited=$(( waited + WAIT_MIN ))
  done
}

queue () {
  if [ -z "$ONLY" ]; then printf '%s\n' $HALOS; return; fi
  for want in $ONLY; do printf '%s\n' "${want#au}"; done
}

# ---- the campaign ----------------------------------------------------------------------------
for H in $(queue); do
  GAL="au${H}"
  RUN="test2_${GAL}_sf_${TAG}"
  MOCK="${GAL}_mock_selection_sigmaDM${SIGMA}_N50k.hdf5"
  echo ""
  echo "########## $RUN  ($MOCK, galaxy $GAL)  $(date '+%F %T') ##########"

  n_done=$(find "$RUNS_DIR/$RUN" -name 'fourier_*.npz' 2>/dev/null | wc -l)
  if [ "$n_done" -ge 257 ] && [ -f "$ROOT/results/$RUN/evals.csv" ]; then
    echo "    already complete (257/257 points); skipping the search"
  else
    wait_for_space
    echo "    $(free_gb "$RUNS_DIR") GB free; $n_done/257 points already done"

    # trial-potential family for this halo (idempotent; seconds once built)
    if ! OMP_NUM_THREADS=$CORES "$PY" -c "
import steady_state_likelihood as ssl
ssl.bar_aligned_potential(1.0, 1.0, '$GAL')" >/dev/null 2>&1; then
      echo "    building the trial-potential family for $GAL ..."
      OMP_NUM_THREADS=$CORES "$PY" -W ignore potential_family.py --galaxy "$GAL" \
        > "/tmp/${RUN}_potfam.log" 2>&1 \
        || { echo "    ERROR: potential_family.py failed (see /tmp/${RUN}_potfam.log); skipping $RUN"; continue; }
      OMP_NUM_THREADS=$CORES "$PY" -c "
import steady_state_likelihood as ssl
ssl.bar_aligned_potential(1.0, 1.0, '$GAL')" >/dev/null 2>&1 \
        || { echo "    ERROR: trial potential for $GAL still unusable; skipping $RUN"; continue; }
      echo "    trial-potential family ready"
    fi

    pids=""; w=0
    while [ "$w" -lt "$WORKERS" ]; do
      OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
        --mock "$MOCK" --run "$RUN" --galaxy "$GAL" --stride "$WORKERS" --offset "$w" \
        --cadence-step "$CADENCE" > "/tmp/${RUN}_w${w}.log" 2>&1 &
      pids="$pids $!"; w=$(( w + 1 ))
    done
    for pid in $pids; do wait "$pid"; done
    OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
      --mock "$MOCK" --run "$RUN" --galaxy "$GAL" --cadence-step "$CADENCE" --merge
    echo "--- $RUN search done $(date '+%F %T') ---"
  fi

  for V in $([ "$EMULATE" = 1 ] && echo sf); do
    OMP_NUM_THREADS=$THREADS "$PY" -W ignore test1_gp_emulator.py \
      --run "$RUN" --column "logLs_${SCHEME}_${V}" --loo-top 32 \
      > "/tmp/${RUN}_emu_${V}.log" 2>&1 \
      && echo "    posterior ok: $RUN / $V" \
      || echo "    posterior FAILED: $RUN / $V (see /tmp/${RUN}_emu_${V}.log)"
    OMP_NUM_THREADS=$THREADS "$PY" -W ignore plot_test1_logl_corner.py \
      --run "$RUN" --column "logLs_${SCHEME}_${V}" \
      --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2 \
      > "/tmp/${RUN}_corner_${V}.log" 2>&1 || true
  done

  if [ "$ARCHIVE" = 1 ]; then
    echo "--- archiving $RUN to the portable disk ---"
    sh ./collect_runs_to_disk.sh --delete "$RUN" || echo "    WARNING: archiving $RUN failed; it stays here"
  fi
  echo "    $(free_gb "$RUNS_DIR") GB free after $RUN"
done

echo ""
echo "=== all done $(date '+%F %T'); $(free_gb "$RUNS_DIR") GB free on the working disk ==="
