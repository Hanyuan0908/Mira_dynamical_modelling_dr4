#!/bin/bash
# Orbital frequencies in the bar frame, for the dynamical bar length. All 25 galaxies.
#
#   sh run_bar_dynamics_cluster.sh --check     # verify prerequisites and exit (seconds)
#   sh run_bar_dynamics_cluster.sh             # truths first, then the posterior draws
#
# WHY A RE-INTEGRATION. The saved orbit library cannot supply this. Its 25/50 Myr cadence puts
# Nyquist (20 / 10 per Gyr) BELOW the median Omega_R of the x1 orbits in every galaxy, slow bars
# included, so 73-99% of the radial frequencies alias and Omega_R/Omega_x comes out wrong by
# 0.4-1.2 against the +-0.2 tolerance of the Zhang et al. (2024) cut. The BASELINE is fine: 1 Gyr
# at 1 Myr keeps 83-93% of the reference x1 set at +-0.06-0.09. See diagnose_x1_frequency_sampling.py
# and results/x1_frequency_sampling.csv. So: same 1 Gyr, same integrator, 25x finer output, and
# only the frequencies are kept.
#
# WHAT IT RUNS, per galaxy:
#   truth   the family at (bar_frac, S) = (1, 1) with the tabulated Omega_b - identically the
#           galaxy's own potential - integrating ORIGINAL SNAPSHOT PARTICLES: no distance error,
#           no selection function. This is the true dynamical length to compare against.
#   draws   N_DRAW potentials from the posterior RE-WEIGHTED BY THE KINEMATIC phi_b PRIOR, each
#           integrating the run's own mock stars as observed (distance errors and selection in).
#
# NO TRAJECTORIES ARE NEEDED OR WRITTEN. The inputs are mock/ and results/ only - this cluster does
# not need the 83 GB of runs/, and the output is ~0.5 MB per draw (results/bar_dynamics/).
#
# Cost, measured: 3,000 orbits in 11 s on 8 threads, so ~41 s per draw at 10k stars and ~3.4 min
# per truth at 50k. 25 galaxies x (30 draws + 1 truth) is ~10 h of 8-thread work, i.e. ~1.3 h on 64
# cores. Raise N_DRAW for a tighter band; it scales linearly.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
ROOT=$(cd .. && pwd)
THREADS=${THREADS:-8}
WORKERS=${WORKERS:-$(( $(nproc) / THREADS ))}
[ "$WORKERS" -lt 1 ] && WORKERS=1
N_DRAW=${N_DRAW:-30}
N_STAR=${N_STAR:-10000}
N_STAR_TRUTH=${N_STAR_TRUTH:-50000}
TAGS=${TAGS:-"t2 t2.375 t2.625 t3 t3.375 t3.625 t4 t4.5 t5 t5.5 t6.25 t7 au5 au9 au10 au12 au14 au17 au18 au20 au22 au24 au26 au27 au28"}
OUT=$ROOT/results/bar_dynamics

CHECK=0
for a in "$@"; do case "$a" in --check) CHECK=1 ;; *) echo "unknown option: $a" >&2; exit 2 ;; esac; done

N=$(printf '%s' "$TAGS" | wc -w)
echo "dynamical bar length: $N galaxies, $N_DRAW draws x $N_STAR stars + truth x $N_STAR_TRUTH,"
echo "  1 Gyr at 1 Myr, $WORKERS workers x $THREADS threads"
"$PY" -c "import agama, sklearn, emcee, h5py; print('environment ok')" || {
  echo "ERROR: astro312 not usable here; set PYTHON=..." >&2; exit 1; }

TAGS="$TAGS" "$PY" - <<'PRE' || exit 1
import json, os, sys
from common import RESULTS_DIR, MOCK_DIR
from galaxies import galaxy
from scan_bar_angle_windows import CACHE
B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
bad = []
if not CACHE.exists():
    bad.append(f"ring cache missing: {CACHE}")
for t in os.environ["TAGS"].split():
    r = f"test2_{t}_sf_12p5pct"
    d = RESULTS_DIR / r
    if not (d / "run_meta.json").exists():
        bad.append(f"{t}: no run_meta.json"); continue
    m = json.load(open(d / "run_meta.json"))
    if not (d / f"emulator_results_logLs_{B}_sfratio.json").exists():
        bad.append(f"{t}: no scheme B + sfratio posterior")
    if not (d / "evals.csv").exists():
        bad.append(f"{t}: no evals.csv (the emulator needs it)")
    if not (MOCK_DIR / m["mock"]).exists():
        bad.append(f"{t}: mock {m['mock']} missing")
    if not galaxy(m["galaxy"])["stars"].exists():
        bad.append(f"{t}: snapshot particles {galaxy(m['galaxy'])['stars'].name} missing "
                   "(needed for the TRUE dynamical length)")
    # CALL the truth function, do not just look for its file. The N-body snapshots have no
    # <gal>_ground_truth.json - their pattern speed is in nbody_pattern_speed_evolution.csv - and
    # checking for files instead of calling let 11 of 25 truth integrations die silently on
    # 2026-11-07. The same trap is recorded for 2026-11-05.
    try:
        from test2b_qmc_distance_error_run import truth_omega
        truth_omega(m["galaxy"])
    except Exception as e:
        bad.append(f"{t}: no truth pattern speed for {m['galaxy']} ({type(e).__name__})")
if bad:
    sys.stderr.write("ERROR:\n  " + "\n  ".join(bad) + "\n"); sys.exit(1)
print(f"    {len(os.environ['TAGS'].split())} galaxies: posterior, evals, mock and snapshot "
      "particles all present")
PRE
[ "$CHECK" = 1 ] && { echo "=== preflight passed; nothing launched (--check) ==="; exit 0; }
mkdir -p "$OUT"

# ---- pass 1: the true dynamical length, one process per galaxy, WORKERS at a time ------------
echo ""
echo "########## pass 1/2: truth potentials, snapshot particles  $(date '+%F %T') ##########"
k=0
for T in $TAGS; do
  GAL=$("$PY" -c "import json;print(json.load(open('$ROOT/results/test2_${T}_sf_12p5pct/run_meta.json'))['galaxy'])")
  if [ -f "$OUT/truth_${GAL}.npz" ]; then echo "    $GAL already done"; continue; fi
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore integrate_orbit_frequencies.py --truth \
    --galaxy "$GAL" --n-star "$N_STAR_TRUTH" > "/tmp/bardyn_truth_${GAL}.log" 2>&1 &
  k=$(( k + 1 ))
  [ $(( k % WORKERS )) -eq 0 ] && wait
done
wait
echo "    truths: $(ls "$OUT"/truth_*.npz 2>/dev/null | wc -l)/$N done  $(date '+%F %T')"

# ---- pass 2: posterior draws, one galaxy at a time, draws split across workers ---------------
echo ""
echo "########## pass 2/2: posterior draws (kinematic phi_b prior)  $(date '+%F %T') ##########"
for T in $TAGS; do
  RUN=test2_${T}_sf_12p5pct
  n=$(ls "$OUT/$RUN"/draw_*.npz 2>/dev/null | wc -l)
  echo ""
  echo "--- $RUN  ($n/$N_DRAW already)  $(date '+%F %T')"
  if [ "$n" -ge "$N_DRAW" ]; then echo "    complete"; continue; fi
  pids=""; w=0
  while [ "$w" -lt "$WORKERS" ]; do
    OMP_NUM_THREADS=$THREADS "$PY" -W ignore integrate_orbit_frequencies.py --run "$RUN" \
      --n-draw "$N_DRAW" --n-star "$N_STAR" --stride "$WORKERS" --offset "$w" \
      > "/tmp/${RUN}_bardyn_w${w}.log" 2>&1 &
    pids="$pids $!"; w=$(( w + 1 ))
  done
  for p in $pids; do wait "$p"; done
  echo "    $(ls "$OUT/$RUN"/draw_*.npz 2>/dev/null | wc -l)/$N_DRAW draws  $(date '+%F %T')"
done

echo ""
echo "=== all done $(date '+%F %T'):  $(ls "$OUT"/truth_*.npz 2>/dev/null | wc -l)/$N truths, "
echo "    $(find "$OUT" -name 'draw_*.npz' | wc -l)/$(( N * N_DRAW )) draws in results/bar_dynamics ==="
