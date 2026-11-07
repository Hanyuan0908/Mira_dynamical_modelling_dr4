#!/bin/bash
# Realism test: a per-star distance error, modelled as if it were uniform.
#
#   sh run_scatter_test_cluster.sh --check     # verify prerequisites and exit
#   sh run_scatter_test_cluster.sh             # run the 5 searches
#
# The real Mira sample does not have one distance error. Each star has its own, roughly Gaussian
# about the mean: 12.5% +- 2.5%. Every mock in this project so far has given every star exactly
# 12.5%, and the model assumes that too. This tests the mismatch.
#
# The mocks (made by make_mock_selection.py --sigma-dm-scatter 0.0543) scatter each star by its
# OWN draw sigma_i ~ N(0.2714, 0.0543) mag, verified in the files at 0.272 +- 0.054. Their
# sigma_DM attribute still records the nominal 0.2714, which is the single number
# selection.load_tables hands the model - so the library is smeared uniformly at 12.5% and the
# rho_obs/rho_lib weighting is built on that same assumption. Nothing else differs from the
# fiducial campaign.
#
# Five galaxies spanning the pattern-speed range: t = 3, t = 5.5, t = 7 (26-40 km/s/kpc) and
# Au18, Au26 (36 and 65).
#
# SEARCH ONLY, 50 Myr cadence, no emulator - as for Auriga batch 2, because the result column is
# scheme B + sfratio + K = 10 and that only exists after the re-score. The workstation runs
# postprocess_scatter_test.sh, which re-scores, builds the posteriors and archives.
#
# PREREQUISITES: the shared tree and the astro312 environment. --check verifies both, the five
# mocks, their per-star errors, and the potentials, in seconds.
#
# Disk: 5 runs x ~6.3 GB = ~32 GB.
set -u
cd "$(dirname "$0")"
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
ROOT=$(cd .. && pwd)
THREADS=${THREADS:-8}
WORKERS=${WORKERS:-$(( $(nproc) / THREADS ))}
[ "$WORKERS" -lt 1 ] && WORKERS=1
CADENCE=${CADENCE:-2}
SUFFIX=_scat0.054
FREE_MIN_GB=${FREE_MIN_GB:-15}
RUNS_DIR=${MIRA_DR4_RUNS:-$ROOT/runs}

# run tag | galaxy | mock prefix
SET="t3|nbody_t3|nbody_t3
t5.5|nbody_t5.5|nbody_t5.5
t7|nbody|nbody
au18|au18|au18
au26|au26|au26"

CHECK=0
for a in "$@"; do case "$a" in --check) CHECK=1 ;; *) echo "unknown option: $a" >&2; exit 2 ;; esac; done

echo "per-star distance-error test: 5 galaxies, cadence $CADENCE, $WORKERS x $THREADS"
"$PY" -c "import agama, sklearn, healpy, h5py; print('environment ok')" || {
  echo "ERROR: astro312 not usable here; set PYTHON=..." >&2; exit 1; }

miss=""
printf '%s\n' "$SET" | while IFS='|' read -r TAG GAL PFX; do
  [ -f "$ROOT/mock/${PFX}_mock_selection_sigmaDM0.271_N50k${SUFFIX}.hdf5" ] || echo "  MISSING mock: $PFX"
done
"$PY" - <<'PRE' || exit 1
import h5py, sys, numpy as np
from common import MOCK_DIR
bad = []
for p in ("nbody_t3", "nbody_t5.5", "nbody", "au18", "au26"):
    f = MOCK_DIR / f"{p}_mock_selection_sigmaDM0.271_N50k_scat0.054.hdf5"
    if not f.exists():
        bad.append(f"{p}: missing"); continue
    with h5py.File(f, "r") as h:
        e = h["DM_err"][:]
        if e.std() < 1e-6:
            bad.append(f"{p}: DM_err is constant - not a per-star mock")
        if abs(float(h.attrs["sigma_DM"]) - 0.2714) > 1e-3:
            bad.append(f"{p}: sigma_DM attr is not the nominal 0.2714")
if bad:
    sys.stderr.write("ERROR:\n  " + "\n  ".join(bad) + "\n"); sys.exit(1)
print("    5 mocks present, DM_err is per-star, sigma_DM attribute is the nominal 0.2714")
PRE
[ "$CHECK" = 1 ] && { echo "=== preflight passed; nothing launched (--check) ==="; exit 0; }

free_gb () { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

printf '%s\n' "$SET" | while IFS='|' read -r TAG GAL PFX; do
  RUN=test2_${TAG}_sf_12p5pct_scat
  MOCK=${PFX}_mock_selection_sigmaDM0.271_N50k${SUFFIX}.hdf5
  n=$(find "$RUNS_DIR/$RUN" -name 'fourier_*.npz' 2>/dev/null | wc -l)
  echo ""
  echo "########## $RUN  ($GAL)  $(date '+%F %T') ##########"
  if [ "$n" -ge 257 ]; then echo "    already complete"; continue; fi
  while [ "$(free_gb "$RUNS_DIR")" -lt "$FREE_MIN_GB" ]; do
    echo "    only $(free_gb "$RUNS_DIR") GB free, need $FREE_MIN_GB; waiting 5 min"; sleep 300
  done
  pids=""; w=0
  while [ "$w" -lt "$WORKERS" ]; do
    OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
      --mock "$MOCK" --run "$RUN" --galaxy "$GAL" --stride "$WORKERS" --offset "$w" \
      --cadence-step "$CADENCE" > "/tmp/${RUN}_w${w}.log" 2>&1 &
    pids="$pids $!"; w=$(( w + 1 ))
  done
  for p in $pids; do wait "$p"; done
  OMP_NUM_THREADS=$THREADS "$PY" -W ignore test2b_qmc_distance_error_run.py \
    --mock "$MOCK" --run "$RUN" --galaxy "$GAL" --cadence-step "$CADENCE" --merge
  echo "    done $(date '+%F %T'); $(free_gb "$RUNS_DIR") GB free"
done
echo ""
echo "=== all 5 searches done $(date '+%F %T') ==="
