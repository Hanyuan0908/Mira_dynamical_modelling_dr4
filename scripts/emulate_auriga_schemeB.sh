#!/bin/bash
# GP posteriors for the Auriga runs on SCHEME B (sigma ring averages dropped), the scheme the
# N-body reference runs and the paper figure use. The campaign's own emulator calls used scheme A
# by mistake; the scheme-B chi2 is already in each run's evals.csv (scoring computes every scheme),
# so this is only the emulator + corner plot -- no rescoring, no integration.
#
#   sh emulate_auriga_schemeB.sh                 # all 7, waiting for any still scoring
#   HALOS="9 10" sh emulate_auriga_schemeB.sh
set -u
cd "$(dirname "$0")"; ROOT=$(cd .. && pwd)
PY=${PYTHON:-/data/hz420-2/astro312/bin/python}
HALOS=${HALOS:-"9 10 17 18 24 26 27"}
TAG=${TAG:-12p5pct}
B=fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024
VAR=${VAR:-sfratio}

for H in $HALOS; do
  RUN=test2_au${H}_sf_${TAG}
  # wait until that galaxy's rescore has produced all 257 points and folded them in
  while : ; do
    n=$(find "$ROOT/runs/$RUN" -name 'sfratio_*.npz' 2>/dev/null | wc -l)
    has=$("$PY" -c "
import pandas as pd,sys
try: print(int('logLs_${B}_${VAR}' in pd.read_csv('$ROOT/results/$RUN/evals.csv').columns))
except Exception: print(0)")
    [ "$n" -ge 257 ] && [ "$has" = "1" ] && break
    echo "$(date '+%H:%M:%S')  $RUN: $n/257 scored, column present=$has; waiting"
    sleep 120
  done
  echo "$(date '+%H:%M:%S')  === $RUN scheme B / $VAR ==="
  OMP_NUM_THREADS=4 "$PY" -W ignore test1_gp_emulator.py --run "$RUN" \
      --column "logLs_${B}_${VAR}" --loo-top 32 > "/tmp/${RUN}_emuB.log" 2>&1 \
      && echo "    posterior ok" || { echo "    FAILED (see /tmp/${RUN}_emuB.log)"; tail -3 "/tmp/${RUN}_emuB.log"; }
  OMP_NUM_THREADS=4 "$PY" -W ignore plot_test1_logl_corner.py --run "$RUN" \
      --column "logLs_${B}_${VAR}" --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2 \
      > "/tmp/${RUN}_cornerB.log" 2>&1 || true
done
echo "=== scheme B posteriors done $(date '+%F %T') ==="
