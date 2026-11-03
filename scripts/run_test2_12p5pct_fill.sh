#!/bin/bash
# Test 2 at 12.5%: the six FILL-IN snapshots, to turn five points on the 1:1 line into eleven.
#
#   sh run_test2_12p5pct_fill.sh                  # the whole queue, in Omega_b order
#   ONLY="test2_t1.5_sf_12p5pct" sh ...           # one run
#   WORKERS=12 THREADS=8 sh ...                   # cluster sizing
#
# Seven points inside 2 <= t <= 8 (user, 2026-11-03): two in each interval between the four fastest
# runs already done (t = 2, 3, 4, 5.5) and one between t = 5.5 and 7, placed by even division of each
# interval in snapshot time.
#
#     t      Omega_b   A2      interval
#     2.375    40.87  0.416    t = 2 .. 3
#     2.625    40.21  0.442    t = 2 .. 3
#     3.375    36.83  0.459    t = 3 .. 4
#     3.625    34.79  0.450    t = 3 .. 4
#     4.5      32.07  0.470    t = 4 .. 5.5
#     5.0      30.60  0.461    t = 4 .. 5.5
#     6.25     26.58  0.433    t = 5.5 .. 7
#
# With the five existing runs: twelve points from 25.2 to 43.9, largest gap down from 5.4 to 3.0.
# The range is unchanged - this fills, it does not extend.
#
#
# All inside 2 <= t <= 8 (user, 2026-11-03): the range is not extended, it is filled. With the five
# existing runs: 25.2, 26.0, 28.2, 28.2, 30.6, 32.1, 33.3, 35.2, 38.7, 41.7, 43.9 - eleven points
# over the same factor 1.74, largest gap down from 5.4 to 3.5 km/s/kpc.
#
# The t = 5.5 / t = 6.0 pair is the point of the design: near-identical pattern speed with a 3.6%
# difference in bar strength, which is the only handle in this run on whether the recovery bias
# tracks Omega_b or tracks the bar. Every other pair has the two confounded, because these snapshots
# are one secular sequence of a single galaxy.
set -u
cd "$(dirname "$0")"
export SIGMA=0.271
export TAG=12p5pct
export ONLY=${ONLY:-"test2_t2.375_sf_${TAG} test2_t2.625_sf_${TAG} test2_t3.375_sf_${TAG} \
test2_t3.625_sf_${TAG} test2_t4.5_sf_${TAG} test2_t5_sf_${TAG} test2_t6.25_sf_${TAG}"}
echo "fill-in run: $(nproc) cores; runs: $ONLY"
exec bash ./run_test2_posteriors.sh
