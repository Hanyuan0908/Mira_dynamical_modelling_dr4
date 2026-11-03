#!/bin/bash
# Test 2 at a 12.5% DISTANCE ERROR -- the workstation half.   Run it with:  sh run_test2_12p5pct_local.sh
#
# Companion to run_test2_12p5pct_cluster.sh; see that file for what the campaign is and how the split
# works. This half does the two t = 7 runs, in this order:
#   test2_t7_sf_12p5pct   (with the selection function -- the one that matters)
#   test2_t7_nosf_12p5pct (the no-selection benchmark)
# ONLY sets the order as well as the queue, so t7_sf is worked first.
#
# Override the queue the same way:
#   ONLY="" sh run_test2_12p5pct_local.sh        # work the entire six-run queue
set -u
cd "$(dirname "$0")"
export SIGMA=0.271
export TAG=12p5pct
export ONLY=${ONLY:-"test2_t7_sf_${TAG} test2_t7_nosf_${TAG}"}
echo "workstation run: $(nproc) cores; runs: $ONLY"
exec bash ./run_test2_posteriors.sh
