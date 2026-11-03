#!/bin/bash
# Test 2 at a 12.5% DISTANCE ERROR -- the cluster half.   Run it with:  sh run_test2_12p5pct_cluster.sh
#
# This is the 10% campaign repeated at 12.5% (sigma_DM = 0.125 * 5 / ln 10 = 0.2714 mag), with the same
# 257-point focused Sobol design, the same DR3-Mira-like selection function and the same chi2 schemes.
# Nothing about the method changes; only the mock's distance error.
#
# The work is split between this cluster and the workstation:
#   cluster      test2_t2_sf_12p5pct, test2_t3_sf_12p5pct,
#                test2_t4_sf_12p5pct, test2_t5.5_sf_12p5pct            (4 x 257 points)
#   workstation  test2_t7_sf_12p5pct, test2_t7_nosf_12p5pct            (2 x 257 points)
#
# The cluster gets three because it is the faster of the two in wall-clock: on the 10% campaign it
# turned a 257-point run round in 2.6-2.7 h against 8.2-8.4 h here. Per core it is actually the
# slower machine (399 s per point against 229 s here); it wins on width, running ~11 workers to our 2.
# The split is only to keep the two machines out of each other's way. Every design point is claimed
# with an atomic O_EXCL lock file and anything already finished anywhere is skipped, so it is safe to
# run this twice, or to take over the whole queue if the workstation falls behind:
#
#   ONLY="" sh run_test2_12p5pct_cluster.sh                 # work the entire six-run queue
#   ONLY="test2_t2_sf_12p5pct" sh run_test2_12p5pct_cluster.sh    # pick up one specific run
#
# Requirements: the shared filesystem (/data/hz420/...) and the astro312 environment -- nothing else.
# The six mocks and all five trial-potential families are already built and committed to the shared
# tree, so there is no set-up step. Outputs land in the shared tree, NOT on any removable disk:
#   /data/hz420/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs/<run>/{traj_NNN.npy,fourier_NNN.npz}
#   .../results/<run>/{evals.csv, run_meta.json, emulator_results_*.json}
# They are collected and moved to the portable disk from the workstation afterwards.
#
# Knobs, if the defaults are wrong here:
#   PYTHON=/path/to/python   THREADS=8   WORKERS=12   sh run_test2_12p5pct_cluster.sh
#
# Expect roughly 235 s of single-core integration per design point, so ~17 h of serial work per run;
# with WORKERS x THREADS cores it divides by WORKERS. Each run writes ~12.3 GB.
set -u
cd "$(dirname "$0")"
export SIGMA=0.271            # sigma_DM in the mock file names -> 12.5% distance error
export TAG=12p5pct            # suffix of the run names
export ONLY=${ONLY:-"test2_t2_sf_${TAG} test2_t3_sf_${TAG} test2_t4_sf_${TAG} test2_t5.5_sf_${TAG}"}
exec bash ./run_test2_posteriors_cluster.sh
