#!/bin/bash
# Post-processing of the focused global search (test1_nbody_qmc_focused_run.py):
#   1. Fourier coefficients on the R-|z| grid and the chi2 / tempered logL table (test1_nbody_qmc_focused_fourier.py)
#   2. GP-emulator posteriors (leave-one-out on the 32 best points) and log L corners of the evaluated points for:
#      the proposed Fourier scheme, the same without the sigma ring averages, and the current sector likelihood.
# Figures: figures/test1_nbody_qmc_focused_{emulator_corner,emulator_loo,emulator_slices,logl_corner}_<column>.png
set -e
cd "$(dirname "$0")"
PY=/Users/hanyuan/miniforge3/envs/astro312/bin/python
$PY -W ignore test1_nbody_qmc_focused_fourier.py
for c in logLs_fourier_R0to6_vR_m24_vphi_sigR_sigphi_m024 logLs_fourier_R0to6_vR_sigR_sigphi_m24_vphi_m024 \
         logLs_sectors_meansdisp; do
  $PY -W ignore test1_gp_emulator.py --run nbody_qmc_focused --column $c --loo-top 32
  $PY -W ignore plot_test1_logl_corner.py --run nbody_qmc_focused --column $c \
      --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2
done
$PY -W ignore plot_test1_logl_corner.py --run nbody_qmc_focused --column logLs_fourier_R0to8_vR_m24_vphi_sigR_sigphi_m024 \
    --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2
echo "postprocess done"
