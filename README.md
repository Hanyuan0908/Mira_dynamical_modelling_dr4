# Mira_dynamical_modelling_dr4

Mock tests, ahead of Gaia DR4, of a steady-state test of the Galactic bar potential
with Mira variables. The idea: integrate the present-day 6D Mira sample forward in a
trial bar potential (built by rescaling the m != 0 part of the Portail+17 / Sormani+22
AGAMA model in amplitude, length and pattern speed). Only in the correct potential is
the sample stationary in the bar frame, so its kinematics stay unchanged.

## Layout

- `scripts/` - analysis code (`common.py` holds paths, constants and plot style)
- `figures/` - saved plots
- `data/` - observational catalogues (not tracked), e.g. `gaia_dr3_lpv_mira.csv`
- `mock/` - simulation snapshots and mock catalogues (not tracked), e.g. Auriga Au18

## Environment

Use the `astro312` conda environment, e.g.

```
conda activate astro312
python scripts/plot_dr3_mira_period_vs_R.py
python scripts/plot_dr3_mira_spatial_period_bins.py
python scripts/prepare_au18_mwframe.py      # needs mock/halo_18/{snapdir,groups}_127
python scripts/plot_au18_density.py
python scripts/plot_au18_velocity_maps.py
python scripts/plot_au18_vR_distance_error.py
```

## Mock galaxy

Auriga Au18 (`Original/4`, MHD, level 4, z = 0 snapshot 127) from the Auriga public
release (Grand et al. 2024), downloaded via Globus into `mock/halo_18/`.
`prepare_au18_mwframe.py` centres it and spins the disc to the MW sense (L_z < 0 in the
astropy Galactocentric frame). It then rotates the bar (m=2 phase) to 27 deg from the
Sun-GC line and writes `mock/au18_stars_mwframe.hdf5`. The transformation is stored in
the file attributes so it can be applied to the dark matter and gas.
