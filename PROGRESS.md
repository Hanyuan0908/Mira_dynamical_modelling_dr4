# Progress log

Mock tests of the steady-state bar-potential method for Gaia DR4 Miras. Each entry lists what
was done, the scripts that reproduce it (run from the repo root with the `astro312` env), the
figures, and the key numbers.

## Conventions

- **Frame:** astropy-default Galactocentric axes. The Sun is at x = −R0, y points to l = 90°,
  z points to the NGP. The disc rotates clockwise seen from +z (L_z < 0).
  v_φ is quoted positive along rotation, i.e. v_φ = −(x v_y − y v_x)/R.
- **Sun:** R0 = 8.122 kpc, z☉ = 20.8 pc, v☉ = (12.9, 245.6, 7.78) km/s, as set in
  `scripts/common.py` (the DR3 Mira catalogue uses the same values).
- **Bar:** MW-like angle of 27° to the Sun–GC line, near end at l > 0 (y > 0).
- **Galactic coordinates:** AGAMA `getGalacticFromGalactocentric` and its inverse
  (`scripts/observe.py`). Proper motions are in mas/yr.

## Data products (not in git)

| File | Made by | Content |
|---|---|---|
| `data/gaia_dr3_lpv_mira.csv` | copied from `~/Desktop/PhD_projects/spectroscopic_catalogues/` | 55,074 DR3 Miras with PL distances, x/y/z, periods, ages |
| `mock/halo_18/{snapdir,groups}_127/` | Auriga public release (Globus) | Au18 `Original/4` z = 0 snapshot and Subfind catalogue |
| `mock/au18_stars_mwframe.hdf5` | `scripts/prepare_au18_mwframe.py` | 2.05M Au18 star particles within r < 50 kpc in the MW-like frame; transformation stored in the attributes |
| `mock/mock_perfect_R8_N50k.hdf5` | `scripts/make_mock_perfect.py` | Test 1 mock: 50k stars, R < 8 kpc, true 6D + error-free observables |

---

## 2026-10-26 — DR3 Mira overview

- Scripts: `scripts/plot_dr3_mira_period_vs_R.py`, `scripts/plot_dr3_mira_spatial_period_bins.py`
- Figures: [P(period | R)](figures/dr3_mira_period_vs_R.png),
  [x–y / x–z maps in 5 period bins](figures/dr3_mira_spatial_period_bins.png)
- Results:
  - The median period is ≈ 270 d at R < 1 kpc, rises to ≈ 320–325 d at R ≈ 4–9 kpc,
    then falls slowly outward. The bar/bulge has relatively more short-period (old) Miras.
  - Short-period bins (median age 10.5 Gyr) are thick and round. The long-period bin
    (5.2 Gyr) is thin and disc-like.
  - The selection is visible in the maps. There are streaks along lines of sight from the
    PL-distance errors (median σ_d ≈ 1 kpc), an extinction gap at z ≈ 0 towards the GC, and
    a sparse far side.
  - Only 3,284 of the stars have DR3 radial velocities.

## 2026-10-26 — Au18 preparation (bar placed at the MW angle)

- Scripts: `scripts/auriga.py` (reader and frame utilities), `scripts/prepare_au18_mwframe.py`,
  `scripts/plot_au18_density.py`
- Figures: [m = 2 Fourier profile](figures/au18_bar_fourier.png),
  [face-on / edge-on stellar density](figures/au18_stellar_density_mwframe.png)
- Results:
  - M200 = 1.22×10¹² M☉ and R200 = 225 kpc. The stellar centre agrees with SubhaloPos to
    about 10 pc.
  - Mean v_φ ≈ 245 km/s at R = 8 kpc (|z| < 0.5 kpc), so Au18 is MW-like.
  - Bar: A2_max = 0.36 at R = 1.4 kpc; half-length R_bar = 5.75 kpc (where A2 falls to half
    its peak). The phase is rotated from 40.1° to −27.0° and stays flat out to R_bar.
    The edge-on view shows a clear boxy/peanut bulge.
  - Peak memory was about 1 GB.

## 2026-10-26 — Au18 face-on kinematics

- Script: `scripts/plot_au18_velocity_maps.py`
- Figure: [⟨v_R⟩, ⟨v_φ⟩, ⟨v_φ⟩ − azimuthal mean](figures/au18_velocity_maps_faceon.png)
  (|z| < 1 kpc, mass-weighted, 200 pc pixels)
- Results:
  - ⟨v_R⟩ shows a clean quadrupole of about ±40 km/s inside R ≈ 4 kpc. Its zero lines lie
    along the bar's major and minor axes, with the signs expected for bar-aligned orbits in
    clockwise rotation.
  - The non-axisymmetric part of v_φ is a quadrupole of about ±50 km/s, rotated 45° from the
    v_R one: stars are faster along the minor axis and slower along the major axis.
  - Outside the bar (R ≈ 6–10 kpc) there is a weaker v_R pattern of about ±20–30 km/s,
    probably from spirals or lopsidedness. It is not steady in the bar frame.

## 2026-10-26 — Distance-error test (frame and transform validation)

- Scripts: `scripts/observe.py`, `scripts/plot_au18_vR_distance_error.py`
- Figure: [⟨v_R⟩ true vs with σ_d/d = 0.2](figures/au18_vR_distance_error.png)
- Validation:
  - AGAMA agrees with astropy to < 3″ in position, < 10⁻⁴ mas/yr in proper motion and
    10⁻¹² km/s in v_los.
  - A noise-free round trip returns the input to 6×10⁻¹².
- Results (m = 2 part of ⟨v_R⟩, |z| < 1 kpc, 1 < R < 4 kpc):

  | | amplitude | position angle of the zero-v_R line |
  |---|---|---|
  | true | 34.9 km/s | −24.7° (bar axis −27°) |
  | σ_d/d = 0.2 | 43.6 km/s | −6.6° |

  - Distance errors slide stars along their line of sight. This changes their azimuth, so
    rotation leaks into v_R: v_R,obs ≈ v_φ sin(φ_obs − φ_true).
  - The result is a fake quadrupole aligned with the Sun–GC line, with the same signs as the
    bar's real one. The two add, so the inferred bar angle is pulled towards 0°.
  - **Implication:** the steady-state comparison must forward-model the errors. One option is
    to integrate, re-observe the evolved sample with the same noise and compare in
    observable space.

## 2026-10-26 — Test 1: perfect mock (no errors, no selection)

- Scripts: `scripts/make_mock_perfect.py` (default N = 50,000, R < 8 kpc, seed 2026),
  `scripts/plot_mock_overview.py mock/mock_perfect_R8_N50k.hdf5`
- Figure: [counts x–y / x–z, ⟨v_R⟩, ⟨v_φ⟩](figures/mock_perfect_R8_N50k_overview.png)
- Mock:
  - 50,000 of the 1.48M Au18 star particles with cylindrical R < 8 kpc, drawn uniformly by
    number. No z cut: 21.8% have |z| > 1 kpc and 4.0% have |z| > 3 kpc.
  - Heliocentric distances run from 0.39 to 48 kpc.
- Results:
  - The bar and boxy/peanut bulge are visible in the counts.
  - The v_R quadrupole survives the down-sampling: amplitude 33.8 ± 1.2 km/s, angle
    −24.0 ± 1.0° (200 bootstrap resamples). The full particle set gives 34.9 km/s and −24.7°.
  - Per-pixel noise is ≈ ±30 km/s in the outer disc (500 pc pixels).
- Status: catalogue ready. **Next:** integrate in the true and trial potentials and define the
  stationarity statistic.

---

## Next steps

1. Measure Au18's bar pattern speed Ω_p, the "true" value the test must recover
   (Tremaine–Weinberg or continuity equation on one snapshot, or adjacent snapshots).
2. Build Au18's own potential: an AGAMA multipole/CylSpline expansion of the dark matter,
   gas and stars, using the stored transformation.
3. Test 1 itself: integrate the perfect mock in the true potential (the null test) and in
   trial potentials with m ≠ 0 rescaled in amplitude, length and Ω_p. Choose the
   discrepancy statistic and the integration time.
4. Add realism step by step: distance errors (forward-modelled), RV jitter, selection
   function, period–age tracer weighting.
