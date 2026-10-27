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
| `mock/au18_potential/{dm,stars,gas,bh}.ini` | `scripts/build_au18_potential.py` | Au18 ground-truth potential (AGAMA coefficient files, MW-like frame) |
| `mock/au18_potential/Au18_{full,axi,baryon_full,baryon_axi}.ini` | `scripts/potential_family.py` | ground truth, axisymmetric part, full and axisymmetrised baryons, for the trial-potential family |

Small numeric results (tracked): `results/au18_ground_truth.json` (pattern speed, resonances, bar),
`results/au18_potential_checks.json` (potential benchmarks).

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

## 2026-10-26 — Au18 ground-truth potential and pattern speed

- Scripts: `scripts/build_au18_potential.py`, `scripts/measure_au18_pattern_speed.py`,
  `scripts/check_au18_potential.py`
- Figures:
  [stellar density check](figures/au18_potential_density_check.png),
  [rotation curve and asymmetric drift](figures/au18_rotation_curve_check.png),
  [pattern speed and resonances](figures/au18_pattern_speed.png),
  [AGAMA vs Arepo potential](figures/au18_potential_vs_arepo.png)
- **Potential.** AGAMA, units Msun–kpc–km/s, built in the MW-like frame (bar at −27°).

  | component | expansion | hyperparameters | where they come from |
  |---|---|---|---|
  | DM (7.19M particles, r < 500 kpc; no low-res particles inside) | Multipole, `symmetry='spherical'` | `gridSizeR=25, rmin=0.1, rmax=500` | AGAMA `py/example_gizmo_snapshot.py` (FIRE zoom: rmin 0.1, rmax 400–500, default gridSizeR) |
  | stars (2.05M, wind excluded) | CylSpline, `symmetry='none'` | `mmax=8, gridSizeR=30, gridSizez=30, Rmin=0.15, Rmax=50, zmin=0.1, zmax=50` | your `Nexus/get_potential*.py` (30×25, Rmin 0.15, zmin 0.1, Rmax 50, zmax 20, mmax 8); zmax raised to Rmax (see below) |
  | gas (309k cells inside the grid, 2.3×10¹⁰ M☉) | CylSpline, `symmetry='none'` | same as stars | same |
  | central BH (1.38×10⁸ M☉, 41 pc off-centre) | Plummer | `scaleRadius=0.369` kpc (level-4 softening) | added for completeness; contributes 8.6 km/s to v_c at 8 kpc |

  - CylSpline ignores mass outside its grid, and 2.9×10¹¹ M☉ of hot halo gas lies outside it.
    With zmax = Rmax the grid contains the whole sphere r < 50 kpc, so the excluded gas is
    roughly spherical. Its contribution to v_c² at 8 kpc is −1.3 (km/s)², which is negligible.
  - At R = 8 kpc, v_c = 258.5 km/s: DM 181.5, stars 182.5, gas 21.7, BH 8.6 km/s.
  - Bug fixed: AGAMA's `export()` does not store the parameters of analytic potentials.
    `bh.ini` is now written by hand, and the build asserts that every reloaded file reproduces
    its potential.
- **Benchmarks** (`results/au18_potential_checks.json`):
  1. Stellar density implied by the potential vs the particles:
     - face-on median ratio 0.998;
     - 68% of pixels within 0.03 dex (the scatter is Poisson noise);
     - Σ(R) agrees within −3% / +5% over 0.5–15 kpc;
     - ρ(z) at R = 2, 5 and 8 kpc agrees;
     - bar A2 peak 0.355 vs 0.357, phase −26.8° in both.
  2. Rotation curve: the stellar ⟨v_φ⟩ predicted by the Jeans equation from the potential's
     v_c (using measured σ_R, σ_φ, ν and ∂⟨v_R v_z⟩/∂z, |z| < 0.5 kpc) matches the measured
     ⟨v_φ⟩ with median offset −0.1% and rms 1.0% for R_bar < R < 15 kpc.
     - It also holds down to R ≈ 3 kpc.
     - Asymmetric drift at 8 kpc is 13.5 km/s.
     - Gas rotates 1.6% below v_c, consistent with pressure support.
  3. Extra check: the AGAMA potential at the star particles vs the potential Arepo stored for
     them (up to a constant).
     - At r ≲ 2 kpc AGAMA is deeper by up to 1.5×10⁴ (km/s)². This is Arepo's gravitational
       softening (ε ≈ 0.37 kpc); the AGAMA potential is unsoftened.
     - At 3–20 kpc the residual is within ±1000 (km/s)², against a potential range of 2×10⁵.
     - The face-on residual shows an **m = 2 pattern aligned with the bar** (see caveat 1)
       plus a uniform gradient (external tidal field).
- **Pattern speed** (`results/au18_ground_truth.json`):
  - Method: the Dehnen, Semczuk & Schönrich (2023) Fourier method from
    `python_script/patternSpeed-main`, used as in `get_bar_property_ForSingleSnapshot.py`.
    The disc is mirrored to rotate anticlockwise; input is kpc and km/s, output km/s/kpc.
    NumPy 2 removed `np.asfarray`, so the wrapper adds a one-line shim instead of editing the
    package.
  - **Ω_p = 36.27 ± 0.36 km/s/kpc**, using the automatic bar region 0.45–3.70 kpc.
    Fixed radial ranges give 36.13–36.26. In the MW-like frame the bar rotates clockwise:
    φ_bar(t) = −27° − Ω_p t.
  - Resonances: corotation at R_CR = 7.19 kpc, OLR at 11.67 kpc, ILRs at 0.80 and 2.15 kpc.
    R_CR / R_bar = 1.25, a fast bar like the MW's.
- **Caveats for the mock tests:**
  1. *Spherical DM drops the halo's bar response.* A non-spherical DM Multipole (lmax = 8,
     symmetry none) has an m = 2 potential aligned with the stellar bar, at 10–20% of the
     stellar m = 2 amplitude for R = 1–6 kpc (27% at 8 kpc).
     - The spherical-DM ground truth therefore underestimates the non-axisymmetric force the
       stars actually feel.
     - Option: use a non-spherical DM Multipole as the ground truth, or treat the DM's m = 2
       as part of the "bar" when rescaling.
  2. *Softening.* Inside r ≲ 2 kpc the unsoftened AGAMA potential is deeper than the softened
     one the particles moved in. The innermost region may look non-stationary in the null test.
  3. An older script, `get_bar_property_ForSingleSnapshot.py`, multiplies Ω by 0.98. That
     factor is not applied here, because kpc and km/s inputs already give km/s/kpc.

## 2026-10-26 — Trial-potential family (bar strength, length, pattern speed)

- Scripts: `scripts/potential_family.py` (`make_potential`, `to_bar_frame`; run it to build the
  .ini files), `scripts/test_potential_family.py`
- Figure: [potential family checks](figures/potential_family_check.png);
  numbers in `results/potential_family_checks.json`
- Construction: follows `make_potential()` in `python_script/oCen_bar`.
  - The axisymmetrised baryonic potential is rebuilt as a CylSpline from stars + gas with
    `symmetry='axisymmetric'` on the same grid.
  - It is subtracted through a copy with `scale=[-1, 1]`; AGAMA has no `-` operator for
    potentials.

  ```python
  pot = agama.Potential(pot_axi,                                            # DM + BH + axisym. baryons
          agama.Potential(file=Au18_baryon_full.ini, scale=[A, S], rotation=[[0, phi0], [1, phi0 + Omega_b]]),
          agama.Potential(file=Au18_baryon_axi.ini,  scale=[-A, S]))       # A = bar_frac * S
  ```

  - `scale` gives Φ → (A/S) Φ(x/S). With A = bar_frac·S (as in oCen_bar), bar_frac sets the bar
    amplitude and S its length.
  - `rotation` turns the bar by ψ = φ0 + Ω_b t, anticlockwise for Ω_b > 0.
  - Truth: bar_frac = 1, S = 1, Ω_b = −36.27 km/s/kpc (clockwise in the MW-like frame).
  - Time unit: kpc/(km/s) = 0.978 Gyr.
  - Time-dependent tables (bartime, phi_b, bar_frac(t), S(t)) are supported as in oCen_bar.
  - DM stays spherical for now, so the m ≠ 0 part comes only from stars and gas.
- Checks (all pass):
  - bar_frac = S = 1 reproduces `Au18_full.ini` to 4×10⁻¹⁶ at t = 0, and the rotated truth at
    t > 0.
  - v_c changes by ≤ 3×10⁻⁴ for bar_frac ∈ {0, 0.5, 2} and S ∈ {0.8, 1.25}.
  - bar_frac = 0 is axisymmetric to 6×10⁻¹¹ (km/s)².
  - The non-axisymmetric part equals bar_frac·Φ_nonaxi(x/S) to 2×10⁻¹⁴.
  - The m = 2 phase follows −27° + Ω_b t to 10⁻¹³ deg.
  - Orbits integrated in `make_potential` (inertial frame) and mapped with `to_bar_frame` match
    `agama.orbit(static, Omega=Omega_b)` to 3×10⁻⁵ kpc and 0.008 km/s after 0.49 Gyr.
- Conventions confirmed for later steps:
  - `agama.orbit(..., Omega=...)` takes and returns *inertial* velocities expressed along the
    rotating axes, not velocities relative to the rotating frame.
  - `export()` of a composite loses the parameters of analytic components (the Plummer BH).
    The `Au18_*.ini` files are therefore written by concatenating the component files.

## 2026-10-27 — Orbit-integration cost and fast-potential benchmark

- Scripts: `scripts/time_orbit_integration.py`; `scripts/benchmark_fast_potential.py`
  (5000 random mock stars). `make_static_fast()` was added to `scripts/potential_family.py`.
  It re-expands the barred baryons into one CylSpline with a chosen mmax and symmetry, for
  integration in the bar frame with `Omega=Omega_b`.
- Results: `results/orbit_timing.json`, `results/fast_potential_benchmark.json`.
  Machine: Apple M5 Pro, 15 cores (about 13 busy), AGAMA OpenMP.
- **Full run:** all 50k mock stars for 1 Gyr in the true potential (`make_potential` with
  `rotation=`, inertial frame, accuracy 10⁻⁸, snapshots every 10 Myr).
  - Wall time 482 s (8.0 min), 104 orbits/s, peak memory 0.29 GB.
  - Jacobi integral |ΔE_J/E_J|: median 2×10⁻⁷, 99th percentile 1.5×10⁻⁶. One orbit reached
    4×10⁻², probably a close pass by the central black hole.
- **Speed vs accuracy** (5000 random stars, 1 Gyr, projected to 50k; the map change is the rms
  over 1 kpc face-on pixels of ⟨v_R⟩ after 1 Gyr, relative to the exact run):

  | set-up | 50k × 1 Gyr | force error / bar force (rms) | by R: 0–1, 1–3, 3–6, 6–10 kpc | map change ⟨v_R⟩ |
  |---|---|---|---|---|
  | exact, acc 10⁻⁸ | 8.2 min | — | — | — |
  | exact, acc 10⁻⁷ / 10⁻⁶ (chaos floor) | 5.2 / 3.2 min | — | — | 1.5 / 3.1 km/s |
  | merged, mmax 8, no symmetry, 10⁻⁶ | 1.6 min | 0.13% | 0.2, 0.1, 0.1, 0.2% | 4.4 km/s |
  | merged, mmax 8, "reflection", 10⁻⁶ | 1.6 min | 0.13% | same | 4.4 km/s |
  | merged, **mmax 4**, no symmetry, 10⁻⁶ | 0.8 min | **7%** | 3.5, 6, 12, 9% | 16.5 km/s |
  | merged, mmax 8 / 6, bisymmetric, 10⁻⁶ | 0.8 / 0.6 min | 15% | 9, 11, 26, 27% | 21 / 17 km/s |
  | merged, mmax 4, bisymmetric, 10⁻⁶ / 10⁻⁵ | 0.45 / 0.34 min | 17% | 10, 12, 28, 28% | 20 / 22 km/s |
  | merged, mmax 2, bisymmetric, 10⁻⁵ | 0.21 min | 26% | 14, 25, 38, 33% | 25 km/s |

  - Merging into one CylSpline costs 0.1 s per trial potential and loses nothing at mmax 8.
  - Au18's bar has real m = 6–8 content near its ends: mmax 4 loses 12% of the bar force at
    3–6 kpc.
  - "Bisymmetric" also drops the odd-m and z-asymmetric parts, which at 3–10 kpc are a quarter
    of the bar force.
  - Every set-up faster than 1.6 min changes the evolved kinematic maps by 15–25 km/s rms,
    5–8× the chaos floor. That is comparable to the bar's own v_R quadrupole (about 35 km/s).
  - **Conclusion:** at mmax 8 the merged CylSpline at accuracy 10⁻⁶ gives a 5× speed-up and is
    accurate. Lower mmax buys another factor of 2–4 but makes the potential a different galaxy.
    MCMC needs reductions elsewhere: fewer stars, a shorter integration, an emulator
    (see the chat discussion).

## 2026-10-27 — Test 1: steady-state likelihood, first design and diagnosis

- Scripts:
  - `scripts/steady_state_likelihood.py` (per-star KDE likelihood);
  - `scripts/test1_qmc.py` (truth plus 64 Sobol points in |Ω_b|, φ_b, ln S, ln bar_frac);
  - `scripts/test1_gp_emulator.py`;
  - diagnostics: `scripts/diagnose_pattern_speed_maps.py`, `scripts/diagnose_velocity_distributions.py`,
    `scripts/diagnose_moment_likelihood.py`.
- Figures: [⟨v_R⟩/⟨v_φ⟩ maps, observed vs library](figures/test1_omega_maps.png);
  velocity distributions [along the bar axis](figures/test1_vdist_major_axis.png) and
  [at 45°](figures/test1_vdist_45deg.png).
- Likelihood v1 (per-star KDE):
  - Stars are rotated by the trial bar angle so that the assumed bar lies along x; the
    potential's bar is fixed along x.
  - Orbits run for 1 Gyr in the bar frame; snapshots from 0.25 to 1 Gyr, 4-fold folded, give
    124 samples per orbit (6.2M library points), weighted by 1/τ for the R < 8 kpc cut.
  - log p(v_i|x_i) uses the 1000 nearest library points in position (adaptive), leaving out
    star i's own orbit.
  - Cost: about 1 minute per evaluation (50 s integration, 8 s likelihood).
- Pilot (other parameters at the truth): the integrator-noise floor of logL is about ±12.

  | trial | ΔlogL vs truth |
  |---|---|
  | Ω_b = 45 | −30 |
  | Ω_b = 30 | −1094 |
  | φ_b = 40° | −840 |
  | bar_frac = 1.5 | −1085 |
  | S = 1.3 | **+62** |

- Why the KDE likelihood fails:
  - The maps show Ω_b = 45 clearly weakens the v_R quadrupole (24 vs observed 34 km/s) and
    speeds up v_φ; Ω_b = 30 slows v_φ by 10–20 km/s.
  - The per-star KDE nevertheless barely penalises Ω_b = 45, because regional mismatches
    cancel: 0–2 kpc +355, 2–5 kpc −440, 5–6 kpc +223.
  - In the inner kpc every library is 5–8% hotter than the data in all three components,
    because the unsoftened AGAMA centre is deeper than Arepo's. That dispersion mismatch is worth
    about as much as a 10 km/s error in the mean.
  - Per-star Δlog p has an rms of 1.2, and the 1% most extreme (tail) stars swing the total by
    hundreds.
- Moment χ² on the same libraries (bar-frame cells in R, φ mod π and |z|, data folded like the
  library, errors σ/√N, N ≥ 20): ΔlogL vs truth.

  | likelihood | Ω_b = 45 | Ω_b = 30 |
  |---|---|---|
  | means only | **−195** | **−257** |
  | means + dispersions, all R | −45 | −801 |
  | means + dispersions, R ≥ 1.5 kpc | −138 | −494 |

  χ²/dof at the truth: means 1.24, dispersions about 3.6. The dispersion value is inflated
  because σ/√(2N) underestimates the error for non-Gaussian distributions; bootstrap errors
  are needed.
- **KDE log-likelihood surface from the 55 finished Sobol points**
  ([corner plot](figures/test1_logl_corner.png), `scripts/plot_test1_logl_corner.py`;
  [emulator 1D slices](figures/test1_emulator_slices.png)):
  - φ_b is well recovered: the profile peaks at 26–27° (truth 27°).
  - Ω_b is broadly 33–44 km/s/kpc, with bumps at about 37 and 41.
  - S ≥ 1 is preferred, and bar_frac covers 0.6–1.1 with a second mode near 0.65.
  - Two *evaluated* points beat the truth by much more than the ±12 noise:

    | Ω_b | φ_b | S | bar_frac | logL − truth |
    |---|---|---|---|---|
    | 36.9 | 27.7° | 1.33 | 0.81 | +158 |
    | 43.1 | 29.9° | 1.17 | 1.04 | +123 |

    So this likelihood's maximum is not at the truth: it prefers a longer bar, which is
    weaker or faster. That fits the missing DM bar response, which grows outward.
  - The GP with 55 points is too coarse for detail: leave-one-out rms 330 near the maximum and a
    fitted noise term of about 185. The double peak in Ω and the emulator maximum
    (Ω = 41, +373 above the truth) cannot be trusted.
- **Decision (provisional, superseded by the follow-up below):** the primary likelihood becomes a moment-based χ² on cell means (dispersions only
  with bootstrap errors, and R ≥ 1.5 kpc while the centre is unsoftened). The per-star KDE is kept
  only for comparison. The KDE emulator run was stopped at 56/65 points; its totals are in
  `results/test1_qmc/evals.csv`.

---

## Next steps

1. ~~Pattern speed~~, ~~ground-truth potential~~ and ~~trial-potential family~~: done
   (2026-10-26). Still to decide: spherical vs non-spherical DM (caveat 1 above).
2. Test 1 itself: integrate the perfect mock in the true potential (the null test) and in
   trial potentials with m ≠ 0 rescaled in amplitude, length and Ω_p. Choose the
   discrepancy statistic and the integration time.
3. Add realism step by step: distance errors (forward-modelled), RV jitter, selection
   function, period–age tracer weighting.
