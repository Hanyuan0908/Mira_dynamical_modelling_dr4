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
| `mock/Nbody/t_t0_7` | provided by you (NEMO format) | isolated barred N-body galaxy, t = 7 (6.84 Gyr): disc, bulge, DM; includes the simulation's potential and accelerations |
| `mock/nbody_{stars,dm}_mwframe.hdf5` | `scripts/prepare_nbody_mwframe.py` | N-body particles in the MW-like frame (comp 0 disc, 1 bulge, 2 DM), with pot_sim and acc_sim |
| `mock/nbody_potential/{dm,dm_spherical,bulge,disc}.ini` | `scripts/build_nbody_potential.py` | N-body ground-truth potential |
| `mock/nbody_mock_perfect_R8_N50k.hdf5` | `scripts/make_mock_perfect.py --galaxy nbody` | N-body Test 1 mock: 50k stars (40,023 disc, 9,977 bulge), R < 8 kpc |
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

## 2026-10-27 — Switch to an isolated N-body galaxy for first-principle tests

- Why: Au18's ground-truth potential is not the one its stars are in equilibrium with. The
  problems are softening of 0.37 kpc (the inner ~2 kpc is wrong) and missing non-axisymmetric
  DM. Every likelihood is biased by this, which is a generic limitation of hydro simulations.
  The N-body model has better inner resolution and no gas.
  Auriga is kept for later tests of age-dependent selection (hotter vs colder populations).
- Code: `scripts/galaxies.py` is a registry, and the analysis scripts take `--galaxy au18|nbody`:
  `plot_au18_density.py`, `plot_au18_velocity_maps.py`, `plot_au18_vR_distance_error.py`,
  `measure_au18_pattern_speed.py`, `check_au18_potential.py`, `make_mock_perfect.py`.
  "au18" in those file names is historical. New: `scripts/nbody.py` (NEMO reader and units),
  `scripts/prepare_nbody_mwframe.py`, `scripts/build_nbody_potential.py`.
- Snapshot: `mock/Nbody/t_t0_7`.
  - Units are kpc, km/s and G = 1, so the mass unit is 232,508 M☉ and t = 7.0 is 6.84 Gyr.
    The units were checked with the enclosed-mass v_c (218–230 km/s) and the stored
    accelerations.
  - Components are identified by particle mass:

    | component | particles | mass |
    |---|---|---|
    | disc | 1.7M | 5.11×10¹⁰ M☉ |
    | bulge | 0.3M | 1.30×10¹⁰ M☉ |
    | DM | 1.0M | 8.72×10¹¹ M☉ |

  - The file stores the simulation's own Potential and Acceleration per particle (float32).
- MW-like frame (same recipe as Au18):
  - The stellar centre is 1.4 kpc from the mean particle position.
  - The disc spins anticlockwise in the file and was flipped to the MW sense.
  - Bar: A2_max = 0.40 at R = 3.4 kpc; R_bar = 5.00 kpc (A2 half max); phase rotated from 41.6°
    to −27.0°.
  - v_φ = 212 km/s at R = 8 kpc.
  - Figures: [density](figures/nbody_stellar_density_mwframe.png),
    [m = 2 profile](figures/nbody_bar_fourier.png),
    [face-on kinematics](figures/nbody_velocity_maps_faceon.png).
  - The v_R quadrupole is about ±25 km/s in the map and 15.8 km/s in the m = 2 fit, weaker than
    Au18's 34.9 km/s. The outer disc is smooth, with almost no spiral or lopsided pattern.
- Distance-error test ([figure](figures/nbody_vR_distance_error.png)): σ_d/d = 0.2 rotates the
  v_R quadrupole from −25.5° to −6.2° and raises its amplitude from 15.8 to 21.1 km/s, the same
  behaviour as in Au18.
- Potential (`scripts/build_nbody_potential.py`; spheroids as Multipole and the disc as
  CylSpline, as in AGAMA's `py/example_actions_nbody.py`):

  | component | expansion | parameters |
  |---|---|---|
  | DM | Multipole, no symmetry | lmax 8, gridSizeR 30, rmin 0.5, rmax 600 (`dm_spherical.ini`: lmax 0) |
  | bulge | Multipole, no symmetry | lmax 8, gridSizeR 30, rmin 0.03, rmax 20 |
  | disc | CylSpline, no symmetry | mmax 8, gridSizeR 30, gridSizez 30, Rmin 0.1, Rmax 50, zmin 0.05, zmax 50 |

  - rmin ≈ the 0.1% radius and rmax just beyond the outermost particle.
  - lmax 8 is where the AGAMA reference documentation says N-body noise saturates.
  - The disc grid is the Au18/Nexus grid.
  - v_c(8 kpc) = 222 km/s.
- Benchmarks (`scripts/check_au18_potential.py --galaxy nbody`, `results/nbody_potential_checks.json`):
  - Stellar density ([figure](figures/nbody_potential_density_check.png)): face-on median ratio
    1.001 with 0.03 dex Poisson scatter; Σ(R) within −3/+5%; A2 0.402 vs 0.406, phase identical.
  - Rotation curve ([figure](figures/nbody_rotation_curve_check.png)): the Jeans-predicted
    stellar ⟨v_φ⟩ matches the measured one at −0.8% median and 0.9% rms for R > R_bar.
    Asymmetric drift is 10 km/s at 8 kpc.
  - Potential vs the simulation's ([figure](figures/nbody_potential_vs_sim.png)): 400 (km/s)² rms
    over 1–20 kpc, 0.25% of the potential range (Au18: 1700 (km/s)², plus offsets of up to
    1.5×10⁴ in the centre).
  - Forces vs the simulation's accelerations ([figure](figures/nbody_force_check.png)):
    - median error 2% at r = 1–20 kpc and 3% at 0.5–1 kpc; it rises sharply inside 0.5 kpc,
      where the simulation is softened (estimated 0.05–0.1 kpc);
    - bar torque F_φ: simulation = 0.9995 × AGAMA, with binned means equal to ±1%;
    - the AGAMA-on-simulation slope of 0.959 is regression dilution from N-body force noise
      (4% of the variance), not a real deficit;
    - torque shares: disc 78.5%, bulge 14%, DM 3%. Multipole `smoothing` (1 → 0) changes nothing.
- **How non-spherical is the N-body DM halo?** (ad hoc analysis; `dm_axisymmetric.ini` added as
  an alternative component)
  - Non-axisymmetry (the halo's bar response):
    - the DM m = 2 potential is 3–5% of the stellar m = 2 for R = 0.5–6 kpc, roughly aligned with
      the bar but twisted by up to about 26°;
    - it carries 3% of the bar torque;
    - beyond 8 kpc its amplitude is only about 50 (km/s)² with an unrelated phase, i.e. noise;
    - the DM particle density A2 is 0.02 inside 4 kpc (Poisson level about 0.01).
  - Shape: b/a = 0.98–1.00. Slightly oblate: c/a = 0.94 within 5 kpc, 0.96 within 10 kpc,
    0.99 within 100 kpc.
  - Against the simulation's forces:

    | DM symmetry | F_z error, R 2–4 kpc | F_z error, R 6–10 kpc (\|z\| 0.5–1.5) | F_z error, R 6–10 kpc (\|z\| 1.5–3) | bar torque k (sim = k × AGAMA) |
    |---|---|---|---|---|
    | spherical | −0.9% | **−2.5%** | **−2.5%** | 1.029 |
    | axisymmetric | +0.6% | 0.0% | +0.5% | 1.029 |
    | none (current) | +0.8% | 0.0% | +0.5% | 0.9995 |

    DM supplies 17–28% of F_z in the outer disc.
  - Conclusion: a spherical halo underestimates the vertical force in the disc by 2.5%, enough to
    put stars slightly out of vertical equilibrium. An axisymmetric halo is accurate vertically
    but drops the 3% DM share of the bar torque.
- **Decision (2026-10-27, user):** the N-body ground truth uses the **spherical DM halo**
  (`dm_spherical.ini`), because the MW set-up will use a spherical halo. The registry was updated
  (`components = dm_spherical, bulge, disc`) and the checks and resonances rerun.
  - v_c(8 kpc) = 221.6 km/s; the Jeans check and density checks are unchanged.
  - Potential residual vs the simulation: 504 (km/s)² rms.
  - Bar torque: simulation = 1.029 × AGAMA.
  - R_CR = 8.61 kpc.
  - Known, deliberate mismatches with what the stars feel: the vertical force in the outer disc
    is 2.5% low, and the halo's 3% share of the bar torque is missing. Expect bar_frac to come
    out about 3% above 1 and possibly small vertical non-stationarity, a realistic level of
    model mismatch for the MW.
- Pattern speed (Dehnen et al. 2023 method, `results/nbody_ground_truth.json`,
  [figure](figures/nbody_pattern_speed.png)):
  - **Ω_p = 26.02 ± 0.19 km/s/kpc**, using the automatic bar region 0.57–6.28 kpc. Fixed ranges
    give 25.4–25.5.
  - R_CR = 8.63 kpc and OLR = 14.8 kpc. R_CR / R_bar = 1.73 using the A2 half-maximum length, or
    1.37 using the method's bar region: a fairly slow bar.
  - A second snapshot (t_t1) would allow a finite-difference cross-check.
- Basic mock ([overview](figures/nbody_mock_perfect_R8_N50k_overview.png)):
  - 50k random star particles with R < 8 kpc (40,023 disc and 9,977 bulge); no errors and no
    selection; 15.5% have |z| > 1 kpc.
  - The v_R quadrupole is 16.9 ± 1.5 km/s with its zero line at −26.1 ± 3.1°.
  - The weaker bar signal means larger mocks (100–200k) may be needed for the likelihood tests.

## 2026-10-27 — Test 1 on the N-body galaxy: pilot, likelihood choice and global search

- Trial potentials (`scripts/potential_family.py --galaxy nbody`, checks in
  `results/nbody_potential_family_checks.json` and [figure](figures/nbody_potential_family_check.png)):
  - The spherical DM halo is fixed. The bulge and disc are each axisymmetrised with their own
    expansion (Multipole for the bulge, CylSpline for the disc), and their non-axisymmetric parts
    are rescaled together with bar_frac and S and rotated as before.
  - The same six checks as Au18 pass: truth reproduced to 4×10⁻¹⁶, v_c changes ≤ 1.3×10⁻⁴,
    exact scaling and rotation, inertial and rotating-frame orbits agree.
- Speed:
  - The exact N-body potential costs 7.2 min per evaluation of 50k × 1 Gyr. The cuspy bulge needs
    many small steps; orbits inside 2 kpc take 80% of the time.
  - Merging the static barred bulge and disc into one CylSpline with a finer inner grid
    (`NBODY_GRID`: 36×36, Rmin 0.03, zmin 0.02) gives a bar-force error of 0.2% rms and runs
    1.6× faster.
  - With integrator accuracy 10⁻⁵, an evaluation takes about 3.5 min.
- `scripts/steady_state_likelihood.py` now also computes, from the same 1000 neighbours per star:
  - a per-star local Gaussian N(v_i | μ_i, C_i);
  - the moment χ² in polar cells;
  - and saves μ_i, C_i, τ and 0.5 kpc box sums (count, mass, Σv, Σvv for v_R, v_φ, v_z).
- **Pilot** (`scripts/test1_nbody_pilot.py`, `scripts/analyse_test1_nbody_pilot.py`;
  [figure](figures/test1_nbody_pilot.png), `results/test1_nbody/pilot_summary.json`).
  ΔlogL = offset − truth; a negative value means the truth is preferred.

  | offset | KDE | **local Gaussian** | Student-t ν=5 | Gaussian + 2% outliers | moments: means | moments: means+disp (R≥1) |
  |---|---|---|---|---|---|---|
  | noise (accuracy 10⁻⁵ vs 3×10⁻⁶) | +26.5 | −5.7 | +1.2 | −4.9 | +0.7 | −3.0 |
  | Ω_b = 20 | −3394 | −1813 | −1032 | −1865 | −907 | −2931 |
  | Ω_b = 32 | −350 | −596 | −1189 | −602 | −378 | −509 |
  | φ_b = 17° | −375 | −265 | −277 | −271 | −110 | −238 |
  | φ_b = 37° | −448 | −281 | −303 | −277 | −94 | −203 |
  | bar_frac = 0.7 | **+200** | −68 | −227 | −69 | −39 | −50 |
  | bar_frac = 1.4 | −468 | −100 | **+137** | −102 | −41 | −159 |
  | S = 0.8 | **+4** | −68 | −164 | −75 | −51 | −80 |
  | S = 1.25 | −481 | −163 | **+6** | −164 | −86 | −205 |

  - The local Gaussian, the outlier mixture and all moment versions prefer the truth on every
    axis, by ≥ 12× the noise.
  - The local-Gaussian penalties come mostly from the bar region (R = 2–6 kpc).
  - The KDE again prefers a weaker bar, gaining at every radius, and has 5× the noise. The
    Student-t prefers a stronger and longer bar.
  - With the N-body potential the moment dispersions behave, which confirms that the Au18
    failures came from its potential mismatch.
- **Decision:** the per-star local-Gaussian p(v|x) is the primary likelihood. It is a p(v|x) form,
  suitable for selection functions and per-star error convolution later.
- **Global search** (`scripts/test1_nbody_qmc.py`):
  - Points: the truth plus 64 Sobol points over |Ω_b| 18–38 km/s/kpc, φ_b 14–44°,
    S 0.75–1.45 and bar_frac 0.55–1.6 (both log-uniform).
  - Launched 2026-10-27 14:06, expected to finish around 18:00.
  - Saved per point: totals in `results/test1_nbody/qmc_evals.csv`; per-star log p, μ, C, τ,
    polar cell moments and box sums in `~/mira_dr4_runs/test1_nbody/NNN.npz`; full trajectories
    in `traj_NNN.npy` (bar frame, float32), all outside Dropbox.

---

## 2026-10-27 — Test 1b: global search results (N-body, perfect mock)

- The global search finished in 1 h 57 min (109 s per point). Totals and bootstrap moment columns
  are in `results/test1_nbody/qmc_evals.csv`; arrays and trajectories in `~/mira_dr4_runs/test1_nbody/`.
- Moment-map likelihood (`scripts/test1_nbody_moments.py`):
  - Folded 3D bar-frame cells (R 0.5 kpc × φ 15° × |z| 0/0.5/1/2/∞), N_obs ≥ 20, data errors from
    200 bootstrap resamplings.
  - χ²/dof at the truth: 1.11 for means (all R) and 1.29 for means + dispersions (R ≥ 1 kpc).
- **Raw evaluations:** the truth is the **best of all 65 points** for the local Gaussian and for
  every moment-map variant. The best other point is 7 lower (Gaussian), 37 lower (means) and
  64 lower (means + dispersions). For the KDE the truth ranks 2nd, behind a point +2 better
  (noise about 26).
- **Corner plots** ([means + dispersions, R ≥ 1](figures/test1_nbody_logl_corner_chi2b_meansdisp_R1.png);
  also `_chi2b_means`, `_logL_gauss`, and `test1_nbody_logl_corner.png` for the KDE):
  - The Ω_b profile peaks at the truth (26), and φ_b at 27–30°.
  - **S and bar_frac are anti-correlated along a ridge:** a shorter, stronger bar resembles a longer,
    weaker one. The ridge runs from (S 0.75, f 1.55) through the truth to (S 1.3, f 0.65). The best
    non-truth evaluations lie on it, e.g. S 0.82, f 1.50.
- **The 64-point emulator is not adequate:**
  - With about 2.8 points per dimension, it cannot resolve a peak 1–2 km/s/kpc wide in Ω.
  - Its leave-one-out error near the maximum is 85–240 in logL, against noise of 1–6.
  - It under-predicts the truth by 44–108 and places its maximum at the box corner along the ridge
    (S = 0.75, f ≈ 1.56), where the evaluated points are actually 64 below the truth.
  - The emulator "posteriors" (`results/test1_nbody/emulator_results*.json`) are therefore
    artefacts. They must not be read as parameter recovery.
- **Next:** a second, focused design (narrow Ω_b and φ_b, keep the full S–bar_frac ridge, or
  reparametrise along and across the ridge), and fit the GP only to points near the maximum.

## 2026-10-27 — Test 1b: tempered χ² and what bar length and strength change

- **Tempered χ²** (Schwarzschild-style):
  - log L = −χ²_meansdisp,R≥1 / (2 √(2 N_data)), with N_data = 6 moments × 451 cells, so
    √(2 N_data) = 73.6. Column `logLs_meansdisp_R1` in `results/test1_nbody/qmc_evals.csv`.
  - [Corner plot](figures/test1_nbody_logl_corner_logLs_meansdisp_R1.png); emulator posterior in
    `results/test1_nbody/emulator_results_logLs_meansdisp_R1.json`.
  - The emulator is now adequate: leave-one-out error 6.7 in tempered units.

  | parameter | posterior | truth | offset |
  |---|---|---|---|
  | Ω_b | 26.19 +1.33/−1.13 km/s/kpc | 26.02 | −0.14σ |
  | φ_b | 27.7 ± 4.5° | 27 | −0.16σ |
  | S | 0.96 +0.22/−0.16 | 1 | +0.21σ |
  | bar_frac | 0.92 +0.40/−0.23 | 1 | +0.25σ |

  - The S–bar_frac anti-correlation ridge is clearly visible.
  - Caveats:
    - √(2N) tempering is a conservative convention, widening errors by 8.6× over the formal
      ones. Rescaling by χ²_min/dof (1.2–1.3) would widen them by only 1.1×.
    - The right tempering should be calibrated from the scatter of best fits over independent
      mock realisations.
    - The emulator still puts its maximum at the box corner (S 0.75, f 1.57, +1.4 tempered
      units above the emulated truth). This is within its error; the evaluated ridge point
      (S 0.82, f 1.50) is −0.7 below the truth.
- **What S and bar_frac change** (`scripts/diagnose_bar_length_strength.py`; Ω_b and φ_b at the
  truth; [maps](figures/test1_nbody_bar_length_strength_polar_pull_maps.png), [radial profiles](figures/test1_nbody_bar_length_strength_radial_profiles.png),
  `results/test1_nbody/diag_SF_chi2.json`). Untempered ΔlogL vs the truth:

  | model | total | ⟨v_R⟩ | ⟨v_φ⟩ | ⟨v_z⟩ | σ_R | σ_φ | σ_z |
  |---|---|---|---|---|---|---|---|
  | ridge S 0.82, f 1.50 | −52 | −13 | −5 | −3 | −9 | −9 | −13 |
  | ridge S 1.25, f 0.70 | −95 | −7 | −37 | 0 | −13 | −35 | −3 |
  | f 0.7 | −61 | −10 | −23 | 0 | −10 | −5 | −12 |
  | f 1.4 | −192 | −8 | −31 | −4 | −59 | −64 | −26 |
  | S 0.8 | −94 | −22 | −25 | −5 | −16 | −14 | −11 |
  | S 1.25 | −236 | −14 | −70 | +1 | −55 | −70 | −28 |

  - The truth is best. The two ridge points cost much less than single offsets of similar size,
    which is the degeneracy.
  - **Correction (2026-10-27):** an earlier version of this entry said the moments carry no
    information on S or bar_frac. That is wrong: only the azimuthal means are identical. The
    m = 2 parts of all four maps, dispersions included, respond to S and bar_frac (see the next
    entry and the table above: dispersions give −123 of −192 for f 1.4 and −125 of −236 for
    S 1.25).
  - Along the ridge, stretching and weakening the bar keeps the m = 2 amplitude at 2–4 kpc similar.
    S 0.82/f 1.5 under-predicts the v_R quadrupole at 3–4 kpc (18 vs 23 km/s); S 1.25/f 0.70 fails
    mainly at 7–8 kpc (ΔlogL ≈ −20 there).

## 2026-10-27 — Face-on mean and dispersion maps for different bar lengths and strengths

- Script: `scripts/plot_bar_length_strength_polar_maps.py`. No new integrations: it reuses the cached library cell
  moments. Numbers are in `results/test1_nbody/SF_maps_numbers.json`.
- Figures (midplane |z| < 0.5 kpc, bar frame):
  [values](figures/test1_nbody_bar_length_strength_polar_maps_values.png),
  [non-axisymmetric part](figures/test1_nbody_bar_length_strength_polar_maps_nonaxisymmetric.png),
  [model − truth model](figures/test1_nbody_bar_length_strength_polar_maps_model_minus_truth.png),
  [m = 2 amplitude and phase profiles](figures/test1_nbody_bar_length_strength_m2_profiles.png).
- **All four maps carry strong m = 2 structure.** Truth model, m = 2 amplitude at
  R = 2.25 / 3.25 / 4.25 / 5.25 kpc:

  | quantity | m = 2 amplitude [km/s] | phase |
  |---|---|---|
  | ⟨v_R⟩ | 17.2 / 26.5 / 16.8 / 6.7 | ≈ 45° |
  | ⟨v_φ⟩ | 15.8 / 41.1 / 30.3 / 8.9 | ≈ ±90° |
  | σ_R | 20.7 / 23.1 / 9.8 / 3.3 | 0°, peaking along the bar |
  | σ_φ | 13.5 / 6.2 / 20.5 / 9.5 | changes sign near R ≈ 3.2 kpc |

- **The m = 2 amplitudes respond to S and bar_frac.**
  - At R = 4.25 kpc the σ_R m = 2 amplitude is 6.1 (f 0.7), 9.8 (truth), 14.2 (f 1.4),
    5.1 (S 0.8) and 12.8 (S 1.25).
  - In ⟨v_R⟩ at 4.25 kpc: 12.1, 16.8, 21.6, 8.6, 23.3.
  - The model − truth maps show coherent m = 2 residuals for the single offsets: a weaker or
    stronger quadrupole in ⟨v_R⟩, ±5–8 km/s in σ_R at the bar ends, and an outward shift of the
    patterns for S 1.25.
- **The degeneracy is in all four maps at once.** rms(model − truth model) over cells with
  |z| < 1, R ≥ 1 and N ≥ 20, for ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ:

  | model | rms [km/s] |
  |---|---|
  | ridge S 0.82, f 1.50 | 3.2, 2.9, 1.8, 2.1 |
  | ridge S 1.25, f 0.70 | 2.9, 3.0, 1.9, 2.0 |
  | single offsets | 3.3–4.6, 3.4–3.9, 2.3–2.7, 2.1–2.5 |

  - Library noise floor for a difference of two libraries: √2 × the orbit-bootstrap error, i.e.
    2.1–3.0 km/s for means and 1.6–1.8 km/s for dispersions.
  - Typical data error per cell: 9.4, 6.9, 6.7, 4.9 km/s.
  - So the differences are 0.3–0.5σ per cell and are detected only by summing about 450 cells.
- **Means only vs means + dispersions** (emulator posteriors, 64 Sobol points; offsets in σ):

  | likelihood | Ω_b | φ_b | S | bar_frac | LOO near max |
  |---|---|---|---|---|---|
  | means, untempered | 25.03 (+6.7σ) | 25.8 (+1.5σ) | 1.025 (−0.8σ) | 0.677 (+11.9σ) | 85 |
  | means + disp R ≥ 1, untempered | 26.29 (−2.0σ) | 29.4 (−5.1σ) | 0.751 (+301σ) | 1.556 (−18σ) | 241 |
  | means, tempered (√2N = 55.6) | 25.7 +1.4/−1.3 | 26.9 ± 5.8 | 1.00 +0.22/−0.16 | 0.82 +0.42/−0.20 | 2.5 |
  | means + disp R ≥ 1, tempered (73.6) | 26.2 +1.3/−1.2 | 27.2 ± 4.3 | 0.95 +0.23/−0.15 | 0.92 +0.40/−0.24 | 6.7 |

  - Both untempered emulators are biased, because their emulator error ≫ the logL scale.
  - Both tempered versions recover the truth within 0.6σ.
  - Adding the dispersions tightens φ_b. In the pilot it increased the penalties for wrong S and
    bar_frac: S 1.25 from −86 to −205 and f 1.4 from −41 to −159.
  - The posterior medians vary by about 0.1σ between emulator and emcee runs (emcee is not seeded).

## 2026-10-27 — Tempered χ² is now the standard moment likelihood

- From now on every moment-map log-likelihood is reported tempered:
  log L = −χ² / (2 √(2 N_data)), with N_data = (number of moments) × (number of cells) at the truth.
  On this scale ΔlogL = −0.5, −2 and −4.5 mark 1σ, 2σ and 3σ.
- Applied to:
  - all four χ² variants: columns `logLs_means` (√2N = 55.6), `logLs_means_R1` (52.0),
    `logLs_meansdisp` (78.6) and `logLs_meansdisp_R1` (73.6) in `results/test1_nbody/qmc_evals.csv`;
  - the S–bar_frac diagnostic (`results/test1_nbody/diag_SF_chi2.json`, profile figure);
  - the pilot table (`results/test1_nbody/pilot_summary.json`).
- Global search: the truth is still the best of 65 evaluated points for every variant. The second
  best is −0.54 to −0.88 tempered, i.e. about 1σ away.
- S–bar_frac diagnostic, tempered ΔlogL (means + dispersions, R ≥ 1 kpc; Ω_b and φ_b at the truth):

  | model | total | ⟨v_R⟩ | ⟨v_φ⟩ | ⟨v_z⟩ | σ_R | σ_φ | σ_z |
  |---|---|---|---|---|---|---|---|
  | ridge S 0.82, f 1.50 | −0.71 | −0.18 | −0.06 | −0.04 | −0.12 | −0.13 | −0.18 |
  | ridge S 1.25, f 0.70 | −1.29 | −0.10 | −0.50 | 0.00 | −0.18 | −0.48 | −0.04 |
  | f 0.7 | −0.82 | −0.13 | −0.31 | 0.00 | −0.14 | −0.07 | −0.17 |
  | f 1.4 | −2.60 | −0.11 | −0.41 | −0.06 | −0.80 | −0.87 | −0.36 |
  | S 0.8 | −1.27 | −0.30 | −0.34 | −0.07 | −0.22 | −0.19 | −0.15 |
  | S 1.25 | −3.20 | −0.19 | −0.95 | +0.02 | −0.75 | −0.96 | −0.38 |

- Pilot, tempered moment ΔlogL (means all R / means R ≥ 1 / means + dispersions R ≥ 1):

  | offset | means all R | means R ≥ 1 | means + disp R ≥ 1 |
  |---|---|---|---|
  | Ω_b = 20 | −16.3 | −17.4 | −39.8 |
  | Ω_b = 32 | −6.8 | −7.3 | −6.9 |
  | φ_b = 17° | −2.0 | −2.0 | −3.2 |
  | φ_b = 37° | −1.7 | −1.6 | −2.8 |
  | bar_frac = 0.7 | −0.7 | −0.6 | −0.7 |
  | bar_frac = 1.4 | −0.7 | −0.9 | −2.2 |
  | S = 0.8 | −0.9 | −1.0 | −1.1 |
  | S = 1.25 | −1.6 | −1.7 | −2.8 |

  The noise floor is about 0.0 on this scale.

## 2026-10-27 — Tempered corner plots and square-pixel maps for bar length and strength

**Tempered corner plots (all R, 64 Sobol points + truth, GP emulator).**
- Figures: [means](figures/test1_nbody_emulator_corner_logLs_means.png) and
  [means + dispersions](figures/test1_nbody_emulator_corner_logLs_meansdisp.png). Profiles of the
  evaluated points: [means](figures/test1_nbody_logl_corner_logLs_means.png),
  [means + dispersions](figures/test1_nbody_logl_corner_logLs_meansdisp.png).
- Numbers: `results/test1_nbody/emulator_results_logLs_{means,meansdisp}.json`.

  | likelihood | √(2N_data) | Ω_b | φ_b [°] | S | bar_frac | LOO rms |
  |---|---|---|---|---|---|---|
  | means | 55.6 | 25.68 +1.53/−1.30 | 26.7 +6.2/−5.7 | 1.008 +0.216/−0.157 | 0.804 +0.428/−0.183 | 2.5 |
  | means + disp | 78.6 | 26.26 +1.32/−1.13 | 27.7 +4.6/−4.6 | 0.941 +0.229/−0.141 | 0.936 +0.403/−0.251 | 6.4 |

- All four truth values lie within 0.65σ in both cases.
- S and bar_frac follow an anticorrelated ridge that runs into the prior box edges:
  - means: the bar_frac marginal peaks at the lower edge (0.55);
  - means + disp: the S marginal peaks at the lower edge (0.75).
- At the edges the emulator is extrapolating. For means + disp its maximum is at the box corner
  (S 0.75, bar_frac 1.56), logL −25.49. That is *below* the evaluated truth (−25.38). The best
  evaluated non-truth point, (S 0.82, f 1.50), is at −26.20.

**Square-pixel face-on maps.**
- Scripts: `scripts/run_bar_length_strength_models.py` integrates S 1.25, S 0.8, f 1.4 and f 0.7 at the true
  Ω_b and φ_b (about 2 min each). `scripts/plot_bar_length_strength_square_maps.py` makes the maps.
- Figures: [bar length](figures/test1_nbody_square_maps_bar_length_vs_data.png) and
  [bar strength](figures/test1_nbody_square_maps_bar_strength_vs_data.png). Pixels are 0.5 kpc, |z| < 1 kpc, in the
  bar frame; the data are folded by point symmetry.
- Rows:
  1. data;
  2. ground-truth model;
  3. the two alternatives;
  4. bootstrap data error;
  5. each model − data;
  6. each alternative − ground-truth model.
- Numbers: `results/test1_nbody/SF_square_maps_numbers.json`.
- Sample: 42,249 stars with |z| < 1; 788 pixels with ≥ 20 stars; median 37 stars per pixel.
- Median bootstrap error per pixel, for ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ: **9.6, 7.0, 5.8, 4.5 km/s**.
- Ground-truth model vs data: rms 11.0, 9.1, 6.5, 6.8 km/s; χ²/pixel 1.35, 1.23, 1.48, 1.80.
  These use data errors only, without library noise.
- **Alternative − ground-truth model:** rms 2.3–4.7 km/s. The pixel-scale library noise in these
  differences is 1.2–2.1 km/s (neighbour-difference estimate). Region averages, in km/s:

  | model | ⟨v_R⟩ quadrupole, R 1.5–3.5 / 3.5–5.5 | σ_R at bar ends, \|x\| ≈ 3.5 / 4.5 | ⟨v_φ⟩ at \|x\| ≈ 4.5 | ⟨v_φ⟩ on minor axis, \|y\| ≈ 2.5 |
  |---|---|---|---|---|
  | (truth's own quadrupole) | 11.3 / 8.0 | | | |
  | S 1.25 | +2.1 / +3.7 | +3.3 / +4.9 | −5.9 | +3.4 |
  | S 0.8 | −3.4 / −4.5 | −4.7 / −4.2 | +6.5 | −4.2 |
  | f 1.4 | +1.1 / +1.7 | +5.6 / +6.4 | −6.9 | +5.1 |
  | f 0.7 | −1.8 / −2.4 | −4.6 / −4.0 | +7.2 | −3.6 |

  "Quadrupole" means the mean of Δ⟨v_R⟩ × sign(sin 2φ), i.e. projected on the truth's own pattern.
  - A longer bar acts like a stronger one in σ_R and ⟨v_φ⟩: same sign, similar size. This is the
    S–bar_frac ridge.
  - ⟨v_R⟩ separates them. ±25% in S changes the ⟨v_R⟩ quadrupole about twice as much as ±40% in
    bar_frac (3.7–4.5 vs 1.7–2.4 km/s at R 3.5–5.5).
- **These differences are below the per-pixel data error.** Model − truth is 2.3–4.7 km/s against
  4.5–9.6 km/s errors. They are not visible by eye in the model − data rows, whose rms changes by only
  −0.04 to +0.37 km/s relative to the truth's.
- Δχ² vs the ground truth on these 2D maps, tempered by √(2·4·788) = 79.4:

  | model | ⟨v_R⟩ | ⟨v_φ⟩ | σ_R | σ_φ | total Δχ² | total ΔlogL |
  |---|---|---|---|---|---|---|
  | S 1.25 | +52 | +186 | +127 | +295 | +661 | −4.16 |
  | S 0.8 | +41 | +47 | +66 | +66 | +220 | −1.39 |
  | f 1.4 | +9 | +98 | +54 | +196 | +356 | −2.24 |
  | f 0.7 | +2 | +44 | +72 | +39 | +157 | −0.99 |

  - This ordering matches the 3D polar-cell diagnostic of the previous entry: S 1.25 −3.20,
    f 1.4 −2.60, S 0.8 −1.27, f 0.7 −0.82.
  - Shorter or weaker bars are penalised about 2–3× less than longer or stronger ones.
  - The mean-velocity maps alone (⟨v_R⟩ + ⟨v_φ⟩) give f 0.7 only Δχ² = +46, i.e. ΔlogL −0.29.

## 2026-10-27 — Does re-binning tighten the S / bar_frac constraints? (no new integrations)

**Set-up.**
- Script: `scripts/test1_binning_study.py`; numbers in `results/test1_nbody/binning_study.json`.
- Every saved library (65 global-search points and the 4 S/f variants) and the data are binned into
  unions of the saved 0.5 kpc boxes. Each box goes to a bin by its centre, with the same mapping for
  data and library.
- Folding: y ≥ 0, z ≥ 0. Boxes used: wholly inside R < 8 kpc, |z| < 4 kpc. Layers: |z| 0–0.5–1–2–4 kpc.
- Errors from 100 bootstraps; N_MIN 20.
- Statistics: the 6 velocity moments, plus **star counts** (library: 1/τ-weighted mass, normalised
  to the data total).
- Check: the box version of the current cells (polar15) gives χ²/N_data = 1.30 at the truth
  (exact polar likelihood: 1.29).

**Means + dispersions.** Untempered Δχ² vs the truth for S 1.25 / S 0.8 / f 1.4 / f 0.7, with the
tempered ΔlogL in brackets:

| binning | N_data | √(2N) | S 1.25 | S 0.8 | f 1.4 | f 0.7 |
|---|---|---|---|---|---|---|
| polar15 (≈ current) | 2526 | 71.1 | +407 (−2.86) | +138 (−0.97) | +339 (−2.38) | +95 (−0.67) |
| polar30 | 1608 | 56.7 | +340 (−3.00) | +114 (−1.01) | +343 (−3.02) | +116 (−1.03) |
| polar45 | 1152 | 48.0 | +323 (−3.37) | +136 (−1.41) | +320 (−3.33) | +124 (−1.29) |
| R 1 kpc × 15° | 1668 | 57.8 | +354 (−3.06) | +133 (−1.15) | +309 (−2.67) | +109 (−0.94) |
| 0.5 kpc squares | 2742 | 74.1 | +488 (−3.29) | +188 (−1.27) | +484 (−3.27) | +150 (−1.02) |
| 1 kpc squares | 1710 | 58.5 | +380 (−3.25) | +160 (−1.37) | +373 (−3.19) | +128 (−1.09) |
| adaptive (≥ 100 stars) | 1392 | 52.8 | +333 (−3.15) | +107 (−1.01) | +351 (−3.32) | +89 (−0.85) |
| Fourier m = 0, 2, 4 | 954 | 43.7 | +186 (−2.13) | +125 (−1.44) | +146 (−1.68) | +85 (−0.98) |
| Fourier m = 0, 2 | 762 | 39.0 | +166 (−2.13) | +93 (−1.19) | +131 (−1.68) | +56 (−0.72) |

- **Noise.** Over 60 data resamplings the realised Δχ² scatters with std 38–47 (polar15),
  30–34 (polar45) and 35–44 (1 kpc squares). In tempered units that is about ±0.3.
- **Untempered Δχ² (information).** It is largest for the finest bins (0.5 kpc squares). Coarser bins
  never raise it beyond the noise.
  - Fourier compression to m ≤ 4 loses about half for S 1.25 and f 1.4 (407 → 186, 339 → 146). The S/f
    signal has structure finer than m = 4, e.g. at the bar ends.
  - Squares beat polar15 (488 vs 407), because 15° sectors are about 1 kpc wide at R = 4.
- **Tempered ΔlogL.** Coarser bins gain ×1.2–1.9, mostly because the denominator shrinks.
  - polar45 vs polar15: √(2526/1152) = 1.48 at fixed Δχ², times the Δχ² retained (0.79–1.31).
  - The S 0.8 and f 0.7 gains are about 1–2× the noise.
- **Ridge point (Sobol #15: Ω 26.77, φ 27.6°, S 0.82, f 1.50).**
  - Means + dispersions Δχ² vs the truth is +3 (polar15), +59 (polar45), +209 (1 kpc squares).
  - The individual statistics swing by ±50–100 between binnings; σ_R ranges from −75 to +14.
  - So separating the ridge is noise-limited. Choosing a binning post hoc to break it would be tuning
    to noise.

**Star counts (density stationarity), a statistic we saved but never used.**
- [Figure](figures/test1_nbody_density_stationarity.png) (`scripts/plot_density_stationarity.py`): the
  time-averaged library density changes coherently in wrong potentials.
  - S 1.25: up to +20% beyond the bar ends (|x| ≈ 6–7 kpc) and −20% on the minor axis at R 5–7.
  - S 0.8: the opposite.
  - f 1.4 and f 0.7: ±10% along the bar and the minor axis at R 3–5.
  - rms change per 0.5 kpc pixel (R < 6, |z| < 1): 6.0–10.0%, vs 16% median Poisson error; the
    patterns span several kpc.
- Counts alone, Δχ² for S 1.25 / S 0.8 / f 1.4 / f 0.7:
  - polar15: +248 / +33 / +116 / +65 (N_data 421);
  - polar45: +188 / +64 / +75 / +78 (192);
  - 1 kpc squares: +253 / +60 / +98 / +99 (285).
  - χ²/N_data at the truth: 1.09–1.46.
- Means + dispersions + counts, tempered ΔlogL:

  | binning | S 1.25 | S 0.8 | f 1.4 | f 0.7 |
  |---|---|---|---|---|
  | polar15 | −4.27 | −1.12 | −2.96 | −1.05 |
  | polar45 | −4.93 | −1.92 | −3.81 | −1.94 |
  | 1 kpc squares | −5.01 | −1.74 | −3.73 | −1.80 |

- Effect on the 64 Sobol points: with counts, only 1 of them is within ΔlogL −2 of the truth for every
  non-Fourier binning (1–5 without counts), and 4–6 are within −4.5 (7–9 without).
- Counts add a positive amount for the ridge point in all 9 binnings (+23 to +76), but the ridge
  stays close: −0.18 (polar15), −1.08 (polar45), −1.84 (1 kpc squares).
- **Caveats.**
  - Counts depend directly on the selection function, so they need S(x) forward-modelled (Test 2a).
  - Library noise is not in the error budget.
  - Bins are box-centre approximations.
  - Fourier rows use diagonal bootstrap errors.

**Conclusion.**
- Re-binning adds no information: the finest well-populated bins keep the most Δχ².
- Under the √(2N_data) tempering, fewer bins give tighter tempered constraints at the same Δχ². The
  tempered posterior width therefore depends on the binning convention.
- ~~The real gain available from the saved data is the density term.~~
- **Retracted (2026-10-27, user):** counts are not usable. Gaia's selection function cannot be worked
  out, and if it could, the counts divided by the selection function would give the bar potential
  directly. The method exists precisely because velocities at fixed position do not need the spatial
  selection function. The density figure and numbers stay only as a record.

## 2026-10-27 — Fourier compression of the moment maps in φ (no new integrations)

**Method.**
- Script: `scripts/test1_fourier_compression.py`; numbers in
  `results/test1_nbody/fourier_compression.json`; [figure](figures/test1_nbody_fourier_compression_vs_mmax.png).
- Each ring (R 0.5 kpc), |z| layer and field is compressed separately. The fields are
  ⟨v_R⟩, ⟨v_φ⟩, ⟨v_z sign z⟩, σ_R, σ_φ, σ_z.
- Its 15° sector values (N ≥ 20) are fitted with even Fourier modes up to m_max by generalised least
  squares, weighted by the data's bootstrap errors.
- The same linear operator is applied to the data and the library, so χ²(m_max) is the part of the
  current sector χ² that lies in modes ≤ m_max.
- Check: with all modes (m_max = 12) it equals the current sector χ² to 10⁻¹². With all modes the
  compression is lossless.
- Inputs: the exact saved polar-cell moments (pilot offsets, ridge points, accuracy twin, all 65
  global-search points). Columns `chi2_fourier_*`, `logLs_fourier_*` in `qmc_evals.csv` (renamed
  2026-10-27 from the cryptic `chi2F*`/`logLs_F*`).

**Does truncation lose information?** Means + dispersions; Δχ² vs the truth on the data, tempered
ΔlogL in brackets, Ω_b and φ_b at the truth:

| m_max | N_data | S 1.25 | S 0.8 | f 1.4 | f 0.7 | ridge S 0.82/f 1.50 | ridge S 1.25/f 0.70 |
|---|---|---|---|---|---|---|---|
| 0 | 288 | +276 (−5.75) | −44 (+0.92) | +278 (−5.79) | −71 (+1.48) | +54 (−1.12) | +69 (−1.44) |
| 2 | 864 | +388 (−4.67) | +132 (−1.59) | +362 (−4.35) | +74 (−0.90) | +88 (−1.05) | +187 (−2.25) |
| 4 | 1434 | +447 (−4.17) | +205 (−1.91) | +383 (−3.57) | +136 (−1.27) | +130 (−1.21) | +213 (−1.99) |
| 12 (= current) | 3090 | +474 (−3.02) | +183 (−1.16) | +381 (−2.42) | +138 (−0.88) | +120 (−0.76) | +204 (−1.29) |

- **The realised Δχ² plateaus by m_max ≈ 4.** m_max = 4 keeps 94–113% of the full value for all six
  models, within the ±30–47 noise, with 46% of N_data.
- **Above m ≈ 4 the extra separation is library noise.** The noise-free separation (library − truth
  library) keeps rising at 0.10–0.11 per added coefficient, the same for every model.
  - Two libraries with independent sampling noise of 20–24% of the data error (measured earlier)
    would give 2 × 0.20²–0.24² = 0.08–0.12.
  - The accuracy twin (same orbits, accuracy 3×10⁻⁶) differs by only 20 in total, so integrator
    chaos is negligible.

**The ring averages (m = 0) are where even the true model fails.**
- Truth vs data χ² per coefficient: 2.59 at m = 0, 1.37 at m = 2, 1.20 at m = 4, then 1.02–1.13.
- The misfit is in the dispersions. Ring-average χ²/N_data is 3.98 for σ_R, 4.26 for σ_φ and 3.41
  for σ_z, against 1.03–1.48 for the means.
- The truth library is hotter than the data at almost every R and |z|: ring-average pulls are mostly
  −1 to −4.
- Part of it is a bug: data dispersions use 1/N, not 1/(N−1), while the library has thousands of
  samples per cell (median N = 57, bias σ/2N = 0.65 km/s against an error of 6.0). Correcting it
  lowers the m = 0 dispersion χ²/N to 3.28 / 3.15 / 2.72, so it explains about 20–25% of the excess.
- The rest is not yet explained. Untested candidates:
  - the spherical halo (F_z 2.5% low);
  - bias of the 1/τ̂ weights estimated from 40 snapshots;
  - residual non-equilibrium.
- **m = 0 drives the preference for weaker or shorter bars.** In m = 0 alone, f 0.7 and S 0.8 fit
  *better* than the truth (Δχ² −71, −44). In modes 2–4 alone the truth beats all of them: S 1.25 +171,
  S 0.8 +249, f 1.4 +105, f 0.7 +207, ridges +76 and +144. Truth χ²/N in modes 2–4 is 1.28.

**Posteriors** (GP emulator, 64 Sobol points + truth; emcee unseeded, so medians vary by about 0.1σ
between runs):

| likelihood | N_data | Ω_b | φ_b [°] | S | bar_frac | LOO | emulator max (S, f) vs evaluated truth |
|---|---|---|---|---|---|---|---|
| sectors (current) | 3090 | 26.26 +1.32/−1.13 | 27.7 ± 4.6 | 0.941 +0.229/−0.141 | 0.936 +0.403/−0.251 | 6.4 | (0.75, 1.56), −25.49 vs −25.38 |
| Fourier m ≤ 4 | 1434 | 26.39 +1.14/−1.01 | 27.3 ± 3.7 | 0.963 +0.293/−0.177 | 0.957 +0.460/−0.251 | 10.1 | (0.75, 1.60), −19.86 vs −20.69 |
| Fourier m = 2–4 only | 1146 | 26.24 +1.96/−1.79 | 26.3 ± 3.7 | 0.974 +0.201/−0.109 | 1.074 +0.316/−0.293 | 2.7 | (0.91, 1.40), −15.42 vs −15.37 |
| means m ≤ 4 + disp m = 2–4 | 1290 | 25.76 +1.32/−1.14 | 26.9 ± 3.7 | 0.972 +0.200/−0.127 | 1.012 +0.382/−0.294 | 3.9 | (0.83, 1.59), −16.74 vs −16.31 |

Corners: `figures/test1_nbody_emulator_corner_logLs_fourier_mmax4_meansdisp.png` (m ≤ 4),
`figures/test1_nbody_emulator_corner_logLs_fourier_m2to4_meansdisp.png` (m = 2–4) and
`figures/test1_nbody_emulator_corner_logLs_fourier_means_m024_sigmas_m24.png` (hybrid). These were regenerated
after the 1/N fix with a seeded emcee, so they match the numbers in the next entry, not this table.

- All truth values lie within 0.35σ in every variant.
- Dropping the dispersion ring averages centres S and bar_frac on the truth (medians 0.97 and
  1.01–1.07, against 0.94 and 0.94).
  - With m = 2–4 only, the emulator maximum moves off the box corner.
  - φ_b tightens from ±4.6° to ±3.7°.
- Ω_b needs the m = 0 means: without them it widens from ±1.2 to ±1.9 km/s/kpc. The hybrid restores
  ±1.2.
- **No variant breaks the S–bar_frac ridge.** Widths change by ≤ 25%, and in the hybrid the ridge
  still reaches (S 0.75–0.85, f 1.4–1.6).

## 2026-10-27 — 1/N dispersion bug fixed; per-field Fourier schemes

**Bug fix.**
- `ssl.cell_moments` now uses the unbiased (reliability-weights) variance: W²/(W² − Σw²) × the 1/W one,
  i.e. N/(N−1) for unit weights. The same fix was made in the data paths of
  `plot_bar_length_strength_square_maps.py`, `test1_binning_study.py` and `diagnose_moment_likelihood.py`.
- The library does not need recomputing. Its effective number of independent orbits per cell has a
  median of 1145 (minimum 128), so its relative bias is 4×10⁻⁴, against 9×10⁻³ for the data (median N 57).
- Rerun: `test1_nbody_moments.py` (all `chi2b_*`/`logLs_*` columns), `test1_fourier_compression.py`
  and the emulators below.
- Earlier entries (binning study, square maps, S/f diagnostics) still use 1/N.
- Effect at the truth:
  - means + dispersions χ²/N_data goes 1.29 → 1.24;
  - ring-average χ²/N goes 3.98 → 3.20 (σ_R), 4.26 → 3.05 (σ_φ), 3.41 → 2.66 (σ_z);
  - the truth is still the best of the 65 points, with the second best at −0.92 (was −0.83).

**Per field** (m_max = 4; Δχ² vs the truth for S 1.25 / S 0.8 / f 1.4 / f 0.7 / ridge 0.82,1.50 /
ridge 1.25,0.70; "m0" = ring average, "m2–4" = modes 2 and 4 with the ring average projected out):

| field | truth χ²/N m0 | truth χ²/N m2–4 | Δχ² m0 | Δχ² m2–4 |
|---|---|---|---|---|
| ⟨v_R⟩ | 1.35 | 1.03 | −3 / −1 / +3 / −3 / +2 / 0 | +36 / +53 / +12 / +24 / +12 / +19 |
| ⟨v_φ⟩ | 1.48 | 1.13 | +64 / −1 / +22 / +11 / −9 / +25 | +71 / +48 / +37 / +38 / +35 / +43 |
| ⟨v_z⟩ | 1.03 | 1.25 | 0 / +3 / +3 / +2 / +4 / +5 | −8 / +9 / −4 / +6 / +11 / −6 |
| σ_R | 3.20 | 1.35 | +47 / −4 / +58 / +1 / −5 / +17 | +42 / +27 / +31 / +27 / −5 / +21 |
| σ_φ | 3.05 | 1.49 | +93 / −13 / +111 / −34 / +41 / +27 | +28 / +69 / +21 / +68 / +4 / +50 |
| σ_z | 2.66 | 1.39 | +33 / −11 / +34 / −21 / +7 / −5 | +5 / +34 / +19 / +33 / +24 / +11 |

- The ring average of ⟨v_R⟩ carries no S/bar_frac information.
- ⟨v_z⟩ carries almost none.
- The σ_φ and σ_z ring averages are what favour the weaker or shorter bars.

**Schemes** (both cos and sin of every m kept):
- `fourier_vR_m24_vphi_sigR_sigphi_m024` (proposed): ⟨v_R⟩ m = 2, 4; ⟨v_φ⟩, σ_R, σ_φ m = 0, 2, 4. N_data 908.
- `fourier_vR_sigR_sigphi_m24_vphi_m024`: the same, but σ_R, σ_φ m = 2, 4 only. N_data 812.

Tempered ΔlogL of the S/f models:

| likelihood | S 1.25 | S 0.8 | f 1.4 | f 0.7 | ridge 0.82,1.50 | ridge 1.25,0.70 |
|---|---|---|---|---|---|---|
| sectors (current, fixed) | −2.77 | −1.22 | −2.20 | −0.98 | −0.72 | −1.25 |
| vR m24, vphi/σ_R/σ_φ m024 (proposed) | −4.46 | −2.11 | −3.42 | −1.59 | −0.87 | −2.35 |
| vR/σ_R/σ_φ m24, vphi m024 | −2.97 | −2.44 | −1.53 | −2.09 | −0.46 | −1.94 |

**Posteriors** (after the fix; emcee seeded from here on, so reruns reproduce; corner =
`figures/test1_nbody_emulator_corner_logLs_<column>.png`):

| likelihood (column `logLs_...`) | Ω_b | φ_b [°] | S | bar_frac | LOO |
|---|---|---|---|---|---|
| sectors, current (`meansdisp`) | 26.16 +1.37/−1.15 | 27.6 ± 4.6 | 0.976 +0.233/−0.162 | 0.912 +0.390/−0.237 | 6.1 |
| all six, m ≤ 4 (`fourier_mmax4_meansdisp`) | 26.25 +1.14/−1.00 | 27.1 ± 3.7 | 0.986 +0.287/−0.192 | 0.924 +0.490/−0.226 | 9.7 |
| all six, m = 2–4 (`fourier_m2to4_meansdisp`) | 26.19 +1.90/−1.83 | 26.1 ± 3.7 | 0.975 +0.187/−0.111 | 1.052 +0.333/−0.280 | 2.7 |
| means m ≤ 4 + all σ m = 2–4 (`fourier_means_m024_sigmas_m24`) | 25.81 +1.31/−1.22 | 26.9 ± 3.6 | 0.987 +0.200/−0.134 | 0.954 +0.416/−0.257 | 3.7 |
| proposed (`fourier_vR_m24_vphi_sigR_sigphi_m024`) | 26.26 +1.11/−0.96 | 27.6 ± 3.4 | 0.929 +0.312/−0.149 | 1.028 +0.431/−0.321 | **10.8** |
| proposed, σ ring averages dropped (`fourier_vR_sigR_sigphi_m24_vphi_m024`) | 25.66 +1.23/−1.08 | 27.2 ± 3.3 | 0.957 +0.194/−0.134 | 0.978 +0.432/−0.276 | 4.0 |

- **Proposed scheme:** the corner is bimodal along the ridge (bar_frac near 0.7 and 1.5, S piled at 0.75).
  The emulator LOO of 10.8 is ≫ the ΔlogL ≈ 1 scale, so that corner is not reliable.
  - Its maximum (S 0.75, f 1.6, −14.92) is an extrapolation above the evaluated truth (−15.57).
  - Of the evaluated points, the truth is the best; the runner-up is −0.72.
- All truth values lie within 0.31σ.

## 2026-10-27 — How much hotter is the ground-truth library? (ring-averaged dispersions)

- Script: `scripts/plot_ring_average_dispersion_residuals.py`.
- Figure: `figures/test1_nbody_ring_avg_dispersion_residuals_vs_R.png`. It plots model − data per 0.5 kpc ring
  and |z| layer, with the ±1σ data band.
- Numbers: `results/test1_nbody/ring_avg_dispersion_residuals.json`.
- "Ring average" is the m = 0 Fourier term of the χ²: the error-weighted mean of the 15° sectors with ≥ 20 stars.
  The data are unbiased (N/(N−1)). The |z| > 2 kpc layer has no usable sector.
- **Ground-truth library − data**, weighted over the 48 ring-layers of each dispersion:

  | | all R | R > 0.5 kpc | median relative | rings hotter |
  |---|---|---|---|---|
  | σ_R | +1.48 ± 0.25 km/s | +1.38 ± 0.25 km/s | +1.6% | 71% |
  | σ_φ | +1.65 ± 0.20 km/s | +1.64 ± 0.20 km/s | +2.4% | 88% |
  | σ_z | +1.00 ± 0.19 km/s | +0.85 ± 0.20 km/s | +1.9% | 71% |

- **The offset is 5–8σ when pooled, but below the error of any single ring** (1.6–1.9 km/s median).
  That is why it is invisible on `figures/test1_nbody_bar_length_strength_radial_profiles.png`, whose σ axes
  run from 40 to 200 km/s and whose data points are biased low by the 1/N bug.
- **S and bar_frac move the ring averages by only a fraction of a km/s**, model − truth mean:
  - σ_R: f 0.7 −0.60, S 0.8 −0.56, f 1.4 +0.68, S 1.25 +0.68 km/s;
  - σ_φ: −0.55, −0.16, +1.06, +0.51 km/s.
  So no single ring separates the models. The χ² sums 48 rings per dispersion.
- **The offset is what makes the ring averages prefer weaker or shorter bars.** A colder library sits closer
  to data that are already 1–1.6 km/s colder than the truth library. With the truth's offset subtracted
  (R > 0.5), the Δχ² of f 0.7 / S 0.8 becomes:
  - σ_R: −1 / −6 → +19 / +12;
  - σ_φ: −34 / −13 → +4 / 0;
  - σ_z: −16 / −14 → −11 / −11.
- The centre is not the main cause. R < 0.5 kpc is hotter by +6 to +11 km/s (the softening region), but
  excluding it barely changes the offset.
- **Excess scatter.** After removing the offset, χ²/N is still 2.44 (σ_R), 1.64 (σ_φ) and 2.10 (σ_z). Ring to
  ring, the residuals scatter 1.3–1.6× more than the bootstrap errors.
- The cause of the offset is not yet known. Untested candidates:
  - the spherical halo (F_z 2.5% low);
  - bias of 1/τ̂ from 40 snapshots;
  - residual non-equilibrium of the N-body.

**Renamed for clarity (2026-10-27).**
- "SF" meant bar length S and strength f; it read like "selection function".
- Scripts:
  - `plot_SF_square_maps.py` → `plot_bar_length_strength_square_maps.py`
  - `plot_SF_model_maps.py` → `plot_bar_length_strength_polar_maps.py`
  - `run_SF_box_models.py` → `run_bar_length_strength_models.py`
- Figures:
  - `test1_nbody_SF_maps` → `test1_nbody_bar_length_strength_polar_pull_maps`
  - `test1_nbody_SF_profiles` → `..._bar_length_strength_radial_profiles`
  - `test1_nbody_SFmaps_{values,nonaxi,diff,m2}` → `..._bar_length_strength_polar_maps_values`,
    `..._polar_maps_nonaxisymmetric`, `..._polar_maps_model_minus_truth`, `..._m2_profiles`
  - `test1_nbody_maps_{length,strength}` → `test1_nbody_square_maps_bar_{length,strength}_vs_data`
  - `test1_nbody_fourier_compression` → `test1_nbody_fourier_compression_vs_mmax`
- Fourier columns and corners: `logLs_F4_inplane` → `logLs_fourier_vR_m24_vphi_sigR_sigphi_m024`, etc.
  The links above were updated.

## 2026-10-27 — Do we need both cos and sin terms? (symmetry-allowed terms only)

- Script: `scripts/test1_fourier_symmetry_terms.py`; numbers in `results/test1_nbody/fourier_symmetry_terms.json`.
- Corner: `figures/test1_nbody_emulator_corner_logLs_fourier_symm_vR_sin24_vphi_sigR_sigphi_cos024.png`.
- **Symmetry argument.** In the bar frame, a steady bar response symmetric under (y → −y, t → −t) has ⟨v_R⟩ odd
  about the bar axis (sin mφ only) and ⟨v_φ⟩, σ_R, σ_φ even (cos mφ only).
  - "Symmetric" scheme: ⟨v_R⟩ sin 2φ, sin 4φ (ring average projected out); ⟨v_φ⟩, σ_R, σ_φ 1, cos 2φ, cos 4φ.
    N_data 528.
  - "Both": cos and sin of every m, as used so far. N_data 908.
  - "Forbidden": the rest. N_data 380.
- **The symmetry holds in the truth library.** The forbidden term carries 0.8% (⟨v_R⟩), 0.4% (⟨v_φ⟩), 0.3% (σ_R)
  and 0.8% (σ_φ) of the m = 2 power (|z| < 1).
- Δχ² vs the truth (tempered ΔlogL in brackets; other parameters at the truth):

  | offset | both | symmetric only | forbidden only |
  |---|---|---|---|
  | S 1.25 | +380 (−4.46) | +352 (−5.41) | +28 (−0.51) |
  | S 0.8 | +180 (−2.11) | +177 (−2.73) | +2 (−0.04) |
  | f 1.4 | +292 (−3.42) | +275 (−4.23) | +17 (−0.31) |
  | f 0.7 | +135 (−1.59) | +113 (−1.74) | +23 (−0.41) |
  | φ_b 17° | +530 (−6.22) | +30 (−0.46) | +500 (−9.07) |
  | φ_b 37° | +394 (−4.62) | +49 (−0.75) | +345 (−6.25) |
  | Ω_b 20 | +5893 (−69.2) | +5874 (−90.4) | +19 (−0.35) |
  | Ω_b 32 | +1168 (−13.7) | +1190 (−18.3) | −23 (+0.41) |

  - S, bar_frac and Ω_b are measured by the symmetric (amplitude) terms.
  - **The bar angle is measured almost only by the forbidden (phase) terms**, which carry 88–94% of its Δχ².
- Posteriors:
  - Symmetric only: Ω_b 25.99 +0.88/−0.70, **φ_b 25.4 +11.7/−7.9°**, S 1.050 +0.241/−0.268,
    bar_frac 0.856 +0.481/−0.179. Emulator LOO 16.1, so the emulator is unreliable here.
  - Both: φ_b 27.6 +3.4/−3.5°.
- **Conclusion.** The forbidden terms are needed whenever φ_b is fitted. They could be dropped only with φ_b fixed
  externally.
  - Also: distance errors rotate the observed ⟨v_R⟩ quadrupole (−25.5° → −6.2° for σ_d/d = 0.2, see the
    N-body entry). This creates forbidden-term power in the data, which the forward model has to reproduce.
- Figure index: add `test1_nbody_emulator_corner_logLs_fourier_symm_vR_sin24_vphi_sigR_sigphi_cos024.png`
  (posterior, symmetric terms only).
- **Check: where do the patterns point when φ_b is wrong?** (`scripts/check_m2_orientation_vs_trial_angle.py`)
  m = 2 orientation relative to the true angle, R 1.5–5 kpc, |z| < 1, in degrees (data / library):

  | trial φ_b | ⟨v_R⟩ | ⟨v_φ⟩ | σ_R | σ_φ |
  |---|---|---|---|---|
  | 17° | −6.9 / +0.4 | −12.6 / +0.6 | −11.1 / +0.5 | −15.7 / +1.2 |
  | 37° | +11.9 / +0.1 | +7.9 / +0.5 | +8.8 / −0.2 | +7.1 / −1.4 |

  - The library patterns always stay on the potential's bar axis (x), within ±1.4°.
  - The data patterns sit where the true bar is in the trial frame (expected −10° / +10°).
  - The phase terms measure exactly this data–model misalignment, which is why they carry the bar-angle
    information; the cos terms change only as cos 2δ.
- **Square maps for a wrong angle**: `figures/test1_nbody_square_maps_wrong_bar_angle_20deg_vs_27deg.png`
  (`scripts/plot_wrong_bar_angle_square_maps.py`). Stars rotated by 20° instead of 27° are integrated in the true
  potential (bar on x); the library is cached in `~/mira_dr4_runs/test1_nbody/maps_phi20.npz`.
  - The library's patterns stay on x: library 20° − library 27° has rms 2.0–3.0 km/s, at the library-noise level.
  - The data's patterns rotate (m = 2 orientation, data / library): ⟨v_R⟩ −4.8° / +0.9°, ⟨v_φ⟩ −6.9° / +0.3°,
    σ_R −6.6° / −0.1°, σ_φ −8.5° / +0.9° (expected data −7°).
  - Data − library rms at the wrong angle vs the true angle: 11.0 vs 11.0 (⟨v_R⟩), 10.4 vs 9.1 (⟨v_φ⟩),
    6.9 vs 6.5 (σ_R), 7.1 vs 6.8 (σ_φ) km/s. A coherent quadrupole residual is visible in ⟨v_φ⟩ at |y| ≈ 2–3 kpc.

## 2026-10-27 (overnight) — Focused global search with the Fourier R–|z| grid (set-up; results in the next entry)

**Integrations.**
- Script: `scripts/test1_nbody_qmc_focused_run.py`; log in `~/mira_dr4_runs/test1_nbody_qmc_focused/run.log`.
- Design (`results/test1_nbody_qmc_focused/design.csv`):
  - the truth and the 64 old Sobol points, reprocessed from their saved trajectories (no re-integration);
  - 192 new scrambled Sobol points (seed 2028) in a box focused on the posterior in Ω_b and φ_b and wider in S and
    bar_frac: |Ω_b| 20–32 km/s/kpc, φ_b 15–39°, S 0.6–1.65, bar_frac 0.45–2.0 (log-uniform). The old box was 18–38,
    14–44°, 0.75–1.45 and 0.55–1.6; the S–bar_frac ridge ran into its corner.
  - Integration as before: accuracy 10⁻⁵, about 90 s per point, so about 5.5 h; started 23:29.
- Saved per point (`~/mira_dr4_runs/test1_nbody_qmc_focused/lib_NNN.npz`): library sums on a fine folded grid
  (R 0.1 kpc to 8 × φ mod π 7.5° × |z| 0.25 kpc to 4 + one overflow bin) — N, W, Σw², Σw v, Σw v_a v_b. Any coarser
  (R, |z|, φ) grid with edges on these is an exact sum.
  - Check: these reproduce the saved polar moments of the truth to 3×10⁻¹² km/s.
  - Also saved: τ, the polar cell moments of the current likelihood, and the trajectories of the new points.

**Fourier grid** (`scripts/test1_nbody_qmc_focused_fourier.py`; per point
`~/mira_dr4_runs/test1_nbody_qmc_focused/fourier_NNN.npz`).
- Rings: R 0–6 kpc in 20 rings of 0.3 kpc (requested), extended by 4 rings of 0.5 kpc to 8 kpc.
- Layers: |z| 0–1, 1–2, 2–3 kpc (folded).
- Per block, a GLS fit of the sector values (data bootstrap errors, ≥ 20 stars per sector) on
  1, cos 2φ, sin 2φ, cos 4φ, sin 4φ for all six fields. Saved: data and library coefficients plus the 5×5 data
  covariance, so any subset of fields or terms can be scored exactly later.
- Star counts per ring-layer: |z| 0–1: 674–5538; |z| 1–2: 85–513; |z| 2–3: 7–60.
  - The |z| 2–3 layer (820 stars, 1.6%) is too sparse for 0.3 kpc rings and drops out.
  - Sectors are chosen once from the data at the true angle: 15° if ≥ 6 usable, else 30°, else 45° (45°: m = 0, 2 only).
  - |z| 0–1 uses 15° everywhere. |z| 1–2 uses 15° at R 0.6–3.0, 45° at 3.0–5.7, and 30° or 15° beyond.
- Columns in `results/test1_nbody_qmc_focused/evals.csv` (tempered `logLs_*`, N_data at the truth):
  - `fourier_R0to6_vR_m24_vphi_sigR_sigphi_m024` (proposed, N_data 669);
  - `fourier_R0to6_vR_sigR_sigphi_m24_vphi_m024` (no σ ring averages, 591);
  - `fourier_R0to6_all6_m024` (1062);
  - `fourier_R0to8_vR_m24_vphi_sigR_sigphi_m024` (821);
  - `sectors_means`, `sectors_meansdisp` (current likelihood).
- Truth χ²/N_data on the proposed scheme: 973/669 = 1.45.

**Fourier moments fitted directly to the stars (main method; replaces the sectors).**
- Script: `scripts/test1_nbody_qmc_focused_fourier_stars.py`; per point
  `~/mira_dr4_runs/test1_nbody_qmc_focused/fourier_stars_NNN.npz`.
- In each (ring, layer) block:
  - stage 1: least squares over the individual stars of v_i = Σ a_k f_k(φ_i), with f = 1, cos 2φ, sin 2φ, cos 4φ,
    sin 4φ, giving ⟨v⟩(φ);
  - stage 2: least squares of r_i²/(1 − h_ii) on the same basis (leverage-corrected; equals N/(N−1) for a constant),
    giving σ²(φ).
- The normal equations are sums over stars, so no sectors are needed. Library: all samples t ≥ 0.25 Gyr from the
  trajectories, weights 1/τ. Data covariance: 5×5 per block from 200 bootstraps. A block needs ≥ 50 stars.
- **Why a least-squares fit and not the plain sum (2/N) Σ v_i cos 2φ_i.** The plain sum measures the density-weighted
  field n⟨v⟩. For ⟨v_φ⟩ m = 2 cos at |z| < 1:
  - plain sum +50 to +84 km/s vs least squares −6 to −45 km/s at R = 1.8–4.8 kpc;
  - the plain sum is dominated by ⟨v_φ⟩₀ × density A₂ (A₂ = 0.45–0.77).
  - Fake selection, removing 60% of the stars within 30° of the bar axis, at R = 3.0–3.3: the plain sum goes
    +65.3 → −8.4 km/s, least squares −43.3 → −41.6 km/s.
- Usable blocks now: all of |z| 0–1 and |z| 1–2 (24 rings each), plus 5 rings at |z| 2–3.
- At the truth, proposed scheme: χ²/N_data = 1170/855 = 1.37 (sectors: 669 terms). Ridge point #15: Δχ² = +46.
- Columns (`logLs_fourier_stars_*` in `results/test1_nbody_qmc_focused/evals.csv`):
  `R0to6_vR_m24_vphi_sigR_sigphi_m024` (855), `R0to6_vR_sigR_sigphi_m24_vphi_m024` (765), `R0to6_all6_m024` (1350),
  `R0to8_vR_m24_vphi_sigR_sigphi_m024` (1007).

- **Why not divide the plain sum by Σ cos mφ_i.** For n = n₀(1 + A₂ cos 2φ) and ⟨v⟩ = v₀ + a₂ cos 2φ,
  Σ v cos 2φ / Σ cos 2φ ≈ v₀ + a₂/A₂: it mixes in the ring mean and diverges as A₂ → 0.
  - Measured for ⟨v_φ⟩ (|z| < 1): +40 to +285 km/s, against least-squares a₂ = −43 to +1 km/s.
  - The correct normalisation is the 5×5 matrix of the density's Fourier moments; least squares is exactly
    a = (Σ f fᵀ)⁻¹ Σ v f.
- **Means and dispersions are correlated** (bootstrap, R 3.0–3.3, |z| < 1):
  - ⟨v_φ⟩ vs ⟨v_φ²⟩ terms: +0.93 (m0), +0.96 (cos 2φ);
  - ⟨v_φ⟩ vs residual σ_φ²: −0.42 / −0.50, from the skewed v_φ distribution.
  - So the χ² now uses the **joint** bootstrap covariance of all kept terms of a block (300 bootstraps, Hartlap
    factor). The block-diagonal version is kept as `chi2_<scheme>_blockdiag`.
  - Truth χ²/N_data (proposed scheme): joint 1.50, block-diagonal 1.31. Ridge #15 Δχ²: +66 joint, +53 block-diagonal.
  - With the joint covariance, fitting ⟨v²⟩ or residual σ² carries the same information.
- **Slow run.** The new integrations take ~300 s each instead of ~90 s, because an unrelated process
  (`validate_kernel.py`, started from another session) uses ~4.7 cores.
  - The post-processing therefore runs on the points finished by 05:00 (`postprocess_partial.log`) and again when
    all points are done (`postprocess_final.log`).

**Post-processing.**
- Runs automatically (see above): `scripts/run_test1_nbody_qmc_focused_postprocess.sh`.
  It runs the sector version too, for comparison.
- Emulator corners use leave-one-out on the 32 best points; there are also log L corners. Emulator registered as run
  `nbody_qmc_focused`.
- Figures: `figures/test1_nbody_qmc_focused_emulator_corner_<column>.png` and
  `figures/test1_nbody_qmc_focused_logl_corner_<column>.png`.

## 2026-10-28 — Focused global search: results on the first 173 points (partial)

**Status at 09:15.**
- 196 of 257 integrations are done (the truth, the 64 old points and 131 new ones).
- The post-processing below used the 173 points finished at 05:12. A final pass runs automatically when the
  integrations end (log `~/mira_dr4_runs/test1_nbody_qmc_focused/postprocess_final.log`).
- **Why the run was late: the Mac went to sleep** (`pmset -g log`).
  - 23:35–00:18: ~300 s per point while `validate_kernel.py` (another session) shared the CPUs.
  - Then 150–175 s per point (01:14–01:41) and 83–107 s (03:08–03:24).
  - 03:25: "Entering Sleep state due to 'Idle Sleep'". After that, only dark wakes of 4–63 s every 5–18 min
    (several ended with "Dark Wake Thermal Emergency"), until a user wake at 08:58.
  - So points 171–187 show 324–2805 s "integration" times. The partial post-processing (started 05:12 in a dark wake)
    took until 08:57.
  - From 08:58: 81–101 s per point again.
  - I did not keep the machine awake. That is fixed now: `caffeinate -i -s -t 14400` started at 09:20.
  - (I first blamed CPU contention with the post-processing; the timings and the power log show that was wrong.)
- **Fix for the final pass.**
  - `test1_nbody_qmc_focused_fourier.py` and `..._fourier_stars.py` now re-use the points already analysed
    (`fourier_NNN.npz`, `fourier_stars_NNN.npz`, previous `evals.csv` rows); the final pass only does new points.
  - Check on a scratch copy: all 28 χ²/logL columns of the 173 points reproduced to 1.4×10⁻¹⁴ in 6 s.
  - `scripts/run_test1_nbody_qmc_focused_postprocess.sh` now also runs the control below.
- Slice and leave-one-out plots in `scripts/test1_gp_emulator.py` were fixed: y range −30…+5 for tempered
  likelihoods (was −3000…+300), LOO axes limited to the left-out points, LOO rms printed to 0.01.

**Evaluated points (172 Sobol + truth; tempered logL; ΔlogL relative to the truth).**

| likelihood (column `logLs_…`) | N_data | truth χ²/N | truth rank | best other (idx) | n within −2 / −4.5 |
|---|---|---|---|---|---|
| `fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024` (proposed, per-star, joint cov.) | 855 | 1.50 | 1 | −0.79 (#15) | 5 / 14 |
| `fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024` (no σ m = 0) | 765 | 1.27 | 1 | −0.56 (#15) | 10 / 21 |
| `fourier_stars_R0to6_all6_m024` | 1350 | 1.65 | 1 | −0.33 (#85) | 7 / 22 |
| `fourier_stars_R0to8_vR_m24_vphi_sigR_sigphi_m024` | 1007 | 1.50 | 1 | −0.66 (#15) | 5 / 14 |
| `fourier_R0to6_vR_m24_vphi_sigR_sigphi_m024` (proposed, sector GLS fit) | 669 | 1.45 | 1 | −0.54 (#15) | 6 / 13 |
| `sectors_meansdisp` (current polar-cell likelihood) | 3090 | 1.25 | 1 | −0.92 (#15) | 5 / 30 |

- #15 is the ridge point (Ω 26.77, φ 27.6°, S 0.82, f 1.50); #85 is (26.40, 29.5°, 0.78, 1.29).
- Block-diagonal versions of the per-star χ² give the same ranking (truth first, #15 second in all four).
- Proposed scheme, best evaluated point in bins of bar_frac (ΔlogL):
  - f 0.45–0.6: −4.84;
  - f 0.6–0.75: −1.08;
  - f 0.75–1.33: −1.20;
  - f 1.33–1.6: −0.79;
  - f 1.6–2.0: −6.09.
  - Without σ m = 0, the f 0.45–0.6 bin reaches −2.80, so that posterior is not closed at the low-f edge of the box.

**GP-emulator posteriors** (172 Sobol points; LOO on the 32 best; seeded emcee; reruns reproduce them exactly).

| likelihood | |Ω_b| | φ_b [°] | S | bar_frac | LOO rms | GP noise σ |
|---|---|---|---|---|---|---|
| per-star, proposed | 26.02 +0.93/−0.80 | 27.3 +2.9/−3.1 | 1.04 +0.22/−0.16 | 0.88 +0.33/−0.24 | 0.86 | 1.17 |
| per-star, no σ m = 0 | 25.67 +1.15/−0.96 | 27.5 +3.2/−3.3 | 1.03 +0.23/−0.18 | 0.95 +0.38/−0.30 | 0.33 | 0.72 |
| sector GLS, proposed (control) | 26.19 +0.81/−0.73 | 27.8 +2.8/−2.8 | 1.11 +0.21/−0.20 | 0.85 +0.33/−0.21 | 0.61 | 1.71 |
| current sector likelihood (means + disp) | 26.07 +1.30/−1.09 | 27.2 +4.3/−4.4 | 0.96 +0.21/−0.15 | 0.93 +0.41/−0.31 | 0.43 | 0.45 |
| *old run (65 points), sector GLS, proposed* | 26.26 +1.11/−0.96 | 27.6 ± 3.4 | 0.93 +0.31/−0.15 | 1.03 +0.43/−0.32 | 10.8* | — |

Truth: 26.02, 27°, 1, 1. All truth offsets are within 0.54σ. (*LOO on the old run was over all points, not the 32 best.)
"GP noise σ" is the square root of the fitted white-noise level, in logL units.

**What this shows.**
- **The surface is well behaved.**
  - In every scheme the truth is the best of the 173 evaluated points.
  - The per-star emulator maximum is at (26.20, 27.6°, 1.006, 0.985), within 0.2 km/s/kpc, 0.6°, 0.6% and 1.5% of the
    truth.
  - The slices through the truth peak at the truth in all four parameters
    (`figures/test1_nbody_qmc_focused_emulator_slices_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024.png`).
- **Tighter than the old run.** Same scheme, old 65-point run → this run: σ(Ω) 1.04 → 0.87 (sector GLS 0.77),
  σ(S) 0.23 → 0.19, σ(f) 0.38 → 0.28.
- **The per-star fit is not what tightened it.** On the same 172 points and the same R–|z| grid, the sector GLS
  version of the scheme gives the same widths within ~10% (Ω 0.77 vs 0.87, φ 2.8 vs 3.0, S 0.20 vs 0.19, f 0.27 vs
  0.28). The per-star fit's benefits are that it needs no sectors (all |z| 0–2 blocks and 5 at |z| 2–3 are usable)
  and that it does not depend on the azimuthal density (fake-selection test, previous entry). Its information
  content is the same.
- **Design vs grid: not yet separated** (corrected 2026-10-28; an earlier version of this entry credited the design
  alone).
  - Old → new changes both the design (65 → 173 points, box containing the ridge) and the grid
    (0.5 kpc polar rings, |z| 0/0.5/1/2/∞ → 0.3 kpc rings, |z| 0–1/1–2/2–3).
  - The polar-cell likelihood, identical in both runs, isolates the design. Its widths barely change:
    - Ω 1.26 → 1.20, φ 4.6° → 4.4°, S 0.20 → 0.18, f 0.31 → 0.36;
    - its emulator maximum leaves the box corner: (S, f) = (0.75, 1.56) → (0.91, 1.13).
  - The old Fourier corner was bimodal, with S piled at the box edge 0.75 and f peaks near 0.7 and 1.5. Its
    maximum (0.75, 1.60) was at the box corner.
  - Emulator accuracy near the peak, LOO on the same number of best points (old values recomputed today):
    - best 16: old 1.48, new 0.44;
    - best 32: old 1.90, new 0.86.
    - The 16 best old points span ΔlogL 0 to −12; the new ones span 0 to −3.9.
    - The old LOO of 10.8 was over all 64 points, so it is not comparable.
- **Proposed vs current sector likelihood.** σ(Ω) 0.87 vs 1.20 and σ(φ_b) 3.0° vs 4.4°; S and f are about the same
  (0.19 vs 0.18 and 0.28 vs 0.36).
- **The S–bar_frac degeneracy is still there.** Shorter and stronger, or longer and weaker; see the S–bar_frac panels
  of both corners. The ridge point #15 (S 0.82, f 1.50) is still second best in five of the six likelihoods.

**Caveats.**
- Only 5 evaluated points lie within ΔlogL −2 of the truth (proposed scheme), so the inner contours depend on the
  emulator.
- The fitted GP noise (σ 1.17 in logL for the proposed scheme) is larger than the 1σ step (0.5).
  - The individual evaluations scatter around the smooth surface by about one 1σ step.
  - Library sampling noise is a likely contributor; this has not been tested.
- One LOO outlier: evaluated −5.2, predicted −1.0.

**Figures** (`figures/`, `<column>` as in the table):
- `test1_nbody_qmc_focused_emulator_corner_<column>.png` — posterior corner;
- `test1_nbody_qmc_focused_logl_corner_<column>.png` — evaluated points (upper) and GP profile (lower);
- `test1_nbody_qmc_focused_emulator_slices_<column>.png` — 1D slices through the truth and the emulator maximum;
- `test1_nbody_qmc_focused_emulator_loo_<column>.png` — leave-one-out check.

## 2026-10-28 — Is the per-star Fourier fit right? Check against the true N-body maps

- Script: `scripts/check_fourier_star_fit_vs_true_maps.py`; numbers in
  `results/test1_nbody_qmc_focused/fourier_star_fit_vs_true_maps.json`.
- Figures: `figures/test1_nbody_fourier_star_fit_vs_true_maps_absz_0to1.png` (0.25 kpc pixels) and
  `figures/test1_nbody_fourier_star_fit_vs_true_maps_absz_1to2.png` (0.5 kpc pixels).
  - Columns: true map; m ≤ 4 reconstruction from the 50k mock; their difference in km/s and in units of the error;
    m ≤ 4 fit to all particles minus the truth.
  - Rows: ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ.
- **Ground truth** = all 1,493,610 star particles at R < 8 kpc. The mock is a random 50k subset of these (checked:
  same pool size). Bar frame at 27°.
- A pixel's reconstruction is evaluated at the true particles in that pixel, so the comparison is like for like.
- Re-fitting the mock in this script reproduces the saved coefficients exactly (max difference 0).

**1. The reconstruction agrees with the truth within its errors.**

| layer | field | mock − true rms [km/s] | pulls: mean, std | \|pull\| > 3 |
|---|---|---|---|---|
| \|z\| < 1 | ⟨v_R⟩ / ⟨v_φ⟩ / σ_R / σ_φ | 7.1 / 6.0 / 4.3 / 4.4 | +0.06, 1.07 / +0.01, 1.09 / +0.07, 0.93 / +0.11, 1.09 | 0.7 / 0.8 / 0.2 / 0.9% |
| 1 < \|z\| < 2 | same | 9.1 / 7.8 / 6.0 / 5.3 | +0.13, 1.03 / −0.01, 1.10 / +0.16, 0.94 / +0.11, 1.01 | ≤ 0.4% |

- No bias in the ring averages. Mean pull of the m = 0 coefficients (mock − true) over the 45 blocks:
  - means: ⟨v_R⟩ −0.00, ⟨v_φ⟩ −0.14, ⟨v_z⟩ −0.07;
  - variances: σ_R² +0.08, σ_φ² −0.26, σ_z² −0.02;
  - each ±0.13–0.16.
- The reconstructed maps show ring-to-ring noise (each 0.3 kpc ring is fitted independently); it is consistent with
  the errors (pulls above).

**2. m ≤ 4 truncation is small.**
- Fit to all particles − truth, rms against the truth-pixel noise alone (|z| < 1):
  - ⟨v_R⟩ 5.8 vs 5.2;
  - ⟨v_φ⟩ 4.9 vs 4.3;
  - σ_R 3.4 vs 3.7;
  - σ_φ 3.2 vs 3.0 km/s.
  - So ~2.4 km/s (in quadrature) of unresolved structure in the means, and ≤ 1.2 km/s in the dispersions.
- What is left out, as mode amplitudes of the true maps per ring (|z| < 1, R < 6; median and max; pool noise ≈ 0.9):
  - m = 1: 2.3 (4.5) for ⟨v_R⟩ and 2.1 (4.5) for ⟨v_φ⟩ — a lopsided part of the N-body. An even-m bar model cannot
    produce it, and the even-m fit ignores it by construction. It is visible as a top/bottom pattern in the ⟨v_φ⟩
    truncation panel.
  - m = 6: 1.5 (2.9) and 1.6 (4.2).
  - m = 3 and m = 8: ≈ 1, i.e. at the noise level.
  - All are ≲ 0.5× the mock's per-ring coefficient error (~5 km/s).

**3–4. The joint bootstrap covariance is miscalibrated in sparse blocks.**
- Test: 200 independent 50k mocks drawn from the same particles give the true sampling covariance.
- Proposed scheme, χ²/N of mock − true coefficients (expected 0.97 with a correct covariance):

  | layer (stars per block) | real mock, bootstrap cov | real mock, true cov | 200 independent mocks, bootstrap cov (mean, 16–84%) |
  |---|---|---|---|
  | \|z\| 0–1 (674–5538) | 1.09 | 1.07 | 1.02 (0.94–1.12) |
  | \|z\| 1–2 (85–513) | 1.49 | 1.02 | 1.25 (1.12–1.37) |
  | \|z\| 2–3 (50–60) | 1.52 | 0.64 | 2.19 (1.56–2.90) |
  | all R < 6 | **1.31** | **1.00** | — |

- Diagonal variances (bootstrap / true):
  - 0.96–1.04 at |z| < 1;
  - 0.92–1.10 at |z| 1–2;
  - 0.69–1.02 at |z| 2–3.
- So at |z| 1–2 the excess comes from the estimated correlations. The joint covariance (19 terms) is estimated from
  85–513 stars, and its inverse is too large. Hartlap only corrects for the 300 bootstraps, not for the few stars.
- Same independent mocks:
  - per-field covariance: 0.99 / 1.13 / 1.75;
  - diagonal only: 0.95 / 1.01 / 1.24.
  - Diagonal-only is calibrated on average but ignores the mean–dispersion correlations (−0.4 to −0.5), so it
    double-counts.
- **Consequence 1: the |z| ≥ 1 blocks currently get ~25–50% too much weight in the likelihood.**
- **Consequence 2: the truth's χ²/N = 1.50 is mostly the error model, not the model.**
  - Truth library vs mock 1.50 ≈ mock vs true N-body 1.31 + truth library vs true N-body 0.18 (all with the
    bootstrap covariance).
  - With the true covariance: mock vs true 1.00, truth library vs true 0.16.
  - The library's misfit is concentrated in the m = 0 dispersion terms (σ_R² 0.56, σ_z² 0.50, σ_φ² 0.40 per term,
    against 0.07–0.13 for m = 2, 4). That is the "hot library".
**Why the |z| 1–2 ⟨v_R⟩ reconstruction from the mock looks broken (2026-10-28, user question).**
- Two map columns were added to the figures:
  - the 50k mock binned directly in the same pixels: median 4 stars per 0.5 kpc pixel at |z| 1–2, and 5 per
    0.25 kpc pixel at |z| < 1;
  - the m ≤ 4 reconstruction from all 1.49 M particles.
- From all particles the quadrupole comes out cleanly at both layers. So the method and the phases are right; the
  mock is noise-limited.
- Signal-to-noise per ring (R < 6; true amplitude / mock 1σ error; median, range):

  | term | \|z\| 0–1 | \|z\| 1–2 | \|z\| 2–3 |
  |---|---|---|---|
  | ⟨v_R⟩ m = 2 | 3.9 (0.6–4.8) | **1.1 (0.2–1.7)** | 0.3 |
  | ⟨v_φ⟩ m = 2 | 3.8 (0.4–11) | 1.5 (0.2–7.0) | 0.2 |
  | σ_R² m = 0 | 23 (17–60) | 10 (5.9–18) | 3.7 |
  | σ_R² m = 2 | 6.5 (1.1–8.7) | 1.0 (0.1–4.0) | 0.4 |
  | σ_φ² m = 2 | 3.7 (0.4–11) | 1.9 (0.3–3.1) | 0.6 |

- At |z| 1–2 each 0.3 kpc ring's ⟨v_R⟩ quadrupole is measured at about 1σ, so the per-ring reconstruction is mostly
  noise. Because each ring is fitted independently, that noise shows up as arcs.
- The dispersion maps look good because their contrast is the m = 0 radial profile, measured at S/N 10–23 per ring.
  Their m = 2 parts are as noisy as ⟨v_R⟩.
- The pulls (std 0.94–1.10) confirm these are noise of the stated size, not a bug.
- Options, not yet done:
  - bigger blocks at |z| ≥ 1 (e.g. ~1 kpc rings at |z| 1–2 for ≳ 700 stars per block, drop |z| 2–3);
  - a model-based covariance (mock-sized draws from a reference library, one random snapshot per orbit);
  - shrinking the joint covariance towards per-field.

## 2026-10-28 — Focused global search: final results on all 257 points

- All 257 integrations finished at 11:03. The final post-processing (`postprocess_final.log`) re-used the 173 points
  already analysed.
- Truth still ranks first on the proposed per-star scheme.
  - The best other point is now #235 (Ω 25.71, φ 27.7°, S 1.17, f 0.69) at ΔlogL −0.13.
  - Then #245 (26.04, 28.1°, 1.32, 0.59) at −0.73, and #15 (26.77, 27.6°, 0.82, 1.50) at −0.80.
  - 9 points lie within −2 and 20 within −4.5.

| likelihood (256 Sobol points) | \|Ω_b\| | φ_b [°] | S | bar_frac | LOO (32 best) |
|---|---|---|---|---|---|
| per-star, proposed | 26.06 +0.86/−0.78 | 27.3 +3.0/−3.1 | 1.06 +0.25/−0.17 | 0.87 +0.31/−0.25 | 0.45 |
| per-star, no σ m = 0 | 25.61 +1.17/−0.95 | 27.8 ± 3.3 | 1.01 +0.23/−0.16 | 0.94 +0.38/−0.31 | 0.33 |
| sector GLS, proposed | 26.18 +0.81/−0.76 | 27.8 ± 2.8 | 1.13 ± 0.21 | 0.84 +0.34/−0.21 | 0.57 |
| polar-cell, means + disp | 26.07 +1.25/−1.05 | 27.3 ± 4.3 | 0.98 +0.22/−0.16 | 0.92 +0.41/−0.29 | 0.36 |

- These are within 0.05σ of the 173-point posteriors; the extra 84 points mainly improved the emulator (LOO
  0.86 → 0.45).
- These still use the bootstrap covariance that over-weights |z| ≥ 1 (see the ground-truth check).

## 2026-10-28 — Why is Ω_b well constrained but S and bar_frac are not?

- Script: `scripts/analyse_bar_length_strength_degeneracy.py`.
  - It uses the same emulator and seeds as the corner. Run on the 173-point table, it reproduces the reported
    posterior exactly.
- Numbers: `results/test1_nbody_qmc_focused/bar_length_strength_degeneracy_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024.json`.
- Figure: `figures/test1_nbody_qmc_focused_bar_length_strength_degeneracy_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024.png`.
  - Left: posterior in (S, bar_frac) with curves of fixed bar m = 2 potential.
  - Middle: fractional width vs R.
  - Right: Q_T(R) posterior.
- All on the 256 Sobol points, proposed per-star scheme.

**1. S and bar_frac trade off: corr(ln S, ln f) = −0.65.**

| parameter | marginal 1σ | conditional 1σ (others fixed at the max) | ratio |
|---|---|---|---|
| \|Ω_b\| | 0.82 km/s/kpc (3.1%) | 0.58 (2.2%) | 1.42 |
| φ_b | 3.0° | 3.1° | 0.99 |
| S | 19% | 9.7% | 1.97 |
| bar_frac | 32% | 17% | 1.90 |

- The degeneracy doubles the S and bar_frac widths.
- Best-measured combination: **bar_frac · S^2.06, to 12.7%**. The orthogonal direction is uncertain to 33%.
- The evaluated points (not the emulator) lie on the same ridge. The best four after the truth have
  f · S^2.06 = 0.95, 1.04, 1.00, 0.95 (#235, #245, #15, #42).
- Ω_b correlates with bar_frac (+0.49) but not with S or φ_b (|r| ≤ 0.04). That is why its marginal width is 1.4× its
  conditional one.

**2. Even with the others fixed, the likelihood is much less sensitive to S and bar_frac than to Ω_b.**
- From the conditional widths, a 10% change gives a tempered ΔlogL of about:
  - −10 for Ω_b;
  - −0.5 for S;
  - −0.2 for bar_frac.

**3. The data do not pin the bar force at any single radius** (my first guess, now ruled out).
- The bar part of a trial potential is bar_frac · Φ_bar(r/S). The AGAMA scale convention was checked to 6×10⁻⁶.
- So Q_T(R) = 2 bar_frac A₂₀(R/S) / (R |dΦ₀/dR|). Its posterior width is ≥ 24% at every R (minimum at R = 4.0 kpc),
  worse than S alone.
- The ridge f ∝ S^−2.06 has the slope of the fixed-potential curve at R ≈ 5–5.5 kpc near S = 1. But those curves bend
  (slope −3.1 to −1.4 over S 0.8–1.6), so no single-radius force is constant along the ridge.
- True Q_T peaks at 0.10 at R = 2.5 kpc. Posterior at 4 kpc: 0.079 +0.020/−0.017 (truth 0.087).

**Context.**
- Tempering inflates every width by (2 N_data)^¼ = 6.4×; the ratios between parameters are unaffected.
- Posterior width relative to the flat prior's in the box: Ω_b 4.2×, φ_b 2.3×, S 1.5×, bar_frac 1.4× narrower.
- Untested idea for the physics: the bar's non-axisymmetric force is at most 10% of the radial force (Q_T ≤ 0.10), so a
  17% change in bar_frac changes the in-plane force by ≤ 1.7%, while Ω_b enters every star's motion in the bar frame.
- A local Fisher analysis would test this. It would use the library coefficients' derivatives with respect to the
  parameters, estimated from the Sobol points near the truth, and show which fields, terms and radii carry the
  information on each parameter.

## 2026-10-28 — Data vs model Fourier coefficients (do the models fit?)

- Script: `scripts/plot_fourier_coefficients_data_vs_models.py`; numbers in
  `results/test1_nbody_qmc_focused/fourier_coefficients_data_vs_models.json`.
- Figures: `figures/test1_nbody_qmc_focused_fourier_coefficients_data_vs_models_absz_0to1.png` and `..._absz_1to2.png`.
  - Rows: ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ. Columns: m = 0, cos 2φ, sin 2φ, cos 4φ, sin 4φ, per 0.3 kpc ring.
  - Data: the 50k mock with bootstrap errors.
  - Models:
    - #0, the truth (best of 257);
    - ridge points #235 (S 1.17, f 0.69) and #15 (S 0.82, f 1.50);
    - #17, rejected: Ω 33.4, near truth otherwise, ΔlogL −26.
  - Grey: the same fit to all 1.49 M N-body particles.
  - σ terms are shown as σ₀ = √c₀ and c_m/2σ₀ in km/s.
- **The truth model fits.** Mean per-panel χ²/N (diagonal errors, R < 6, the 19 panels in the likelihood): 1.11 at
  |z| < 1 and 1.08 at |z| 1–2. It tracks both the data and the N-body line in every term.
- **The symmetry holds in the data.** ⟨v_R⟩ cos terms and the ⟨v_φ⟩/σ sin terms scatter around zero. The bar signal
  is in ⟨v_R⟩ sin 2φ (peak ≈ 22 km/s at R ≈ 3.5 kpc) and in the cos 2φ / cos 4φ terms of ⟨v_φ⟩, σ_R, σ_φ.
- **The rejected pattern speed fails visibly**: mean χ²/N 3.65 (|z| < 1) and 1.73 (|z| 1–2). Worst panels at |z| < 1:
  - ⟨v_φ⟩ cos 2φ: 12.9 (its trough is −22 km/s, against −45 in the data);
  - σ_R cos 2φ: 10.2 (≈ 12 against ≈ 24 km/s);
  - σ_φ cos 2φ: 9.4;
  - ⟨v_φ⟩ m = 0: 8.8 (rotation curve ~10 km/s too high at R 1–4);
  - ⟨v_R⟩ sin 2φ peak: ≈ 10 against 22 km/s.
- **The ridge models are almost indistinguishable from the truth**: mean χ²/N 1.18 and 1.18 at |z| < 1, and 1.07 and
  1.06 at |z| 1–2.
  - Their departures are within the error bars: ⟨v_R⟩ sin 2φ peak ≈ 18–20 against 22, and ⟨v_φ⟩ cos 2φ / cos 4φ near
    R ≈ 1–2 and 4 kpc.
  - In the σ ring averages they do *better* than the truth: σ_R m = 0 1.51 (#15) vs 1.87; σ_φ m = 0 0.84 (#235)
    vs 1.18. That is the hot-library effect.
  - This is the S–bar_frac degeneracy seen directly in the data space.
- The data are plotted at 27°. The models' own angles move the data coefficients by ≤ 0.29σ (#235, #15) and 0.49σ
  (#17), and the per-panel χ² uses each model's own data.

## 2026-10-28 — Test 2b: Mira-like distance errors (mock and first look)

**Error level from the DR3 Mira catalogue** (`data/gaia_dr3_lpv_mira.csv`, column `DM_err`).
- With the usual cut σ_DM < 0.4: 54,573 of 55,074 stars (99.1%).
- **Median σ_DM = 0.2654 mag** (16–84%: 0.233–0.292). That is a fractional distance error of
  ln(10)/5 × 0.2654 = **12.2%**.
- The catalogue's median `luminosity_distance_err / luminosity_distance` is also 0.1222.
- R < 8 kpc subset: 0.2645 mag (12.2%). Stars with a Gaia RV: 0.2115 mag.

**Mock.**
- Script: `scripts/make_mock_distance_error.py` (new helper `observe.scatter_distance_modulus`).
- File: `mock/nbody_mock_distance_error_sigmaDM0.265_R8_N50k.hdf5`.
- Same 50k stars as the perfect mock (selected on true R < 8 kpc).
- Error model: DM_obs = DM_true + 0.2654 N(0, 1) for every star (log-normal in distance); l, b, proper motions and
  v_los are exact.
- Observed x..vz, d are stored under the usual names, and true values as `<name>_true`. Seed 2029.
- Converting back at the true distance reproduces the perfect mock to 6×10⁻¹³.
- Size of the smearing:
  - d_obs / d_true − 1: std 0.124 (16–84%: −0.114 to +0.129);
  - 3D position shift: median 0.66 kpc (84%: 1.48);
  - |R_obs − R_true|: median 0.37 kpc, more than one 0.3 kpc ring;
  - observed R > 8 kpc: 2.7%.

**What the errors do to the Fourier coefficients.**
- Script: `scripts/plot_distance_error_effect_on_fourier_coefficients.py`; numbers in
  `results/test1_nbody_qmc_focused/distance_error_effect_on_fourier_coefficients.json`.
- Figures: `figures/test1_nbody_distance_error_fourier_coefficients_perfect_vs_smeared_absz_{0to1,1to2}.png`.
- χ²/N against the truth model (proposed scheme, joint covariance):
  - perfect data 1.50;
  - **smeared data 4.61** with no correction;
  - **1.41** when the truth library is "observed" with the same errors.
- To observe the library: each sample with true R < 8 is rotated to the MW frame, seen from the Sun, given
  σ_DM = 0.265 and converted back; samples landing at R > 8 drop out, as the data do.
- Without correction, per field and term (smeared / forward-modelled):
  - σ_φ² m = 0: 20.0 / 1.74;
  - σ_z² m = 0: 17.4 / 1.62;
  - σ_R² m = 0: 10.6 / 1.45;
  - ⟨v_φ⟩ m = 0: 7.4 / 1.06 (rotation curve ~10 km/s low at R 1–4).
- The smearing also creates terms that the bar's symmetry forbids, because it runs along lines of sight from the Sun,
  27° from the bar axis:
  - ⟨v_R⟩ cos 2φ down to −20 km/s at R 1–2;
  - σ_R and σ_φ sin 2φ and sin 4φ of +10 to +25 km/s.
  - The forward-modelled truth reproduces all of them.
- Data errors grow by 5% (⟨v_R⟩) to 47–55% (σ_φ², σ_z²).
- **Caveat.** This forward model starts from the true positions (the truth library). Real data give only smeared
  initial conditions, so the orbits themselves start from a smeared, non-equilibrium distribution. That part is not
  yet tested.

## 2026-10-28 — Test 2b: orbits from distance-smeared stars (one orbit per star fails)

- Script: `scripts/test2b_smeared_initial_conditions.py`; numbers in `results/test2b/smeared_initial_conditions.json`.
- Figures: `figures/test2b_smeared_initial_conditions_fourier_coefficients_absz_{0to1,1to2}.png`.
- Set-up:
  - The smeared mock (σ_DM = 0.265) was integrated in the true potential at the true Ω_b and φ_b: one orbit per star,
    accuracy 10⁻⁵, 1 Gyr, 65 s.
  - Trajectory: `~/mira_dr4_runs/test2b/traj_smeared_ics_true_potential.npy`.
  - Everything is compared with the smeared data (per-star Fourier fit, 300 bootstraps).

| library (true potential) | χ²/N, proposed scheme | σ_R² m = 0 | σ_φ² m = 0 | σ_R² m = 2, 4 | ⟨v_φ⟩ m = 0 |
|---|---|---|---|---|---|
| (a) smeared initial conditions, as is | **11.99** | 43.1 | 20.0 | 8.4 | 7.0 |
| (b) smeared initial conditions, observed again | **13.43** | 58.5 | 47.3 | 3.4 | 10.0 |
| (c) true initial conditions, observed (target) | **1.41** | 1.4 | 1.7 | 0.96 | 1.06 |
| (d) true initial conditions, as is (no correction) | 4.61 | 10.6 | 20.0 | 4.5 | 7.4 |

- **(a) is worse than doing nothing (d).** Started from the smeared stars, the library phase-mixes in the bar potential
  into a bar-symmetric state.
  - The terms the symmetry forbids, which the data keep, vanish: σ_R and σ_φ sin 2φ ≈ 0 against +15 to +25 km/s, and
    ⟨v_R⟩ cos 2φ.
  - The bar response weakens: σ_R cos 2φ peaks at ≈ 13 against ≈ 25 km/s, and ⟨v_R⟩ sin 2φ at ≈ 12 against ≈ 20.
  - It is hotter than the data in σ_R at R > 2.5 kpc (by 5–10 km/s).
- **(b) puts the forbidden terms back but applies the errors twice.** The ring-averaged dispersions are far too hot.
- **The drift is fast.** Proposed-scheme χ²/N of (a) by time window:
  - 8.2 in 0–0.1 Gyr;
  - 11.7 in 0.25–0.5 Gyr;
  - 13.2 in 0.5–0.75 Gyr;
  - 13.4 in 0.75–1 Gyr.
  - For (b): 15.8, 13.2, 15.2, 15.7.
  - At t = 0, (a) is the data themselves, so the smeared stars leave their starting distribution within tens of Myr.
- **Conclusion.** One orbit per smeared star does not work: the smeared distribution is far from a steady state of the
  true potential. The orbits need initial conditions whose distribution matches the true one; then (c) shows that
  observing the library with the errors fixes the comparison.
- Candidates, not yet run:
  - several orbits per star, with distances drawn from each star's distance posterior (a prior along the line of
    sight is needed);
  - iterating with the library's own density as that prior;
  - a library that does not start from the data stars.

## 2026-10-28 — Test 2b: posterior distance draws as initial conditions (best case; does not solve it)

- Script: `scripts/test2b_posterior_distance_draws.py`; numbers in `results/test2b/posterior_distance_draws.json`.
- Figure: `figures/test2b_posterior_distance_draws_fourier_coefficients_absz_0to1.png`.
- **Method.** Five distances per smeared star, drawn from its distance posterior on a ±4σ grid in DM. The prior is the
  true N-body: 1.44 M star particles at R < 8, without the mock's own 50k stars; k = 64 nearest-neighbour estimates.
  - A, density prior: likelihood × ρ(x(d)) d².
  - B, full phase-space prior: likelihood × f(x(d), v(d)) d⁴, using the measured proper motions and v_los; 6D metric
    40 km/s ≙ 1 kpc. Exact draws from B would reproduce the true distribution.
  - 250k orbits each in the true potential: 427 s (A) and 402 s (B). Trajectories in `~/mira_dr4_runs/test2b/`.
- Distances: the posterior means are closer to the truth than the observed distances (median |error| 6.0% for A and
  5.4% for B, against 8.2%). Single draws scatter like the errors, as they should.

| initial conditions | t = 0 draws vs perfect data, χ²/N | library observed vs smeared data, χ²/N | by window (0–0.1 / 0.25–0.5 / 0.5–0.75 / 0.75–1 Gyr) |
|---|---|---|---|
| smeared stars, 1 per star | 15.6 | 13.4 | 15.8 / 13.2 / 15.2 / 15.7 |
| draws A, density prior | **19.7** | **8.7** | 12.8 / 8.9 / 9.1 / 9.5 |
| draws B, full phase-space prior | **4.5** | **4.1** | 3.3 / 4.0 / 3.9 / 4.0 |
| true positions (target) | — | 1.41 | 1.57 / 1.52 / 1.48 / 1.40 |

- **Even the best case fails.** B leaves χ²/N ≈ 4, against 1.4 for the true positions.
  - The m = 2, 4 terms are close (1.0–1.3 per field).
  - The ring averages are too hot: σ_R² m = 0 10.6, σ_φ² m = 0 12.6, ⟨v_φ⟩ m = 0 3.7. σ_R and σ_φ are about 5 km/s
    high at R > 2.5 kpc, and the σ_R cos 2φ bar response is weaker.
  - Likely reason (not tested): the k = 64 6D estimate of f is too smooth to pin each star's distance, so draws at wrong
    distances get wrong tangential velocities.
- **A, the version we could use with a Galaxy density model, is worse at t = 0 than the smeared stars themselves**
  (19.7 vs 15.6). It ignores how f depends on velocity.
- **Conclusion.** De-smearing the initial conditions star by star needs an accurate 6D DF prior, which is the thing we
  are modelling. Even with the true N-body it leaves χ²/N ≈ 4. As the user expected, this does not solve the problem.

## 2026-10-28 — Test 2b: choosing among candidate distances by fitting the kinematic field (no discrimination)

- Script: `scripts/test2b_distance_draw_selection.py`; numbers in `results/test2b/distance_draw_selection.json`; log
  `~/mira_dr4_runs/test2b/draws_K10_fit.log`.
- **Idea (user):** integrate several candidate distances per smeared star and choose among them so that the observed
  model reproduces the observed field.
- **Set-up.**
  - Candidates: K = 10 per star at fixed quantiles of the DM error (z = ±0.13 … ±1.64σ). 500k orbits per trial
    potential, 650–750 s each. Potentials: #0 (truth), #17 (|Ω_b| 33.4; perfect data reject it at ΔlogL −26) and #15
    (S 0.82, f 1.50; ridge).
  - Library: every sample (t ≥ 0.25 Gyr, true R < 8, 1/τ) is observed with the distance errors.
  - Stored per (orbit, block): the 35 sums of G and of H for ⟨v_R⟩, ⟨v_φ⟩, ⟨v_R²⟩, ⟨v_φ²⟩.
  - χ²: our likelihood on c = G⁻¹H (19 terms per block, joint bootstrap covariance with Hartlap). 44 blocks,
    N = 836.
- **Fits (no regularisation).**
  - pick1 (user's proposal): one candidate per star; greedy coordinate descent from a random start with exact local
    χ² updates.
  - simplex: weights summing to 1 per star; softmax + L-BFGS with the analytic gradient.
  - free: positive weights (Schwarzschild-like).
- **Checks before the run:**
  - ⟨v⟩ coefficients equal the likelihood code's (difference 0);
  - analytic gradient matches finite differences (relative ≤ 2×10⁻⁶);
  - greedy running χ² equals a full recomputation.

| χ²/N | #0 truth | #17 wrong Ω_b | #15 ridge |
|---|---|---|---|
| uniform weights (all 10 candidates) | 32.4 | 31.3 | 30.7 |
| pick1: random start → final | 35.6 → **0.0002** (22 sweeps) | 31.4 → **0.0002** (16) | 31.9 → **0.0001** (9) |
| simplex | **0.0000** (159 iterations) | **0.0000** (126) | **0.0000** (100) |
| free | **0.0001** (300) | **0.0004** (300) | **0.0000** (300) |
| reference: truth library from true positions | 1.62 | — | — |

- **Every method fits every potential to χ²/N ~ 10⁻⁴**, including a pattern speed that perfect data reject at
  ΔlogL −26. So there is no discriminating power. This confirms the capacity argument: 50k choices against 836 numbers.
  Without fitting, the uniform-weight library does not separate them either (truth 32.4 vs 31.3 and 30.7).
- **The chosen distances are not the true ones.**
  - pick1: median |d/d_true − 1| = 9.5% (#0), 9.4% (#17), 9.4% (#15). That is worse than the observed 8.2%, although
    the nearest candidate is within 1.1% of the truth (median).
  - pick1 slightly prefers the nearest candidate (−1.64σ: 6,247 stars at #0, against 4,600–5,500 for the others).
  - simplex: weighted-mean errors of 7.1%, 7.7% and 7.2%; 4.6, 2.5 and 4.4 effective candidates per star.
- **Conclusion.** Free per-star choices among candidate distances have far more freedom than the moment data constrain.
  Any useful version needs far fewer free parameters than data terms (e.g. weights as a smooth function of orbit
  properties), or many more constraints.
- The per-orbit sums are saved (`~/mira_dr4_runs/test2b/draws_K10_features_NNN.npz`; trajectories
  `draws_K10_traj_NNN.npy`), so other weighting schemes on these libraries need no new integrations.

## 2026-10-28 — Test 2b: how large a distance error can the current method take?

- Scripts: `scripts/test2b_distance_error_truth_only.py` (cheap check, used) and `scripts/test2b_distance_error_tolerance.py`
  (bias scan with ± offsets; stopped after level 0 and three 2.5% points; trajectories kept in
  `~/mira_dr4_runs/test2b/tolerance/`).
- Method: the current method unchanged (orbits from the observed stars, proposed Fourier scheme, joint bootstrap
  covariance), in the true potential. Every level uses the same normal deviates, scaled.

| fractional distance error (σ_DM) | χ²/N at the truth | extra χ² vs perfect data | tempered ΔlogL (1σ = −0.5) |
|---|---|---|---|
| 0% | 1.50 | 0 | 0 |
| 2.5% (0.054 mag) | 1.64 | +89 | **−1.08** |
| 5% (0.109 mag) | 2.69 | +1019 | **−12.3** |
| 12.2% (0.265 mag, DR3 Miras) | 11.99 | ≈ +9000 | ≈ −110 |

- The extra misfit rises steeply: from 2.5% to 5% it grows by 11.5×, about the power 3.5 of the error.
- It reaches the 1σ step of the tempered likelihood (Δχ² = √(2N) ≈ 41) at about **2% distance error
  (σ_DM ≈ 0.04 mag)**, by interpolation.
- Above that, the misfit the distance errors cause at the truth exceeds what a 1σ change of the potential produces.
- This measures the size of the damage, not its direction (bias vs broadening). The bias scan was stopped because the
  machine was shared and the threshold was already clear.
- Perfect-data baseline of the bias scan (level 0, ± offsets): shifts of the maximum from the data noise alone were
  Ω_b +0.30 km/s/kpc (0.46σ), φ_b +0.2° (0.06σ), ln S −0.03 (−0.23σ) and ln bar_frac −0.07 (−0.45σ).
- Side result: the library from the true positions against the 2.5% smeared data gives χ²/N 1.36, lower than 1.50 on
  perfect data, because the smearing heats the data towards the hot library.
- **Conclusion.** With one orbit per observed star, the method needs distances to about 2% (≈ 0.04 mag). DR3 Mira
  distances (12.2%) are about 6× too uncertain.

## 2026-10-28 — Test 2b: are the bootstrap errors still right with distance errors? (7.5%)

- Script: `scripts/check_bootstrap_errors_with_distance_errors.py`; numbers in
  `results/test2b/bootstrap_errors_with_distance_errors.json`.
- Method:
  - Truth: 200 independent 50k mocks drawn from the 1.49 M N-body particles, each with its own σ_DM = 0.163 mag
    errors, fitted like the data at 27°.
  - Population value: all particles, smeared once.
  - Compared with the 300-bootstrap covariance of the 7.5% mock.
- Why the bootstrap should work: each smeared star is an independent draw from the population convolved with the
  error law, so resampling the smeared stars measures sampling noise and distance-error noise together.
- Variance ratio bootstrap/true (median per field):
  - |z| < 1: 0.99–1.04;
  - |z| 1–2: 0.99–1.07;
  - |z| 2–3 (4 blocks, 50–60 stars): 0.58–0.93.
- χ²/N of the mock against the population (proposed scheme):
  - bootstrap covariance: 1.19 / 1.29 / 1.61 for the three layers, 1.28 overall;
  - true covariance: 1.16 / 1.05 / 0.66, 1.06 overall.
  - This is the same sparse-block correlation problem as on perfect data; the distance errors don't cause it.
- Distance errors inflate the coefficient errors (|z| < 1, relative to the perfect mock) by:
  - ⟨v_R⟩ 1.00×, ⟨v_φ⟩ 1.06×, σ_R² 1.08×, σ_φ² 1.20×, σ_z² 1.19×.
- So the data errors include the distance errors correctly. The truth's χ²/N of 5.3 at 7.5% is model misfit (orbits
  from the observed stars), not underestimated errors. These errors cover the data only, not the library's dependence
  on the noisy initial conditions.
- **The "real observation" recipe double-counts** (user's proposal: a Monte Carlo over each star's distance error,
  300 re-draws of DM_obs + σ_DM N(0, 1), co-added to the bootstrap). Variance ratio against the true scatter
  (|z| < 1; |z| 1–2 is similar):

  | field | Monte Carlo alone | bootstrap alone | co-added |
  |---|---|---|---|
  | ⟨v_R⟩ | 0.81 | 1.03 | 1.85 |
  | ⟨v_φ⟩ | 0.88 | 1.04 | 1.96 |
  | σ_R² | 1.00 | 1.02 | 2.05 |
  | σ_φ² | 1.23 | 1.00 | 2.22 |

  - χ²/N of the mock against the population: 0.53 with the co-added covariance, against 1.06 with the true one and
    1.28 with the bootstrap alone.
  - At 7.5% the distance re-draws (median shift 0.41 kpc) move most stars between 0.3 kpc rings. So the Monte Carlo on
    its own already re-samples each ring and reproduces 81–123% of the total variance.
  - A Monte Carlo is the right tool only for errors the bootstrap cannot see: errors correlated between stars, e.g. the
    zero point of the period–luminosity relation (one shift for all stars per draw).

## 2026-10-28 — Test 2b: truth χ²/N at 7.5% with the co-added Fourier error (user's recipe)

- Script: `scripts/test2b_coadded_errors_truth_chi2.py`; numbers in `results/test2b/coadded_errors_truth_chi2.json`.
- Error (user's recipe for a real observation): C_MC + C_boot.
  - C_MC: 300 re-draws of every observed star's distance (DM_obs + σ_DM N(0, 1)) propagated to the Fourier
    coefficients.
  - C_boot: the 300-bootstrap covariance of the unsmeared data, as in Test 1.
  - Hartlap factor with 300 samples.
- Proposed scheme, true potential, library from the observed stars; 817 terms (blocks usable in both data sets).

| error | library observed again | library as is |
|---|---|---|
| **C_MC + C_boot(unsmeared data)** | **2.18** | **2.20** |
| C_MC + C_boot(observed data) | 2.09 | 2.05 |
| C_boot(observed data) only | 5.40 | 5.07 |

- The perfect-data reference is 1.50 (bootstrap). The mock check above found this co-added covariance about 2× the
  true scatter in variance.
- The co-added error needs computing only once. At another trial angle, each mode's (cos, sin) pair rotates by
  m Δφ_b, and rings, layers, v_R and v_φ are unchanged. So c' = T c and C' = T C Tᵀ exactly; checked against the data
  fits saved at each angle (max difference 1.3×10⁻¹⁰). Rescoring: `scripts/test2b_qmc_7p5_rescore_coadded.py` →
  `results/test2b_qmc_7p5/evals_coadded.csv`.
- **First 10 points (truth + 9 old Sobol) with the co-added error:** one point beats the truth.
  - #5 (Ω 29.48, φ 25.8°, S 1.06, f 1.28): χ²/N 2.105 vs 2.176 for the truth, ΔlogL **+0.71** (observed-again
    library); +0.07 with the library as is.
  - On perfect data, #5 was ΔlogL −5.3.
  - The other 8 points are worse than the truth: −3.8 to −26.5.

## 2026-10-29 — Test 2b: the current method on a 7.5% distance-error mock (full focused Sobol grid)

- Scripts: `scripts/test2b_qmc_distance_error_run.py`, `scripts/test2b_qmc_7p5_rescore_coadded.py`.
  - Tables: `results/test2b_qmc_7p5/evals.csv` (bootstrap error) and `evals_coadded.csv` (co-added error).
  - Emulator runs `test2b_qmc_7p5` and `test2b_qmc_7p5_coadded`.
- Set-up:
  - Mock `mock/nbody_mock_distance_error_sigmaDM0.163_R8_N50k.hdf5`: σ_DM 0.163 mag, i.e. 7.5% (std of
    d_obs/d_true 7.54%); median position shift 0.41 kpc.
  - At all 257 points of the focused grid, the observed stars were integrated (9.8 h, finished 07:00).
  - Library: observed again with the same error (user's choice), and as is.
  - Proposed per-star Fourier scheme.
  - Error: co-added (user's recipe: distance Monte Carlo + bootstrap of the unsmeared data), or the bootstrap of the
    observed data.
- Truth: χ²/N 2.18 (co-added, observed again). Points better than the truth: 17 of 256 (observed again) and 22 (as
  is). Best: #248 (Ω 29.1, φ 24.3°, S 0.92, f 1.55) +2.24; #120 (29.0, 19.8°, 0.74, 1.96) +2.03; #180 (27.8, 20.0°,
  0.79, 1.53) +1.88.

| posterior (GP emulator, 256 Sobol) | \|Ω_b\| | φ_b [°] | S | bar_frac | LOO |
|---|---|---|---|---|---|
| perfect data (reference) | 26.06 +0.86/−0.78 | 27.3 ± 3.0 | 1.06 +0.25/−0.17 | 0.87 +0.31/−0.25 | 0.45 |
| **7.5%, co-added error, observed again** | **28.54 +1.95/−1.66** (−1.4σ) | **23.3 +4.7/−4.4** (+0.8σ) | **0.90 +0.23/−0.14** (+0.6σ) | **1.41 +0.38/−0.44** (−1.0σ) | 0.53 |
| 7.5%, co-added error, as is | 27.57 +1.89/−1.40 (−0.9σ) | 22.5 +5.1/−4.3 (+0.9σ) | 0.91 +0.22/−0.15 (+0.5σ) | 1.43 +0.39/−0.46 (−1.0σ) | 0.61 |
| 7.5%, bootstrap error, observed again | 28.52 +1.81/−1.34 (−1.6σ) | 22.7 +3.5/−3.3 (+1.3σ) | 0.84 +0.13/−0.09 (+1.4σ) | 1.63 +0.25/−0.36 (−2.1σ) | 1.31 |
| 7.5%, bootstrap error, as is | 27.39 +1.24/−0.84 (−1.3σ) | 21.4 +3.5/−3.4 (+1.6σ) | 0.87 +0.13/−0.10 (+1.1σ) | 1.55 +0.31/−0.37 (−1.6σ) | 1.55 |

(σ offsets are truth − median in units of the posterior half-width.)
- **All four shift the same way: faster, shorter and stronger bar at a smaller angle.** With the co-added error:
  Ω_b +2.5 km/s/kpc (+10%), φ_b −3.7°, S −10%, bar_frac +41%. The posterior is about 2× wider than on perfect data,
  and the truth lies within 1.4σ in every parameter.
- **The posterior is cut by the grid box.** The Ω_b marginal still rises at the upper edge (32), and bar_frac at 2.0.
  The box was chosen around the perfect-data posterior, so the medians and widths above are box-limited and the real
  shift may be larger.
- Figures: `figures/test2b_qmc_distance_error_7p5_coadded_error_emulator_corner_logLs_{observed,asis}.png`,
  `..._logl_corner_...`, and the bootstrap-error versions
  `figures/test2b_qmc_distance_error_7p5_emulator_corner_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024_{observed,asis}.png`.

## 2026-10-29 — Test 2: mock with distance errors and a DR3-Mira-like selection function

- Recipe (user):
  1. smear the distances of all N-body star particles first;
  2. select them so they follow the DR3 Miras in (l, b);
  3. and in their radial distribution.
- Script: `scripts/make_mock_selection.py`; mock `mock/nbody_mock_selection_sigmaDM0.163_N50k.hdf5`.
- Figure: `figures/test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.163.png`.
- **Steps.**
  - All 2 M star particles are observed with σ_DM = 0.163 mag (7.5%; l, b, proper motions and v_los exact;
    seed 2031). Particles with observed R < 8 kpc enter.
  - Target: DR3 Miras with DM_err < 0.4 and catalogue R < 8 kpc (45,206).
  - Selection function on observed quantities: S(l, b, d) = S_sky[HEALPix nside 32] × S_dist[0.25 kpc bins in
    heliocentric d] × S_R[0.25 kpc bins in Galactocentric R].
  - Fitted by iterative proportional fitting (10 cycles) to the three target margins, scaled to 50k expected stars.
  - Each particle is kept with probability min(S, 1) (seed 2032). The S tables are stored in the file, so S can be
    evaluated anywhere.
- **Result: 50,281 stars; no particle needed S > 1.** Sky-map correlation with the Miras 0.997. Medians:
  - distance: 8.04 vs 8.02 kpc;
  - R: 3.33 vs 3.33 kpc;
  - |b|: 3.41° vs 3.39°.
  - With sky and distance only, R was too concentrated (median 3.08), so R was added as a third margin.
  - True |z| > 1 kpc: 14.6%; disc particles: 89%.
- **Caveat.** The Miras' extinction gap at |b| < 1° is narrower than a HEALPix nside-32 pixel (1.8°), so the mock
  shows only a weak stripe there.
- **Not done yet: the library side.** The library must be weighted by the selection: 1/⟨S⟩ along each orbit and
  S per sample. Because S acts on observed quantities, S has to be evaluated on the library observed with the
  distance errors (or averaged over them).
- Spatial comparison: `figures/test2_spatial_mock_vs_dr3_miras_nbody_mock_selection_sigmaDM0.163_N50k.png`
  (`scripts/plot_mock_vs_dr3_mira_spatial.py`).
  - Fraction at |b| < 1°: Miras 5.1%, mock 9.1% (the extinction lane is under-resolved).
  - Fraction at |z| > 3 kpc: Miras 0.45%, mock 0.06% (absent from the N-body).
  - **Accepted by the user as good enough (2026-10-29).**

## 2026-10-29 — Pattern speed of the N-body run over 2-8 Gyr (snapshot choice for Test 2)

- Script: `scripts/measure_nbody_pattern_speed_evolution.py`; table `results/nbody_pattern_speed_evolution.csv`;
  figure `figures/nbody_pattern_speed_evolution.png`.
- Purpose: pick the snapshots for the planned faithfulness tests (10% distance errors + selection function) at
  several pattern speeds. The whole run `mock/Nbody/t_t0_*` is one secularly evolving galaxy, so different times
  give different Ω_b with the same physics.
- Method, reusing the existing mock machinery: `nbody.read_nemo` / `components` read each NEMO snapshot,
  `prepare_au18_mwframe.shrinking_sphere` centres it, the disc is spun to L_z < 0 (the mocks' MW-like frame), and
  Ω_b comes from the Dehnen, Semczuk & Schönrich (2023) Fourier method exactly as in
  `measure_au18_pattern_speed.py`. The bar is **not** rotated to −27° here, so the m = 2 phase can be compared
  between snapshots. 49 snapshots, ~4 s each.
  - Ω(R) and hence R_CR come from the snapshot's **own accelerations** (azimuthal average of a_R in a |z| < 0.25 kpc
    slab), not a spherical enclosed-mass estimate: at t = 7 this gives R_CR = 8.61 kpc against the 8.61 kpc of the
    AGAMA ground-truth potential, while the spherical estimate is 2 kpc too small.
- **Validation:** at t = 7 the script returns Ω_p = 26.0170 ± 0.1933, R_bar = 5.00 kpc, R_CR = 8.608 kpc, against
  `results/nbody_ground_truth.json`'s 26.0168 ± 0.1926, 5.00 and 8.615 — the same numbers to 4 significant figures.

| t (model) | t [Gyr] | Ω_b [km/s/kpc] | A₂^max | R_bar [kpc] | R_CR [kpc] | R_CR/R_bar |
|---|---|---|---|---|---|---|
| 2.125 | 2.08 | 41.92 | 0.430 | 3.75 | 5.31 | 1.42 |
| 3.125 | 3.06 | 38.23 | 0.455 | 4.00 | 5.84 | 1.46 |
| 4.125 | 4.03 | 33.58 | 0.468 | 4.25 | 6.59 | 1.55 |
| 5.125 | 5.01 | 29.87 | 0.445 | 4.50 | 7.43 | 1.65 |
| 6.125 | 5.99 | 28.25 | 0.423 | 4.75 | 7.86 | 1.65 |
| **7** | **6.84** | **26.02** | **0.397** | **5.00** | **8.61** | **1.72** |
| 8.125 | 7.94 | 24.21 | 0.357 | 5.25 | 9.32 | 1.78 |

- **The bar slows steadily:** 41.9 → 24.2 km/s/kpc over 2.1–7.9 Gyr, a 42% drop, roughly linear at
  −3.1 km/s/kpc per Gyr. It is smooth, so any Ω_b in this range can be had by picking a snapshot.
- As it slows it grows and weakens: R_bar 3.75 → 5.25 kpc, A₂^max peaks at 0.47 near 4 Gyr then falls to 0.36.
  R_CR/R_bar rises 1.42 → 1.78, i.e. the bar goes from fast to slow in the usual classification.
  v_c(8 kpc) stays within 220–227 km/s, so the mass model is near enough unchanged across the range.
- **Consistency check** (`phase_check`): the snapshots are 0.122 Gyr apart, over which the bar turns ~1.25π, so the
  m = 2 phase (defined mod π) is aliased and cannot be differenced to get Ω_b — hence the single-snapshot method.
  Used in reverse, per step, it still tests Ω_b: the median residual is −1.5 ± 6.5°, i.e. a systematic Ω_b offset of
  −0.21 ± 0.91 km/s/kpc, consistent with zero. The per-step scatter (rms 36°) is the bar's own phase wobble, not a
  method error. This must be done per step: integrating across 2–8 Gyr accumulates 60π of phase, where a 1% error in
  Ω_b already wraps.
- **Caveat for the tests ahead.** Ω_b changes by ~0.4 km/s/kpc over the 1 Gyr the steady-state method integrates
  (~1.5% at t = 7). The method assumes a *steady* bar, so this is a real, if small, model mismatch, built into the
  test galaxy at every snapshot and largest at early times where the bar slows fastest.

## 2026-10-29 — Test 2 mocks: 10% distance errors, with and without the selection function, at five pattern speeds

- Scripts: `scripts/prepare_nbody_mwframe.py --time T`, `scripts/make_mock_perfect.py`,
  `scripts/make_mock_distance_error.py --frac 0.10`, `scripts/make_mock_selection.py --galaxy nbody_t<T> --frac 0.10`.
- Purpose: test the faithfulness of the method under a 10% distance error **and** a DR3-Mira-like selection
  function, across the pattern speeds the run provides (previous entry).
- **The snapshot times below are model times** (the `t_t0_<T>` file names), not Gyr: 1 model unit = 0.9778 Gyr,
  so model t = 7 is 6.84 Gyr. Both columns are given.
- **Error level:** 10% in distance, i.e. σ_DM = 0.2171 mag (σ_DM = 0.10 · 5 / ln 10), up from the 7.5%
  (0.1629 mag) of the earlier runs. `make_mock_distance_error.py` and `make_mock_selection.py` now take
  `--frac` directly.
- **Galaxy registry:** each snapshot is its own test galaxy, `nbody_t<T>` (`galaxies.nbody_name(t)`), with all its
  files tagged `nbody_t<T>_...`. `nbody` remains model t = 7 on the original un-tagged names, so every earlier
  result stays valid.

| t (model) | t [Gyr] | Ω_b [km/s/kpc] | A₂^max | R_bar | R_CR | set | N | median R |
|---|---|---|---|---|---|---|---|---|
| 2 | 1.96 | 44.25 | 0.472 | 4.25 | 5.07 | selection | 50,130 | 3.33 |
| 3 | 2.93 | 39.53 | 0.446 | 4.00 | 5.65 | selection | 50,079 | 3.31 |
| 4 | 3.91 | 33.96 | 0.472 | 4.25 | 6.49 | selection | 50,239 | 3.34 |
| 5.5 | 5.38 | 28.84 | 0.448 | 4.50 | 7.72 | selection | 50,077 | 3.32 |
| **7** | **6.84** | **26.02** | **0.397** | **5.00** | **8.61** | **selection** | **50,281** | **3.33** |
| **7** | **6.84** | **26.02** | **0.397** | **5.00** | **8.61** | **benchmark, no selection** | **50,000** | **1.88** |

- Files: `mock/<prefix>_mock_selection_sigmaDM0.217_N50k.hdf5` and, for the benchmark,
  `mock/nbody_mock_distance_error_sigmaDM0.217_R8_N50k.hdf5`
  (prefix = `nbody` for t = 7, `nbody_t<T>` otherwise). Figures
  `figures/test2_<prefix>_mock_selection_vs_dr3_miras_sigmaDM0.217.png`.
- **The selection matters a lot radially:** median R is 3.33 kpc in the selection mocks against **1.88 kpc** in the
  no-selection benchmark. The benchmark draws uniformly from the centrally concentrated R < 8 kpc particles, while
  the Miras' sky + distance selection pushes the sample outwards. This is the main thing the benchmark isolates.
- All five selection mocks reproduce the DR3-Mira target margins: sky-map correlation 0.997–0.998, median distance
  8.02–8.05 vs 8.02 kpc, median R 3.31–3.34 vs 3.33 kpc, median |b| 3.38–3.40 vs 3.39°. Fewer than 0.01% of
  particles needed S > 1.
- **Validation of the refactor.** Re-running the 7.5% selection mock on the new code reproduces the 2026-10-29 entry
  exactly (50,281 stars, sky correlation 0.997, 0.00% capped), and the rebuilt t = 7 products reproduce the recorded
  numbers: pool 1,493,610 particles at R < 8 kpc, mock split 40,023 disc / 9,977 bulge, A₂^max 0.397, R_bar 5.00 kpc.
- **Not done yet:** the potentials for the four new snapshots (`build_nbody_potential.py`, `potential_family.py`) and
  their `*_ground_truth.json`, all of which the orbit-integration stage needs. Ω_b, R_bar and R_CR above come from
  `results/nbody_pattern_speed_evolution.csv` and can serve as the truth values.
- **Still open (unchanged):** the library side of the selection function — orbits weighted by 1/⟨S⟩ and samples by
  S(x_j), with S evaluated on the library *observed* with the distance errors (section 5.2 of PIPELINE.md).

## 2026-10-29 — Test 2 mocks: spatial distribution and the v_R quadrupole

- Script: `scripts/plot_test2_mock_spatial_and_vR.py` (`--frac` picks the error level; `--mocks` any file list).
- Figures: `figures/test2_mock_spatial_10pct.png`, `figures/test2_mock_vR_quadrupole_10pct.png`.
- **Spatial** (x–y, x–z, l–b; 0.2 kpc and 1° pixels), one column per set with the DR3 Miras as reference. The five
  selection mocks track the Miras closely in all three projections. The no-selection benchmark stands out: it is
  visibly more centrally concentrated face-on and much more peaked at l ≈ 0 on the sky, which is the 1.88 vs
  3.33 kpc difference in median R seen before. The Miras' extinction lane at |b| < 1° is sharper than the mocks'
  (the known HEALPix nside-32 resolution limit).
- **v_R quadrupole.** Face-on ⟨v_R⟩ maps (|z| < 1 kpc) show the bar's quadrupole clearly in every set, and the m = 2
  amplitude and node are fitted per ring with `test1_nbody_qmc_focused_fourier_stars.fit` — the same per-star
  least-squares the likelihood uses, which fits ⟨v_R⟩(φ) rather than summing v_R cos 2φ and so is **not** biased by
  the selection function.

| set | N | N(\|z\|<1) | m=2 amp obs | true | change | node obs | node true | shift |
|---|---|---|---|---|---|---|---|---|
| t = 2 | 50,130 | 43,649 | 33.65 | 31.18 | +7.9% | −10.7° | −20.7° | +10.0° |
| t = 3 | 50,079 | 43,347 | 27.43 | 24.94 | +10.0% | −20.9° | −31.7° | +10.8° |
| t = 4 | 50,239 | 43,439 | 23.56 | 24.71 | −4.6% | −15.9° | −25.9° | +10.1° |
| t = 5.5 | 50,077 | 43,007 | 22.04 | 19.60 | +12.4% | −11.7° | −22.9° | +11.2° |
| t = 7 | 50,281 | 42,992 | 21.41 | 19.35 | +10.6% | −16.9° | −24.6° | +7.7° |
| t = 7, no selection | 50,000 | 42,177 | 22.12 | 18.75 | +18.0% | −18.6° | −23.9° | +5.3° |
| DR3 Miras (reference) | 45,206 | — | — | — | — | — | — | — |

(medians over 1 < R < 4 kpc, |z| < 1 kpc; "true" = the same stars at their true positions.)
- **The quadrupole amplitude falls with the bar's pattern speed,** 33.7 km/s at t = 2 (Ω_b 44.3) down to 21.4 at
  t = 7 (Ω_b 26.0), tracking A₂^max.
- **The 10% distance error rotates the measured node by +8 to +11° in every set,** i.e. the quadrupole appears at a
  *smaller* bar angle than the truth. The direction matches the posterior shift already seen at 7.5% (φ_b −3.7°,
  entry of 2026-10-29), so this is the same systematic seen directly in the data rather than through the posterior.
  The no-selection benchmark shifts least (+5.3°), so **the selection function roughly doubles the bias** — worth
  keeping in mind when the library-side weighting is written.
- The amplitude is *not* damped by the distance error; it rises by ~10% (+18% without selection). Smearing along the
  line of sight moves stars in R and φ, which mixes the strong inner quadrupole outwards.
- **Validation of the estimator:** on all 1,493,610 parent particles at their true positions the fit gives a node of
  −26.2° against the bar's −27.0° (the repo's density-weighted `quadrupole_node` gives −25.5°), so the shifts above
  are physical, not an artefact of the fit. On the smeared 50k benchmark the density-weighted estimator moves far
  more (−26.1° → −10.5°) than the least-squares one (−23.9° → −18.6°), another reason to prefer the latter.

## 2026-10-29 — Pipeline rebuilt and run on this machine: ground-truth potential and chi2/N at the truth

- New script: `scripts/run_single_point.py` — the per-point work of the global searches for a single parameter
  point, so one data set can be scored without launching a Sobol design. Integrates, fits the data and library
  Fourier modes, and reports chi2/N for every scheme in `fs.SCHEMES`. Resumable (re-uses a saved trajectory
  unless `--fresh`).
  - Outputs: `$MIRA_DR4_RUNS/single_point/{traj_<tag>.npy, point_<tag>.npz}` and
    `results/single_point/single_point.csv`.
  - Two library variants when the mock has a distance error, as in `test2b_qmc_distance_error_run.py`:
    `asis` (model not smeared) and `observed` (library re-observed with the same error, seed 11).
- **Ground-truth potential rebuilt here** (`build_nbody_potential.py`, `potential_family.py --galaxy nbody`,
  `measure_au18_pattern_speed.py --galaxy nbody`). The disc CylSpline is the expensive part: 172 s of the build.
  - v_c(8 kpc) = 221.1 km/s against the 222 recorded; Omega_p = 26.02 ± 0.19; R_CR = 8.60 kpc (recorded 8.615);
    R_CR/R_bar = 1.72. Resonances: ILR 3.40, CR 8.60, OLR 14.74 kpc.

| data (t = 7, no selection) | N | chi2/N, library as is | chi2/N, library observed | integrate [s] | total [s] |
|---|---|---|---|---|---|
| perfect mock | 50,000 | **1.439** | — | 162.9 | 169.8 |
| 7.5% distance error | 50,000 | **4.612** | 5.265 | 149.0 | 157.0 |
| 10% distance error | 50,000 | **8.411** | 8.910 | 148.2 | 156.1 |

(truth point, |Omega_b| = 26.02, phi_b = 27°, S = bar_frac = 1; scheme
`fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024`, N_data 855 (836 at 7.5%), accuracy 1e-5, 300 bootstraps.)
- **Cost on this machine** (i7-10700KF, 8 cores / 16 threads): ~150 s for 50k orbits x 1 Gyr x 41 snapshots, plus
  4.7 s for the data fit and 2 s for the library fit and chi2 — about **2.6 min per design point**, so a 257-point
  search is roughly **11 h**. The reference Apple M5 Pro took 80–95 s per integration.
- **Thread scaling measured** (10k orbits): 4 threads 178 orbits/s, 8 → 356, 12 → 391, 16 → 426. Hyperthreading is
  worth ~20% over the physical core count, so `common.py` now defaults `OMP_NUM_THREADS` to `os.cpu_count()`
  (an explicit setting or a scheduler's choice still wins). chi2 is bit-identical at 8 and 16 threads.
- **The extra misfit from the distance error** is +2712 chi2 at 7.5% and +5960 at 10% over the perfect-data value,
  a ratio of 2.20 for a 1.33x larger error, i.e. a power of 2.7. That sits between the 3.5 the earlier scan found
  at 2.5–5% and the 2.4 implied by its 5% → 12.2% points, so the steep growth flattens with error, as expected.
- **Smearing the library the same way makes the fit slightly worse, not better** (4.61 → 5.27 at 7.5%,
  8.41 → 8.91 at 10%), for the no-selection data. The 7.5% search saw the opposite ordering on the posterior, so
  this is worth watching when the selection-function runs are scored.
- **Caveat: chi2/N on perfect data is 1.439 here against the 1.498 recorded** (same N_data = 855, same mock — the
  mock reproduces exactly: pool 1,493,610, split 40,023 disc / 9,977 bulge). The likely cause is the AGAMA build:
  1.0 (May 2024) here against 1.0.156 for the recorded runs, which changes the potential expansion coefficients and
  so the orbits slightly. PIPELINE.md section 7 already flagged that bit-identical integration across machines was
  untested; this is the first measurement of the size of the difference, about 4% in chi2/N at the truth.

## 2026-10-29 — Test 2 posterior runs launched: six searches at 10% distance error (in progress)

- Scripts: `scripts/test2b_qmc_distance_error_run.py` (extended), `scripts/selection.py` (new),
  `scripts/run_test2_posteriors.sh`, `scripts/run_test2_posteriors_cluster.sh`.
- Six Sobol searches of 257 points each, then the GP emulator and corner plots, over the 10% mocks:

| run | data | machine |
|---|---|---|
| `test2_t7_nosf_10pct` | t = 7, no selection (benchmark) | workstation |
| `test2_t7_sf_10pct` | t = 7, selection function | workstation |
| `test2_t2_sf_10pct` | t = 2, selection function | workstation |
| `test2_t4_sf_10pct` | t = 4, selection function | workstation |
| `test2_t5.5_sf_10pct` | t = 5.5, selection function | 88-core cluster |
| `test2_t3_sf_10pct` | t = 3, selection function | 88-core cluster |

- **What was extended, following PIPELINE sections 7.7-7.8 rather than writing a new runner:**
  - mock, run directory and galaxy are arguments; `--stride/--offset` give each worker a subset of the design;
    `--merge` builds evals.csv from the per-point files.
  - Points are claimed with an atomic `O_EXCL` lock, so the two machines share one output directory without
    duplicating work; a lock older than 3 h with no result is reclaimed, so a dead worker strands nothing.
  - The data fit and its 300 bootstraps are done **once** at the true bar angle and rotated to every trial angle
    (`fs.rotation`, moved there from the rescore script, which now imports it). Checked on each data set at
    start-up: max |direct - rotated| about 1e-10. Saves a refit at every point.
  - The design's Omega axis is scaled by the galaxy's own pattern speed (`scaled_design`), so the same 257 Sobol
    points cover 20.1-41.8 for t = 5.5 instead of the fiducial 20-32. The other three axes are unchanged.
  - `build_nbody_potential.py` takes `--galaxy`; `potential_family.py` now resolves any `nbody_t<T>` name.
  - The shell script builds a galaxy's trial potential if it is missing, before launching any worker.
- **Library weighting for the selection function (`selection.py`), the piece listed as "not done yet":**
  weight = S(x_j) / tau_i with tau_i = <S> along the orbit, S evaluated on the library *observed* with the mock's
  distance error (bar frame -> MW frame -> l, b, d -> smear -> back). With S replaced by the R < R_SEL indicator
  this is exactly the Test 1 weighting, so the two are one scheme with different S. All three variants
  (`asis`, `observed`, `sf`) are scored from the same trajectory, so the extra ones cost no integration.
  First check at the t = 7 truth: chi2/N 7.82 (`sf`) against 8.17 (`asis`) and 10.23 (`observed`), and a wrong
  design point gives 13.7-17.1, so the scheme discriminates.
- **Cost:** about 150 s per point for the no-selection mock and 110 s for the selection mocks (sparser centre),
  so roughly 10.5 h and 8 h per search. Trajectories are kept (47 MB each, about 12 GB per run); the runner stops
  saving them below 25 GB free rather than filling the disc.
- **Checked before committing two days of compute:** the GP emulator, its corner/slices/LOO figures and
  `plot_test1_logl_corner.py` were all run against the recorded 257-point 7.5% search renamed to the new
  convention. They reproduce the logged posterior (Omega_b 28.52 +1.73/-1.38 against 28.52 +1.81/-1.34;
  bar_frac 1.601 against 1.63), so the whole chain works end to end.
- **Caveat to read the results with: bar_frac will be box-limited.** The focused design caps bar_frac at 2.0 and
  the 7.5% posterior already piles up against that edge; at 10% the shift is larger. Omega was rescaled per
  galaxy but the other axes were left as they are, so the bar_frac medians and widths are lower bounds on the
  shift, not measurements. Widening that axis is a clean follow-up on the same machinery.

## 2026-10-29 — Test 2 results: posteriors at 10% distance error (RUNNING; two of six done)

<!-- TEST2-RESULTS-TABLE -->
Posteriors from the GP emulator of `logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024_<variant>`; `sf` is the
selection-function library weighting, which is the correct one for a selection mock. Offsets are
(truth − median) in units of the posterior half-width.

| run | Ω_b | φ_b [°] | S | bar_frac | LOO |
|---|---|---|---|---|---|
| **t = 3, SF** (truth Ω 39.53) | 42.33 +1.60/−1.54 (−1.8σ) | 19.80 +4.57/−3.24 (+1.8σ) | 0.784 +0.221/−0.080 (+1.4σ) | 1.375 +0.203/−0.497 (−1.1σ) | 2.15 |
| **t = 5.5, SF** (truth Ω 28.84) | 33.66 +1.13/−1.82 (−3.3σ) | 24.87 +4.08/−5.50 (+0.4σ) | 0.949 +0.158/−0.073 (+0.4σ) | 1.405 +0.233/−0.478 (−1.1σ) | 2.84 |
| **t = 7, no SF** `asis` (truth Ω 26.02) | 27.97 +1.18/−0.98 (−1.8σ) | 18.32 (+2.7σ) | 0.708 (+4.5σ) | 1.730 (−3.2σ) | 3.57 |
| **t = 7, no SF** `observed` | 29.49 +1.67/−1.33 (−2.3σ) | 22.53 (+0.8σ) | 0.728 (+2.2σ) | 1.868 (−5.6σ) | 4.09 |
| **t = 7, SF** `sf` | 28.23 (−2.2σ) | 21.57 (+1.2σ) | 0.823 (+2.3σ) | 1.367 (−1.4σ) | 3.00 |
| **t = 2, SF** `sf` (truth Ω 44.25) | 43.73 (−1%) | 17.0 | 0.986 | 0.846 | 2.93 |
| **t = 4, SF** `sf` (truth Ω 33.96) | 38.48 (+13%) | 24.2 | 0.913 | 1.359 | 2.33 |

All six searches finished 2026-10-31 01:47: 257/257 design points each, 1,542 orbit integrations.
The complete per-variant table is in the entry below.

Truth χ²/N at each point, by library variant:

| run | sf | asis | observed | points beating the truth (sf) |
|---|---|---|---|---|
| t = 3, SF | **6.512** | 7.496 | 8.217 | 14/256 |
| t = 5.5, SF | **7.886** | 8.271 | 9.595 | 18/256 |

- **The selection-function weighting works:** `sf` gives the lowest truth χ²/N on both runs, below both `asis`
  (which ignores the selection) and `observed`. That is the first evidence that the new library weighting in
  `selection.py` is doing what it should.
- **The bias has the same direction as the 7.5% no-selection search: faster, shorter, stronger bar at a smaller
  angle.** bar_frac is the most repeatable, +38% and +41% here against the +41% recorded at 7.5%.
- **Ω_b, the parameter the method exists to measure, is biased high and by a snapshot-dependent amount:** +7.1%
  at t = 3 and +16.7% at t = 5.5. Whether that scatter is the pattern speed itself, the bar strength (A₂ 0.446 vs
  0.448) or just realisation noise needs the remaining four runs.
- **At 10% with a selection function the truth is not recovered within the stated errors** (1.8σ and 3.3σ).
- **Emulator quality is much worse than on perfect data:** LOO rms 2.15 and 2.84 against 0.45, so the likelihood
  surface is noisy and the quoted widths should be treated as approximate.
- **Caveats.** The Ω_b posterior at t = 5.5 reaches the top of the emulator box (35.47); the peak itself is
  contained — the 20 design points that extend to Ω = 41.8 score ΔlogL −1.56 at best against +4.81 inside — but the
  upper tail is clipped, so the + error is a lower bound. bar_frac remains box-affected at 2.0 throughout.

## 2026-10-30 — The selection-function library weighting is NOT yet validated (read `sf` and `asis` together)

`selection.py`'s `sf` weighting lowers χ²/N at the truth on every selection mock, but that does **not** carry over
into better parameter recovery. Both variants, same data, same design:

| run | lib | Ω_b bias | S | bar_frac | worst offset | LOO |
|---|---|---|---|---|---|---|
| t = 3 | asis | +8% | 0.798 | 1.707 | **4.9σ** | 2.86 |
| t = 3 | **sf** | +7% | 0.784 | 1.375 | **1.8σ** | 2.15 |
| t = 5.5 | asis | **−2%** | 0.897 | 1.186 | 3.5σ | 2.12 |
| t = 5.5 | **sf** | **+17%** | 0.949 | 1.405 | 3.3σ | 2.84 |
| t = 7 | asis | +3% | 0.889 | 1.296 | **1.4σ** | 2.24 |
| t = 7 | **sf** | +9% | 0.823 | 1.367 | 2.3σ | 3.00 |

- **No consistent winner.** `sf` wins clearly at t = 3 (worst offset 1.8σ against 4.9σ, mostly through bar_frac),
  loses at t = 7, and at t = 5.5 turns a −2% Ω_b bias into +17%.
- **`sf` biases Ω_b high every time** (+7%, +17%, +9%) while `asis` scatters about zero (+8%, −2%, +3%). A one-sided
  bias in the one parameter the method exists to measure looks more like an implementation artefact than noise.
- **Variance inflation is not the cause** — measured on the t = 7 truth trajectory, the `sf` weights have a *higher*
  effective sample size than the `asis` ones (51.3% of 1.51 M samples against 31.2%) and a smaller max/median
  weight ratio (30.9 against 40.0). 10.6% of library samples fall outside the selection (S = 0).
- **Two untested hypotheses for the one-sided Ω_b bias:**
  1. τ_i = ⟨S⟩ is averaged over the t ≥ 0.25 Gyr snapshots only, while the `asis` τ uses all t > 0; the
     inconsistency sits in the 1/τ normalisation.
  2. The data were selected with one realisation of the distance errors, but the library re-draws its own
     (seed `OBS_SEED` = 11) when evaluating S. Correct in expectation, but not the same S field.
- **How to read the results:** quote `sf` and `asis` side by side for every selection mock, not `sf` alone. A lower
  χ²/N at the truth is not evidence that a weighting recovers parameters better, and it was wrong to present it
  that way.

## 2026-10-30 — It is the library *smearing*, not the selection weighting, that biases Ω_b high

Follow-up to the entry above, which wrongly blamed `selection.py`. Running the emulator on all three library
variants of the same searches separates the two effects: `asis` is neither smeared nor S-weighted, `observed` is
smeared only, `sf` is smeared **and** S-weighted.

| run | `asis` (no smear, no S) | `observed` (smear only) | `sf` (smear + S) |
|---|---|---|---|
| t = 3, SF | +8% | +12% | +7% |
| t = 5.5, SF | −2% | +13% | +17% |
| t = 7, no SF | +7% | +13% | — |
| t = 7, SF | +3% | +7% | +9% |

(Ω_b bias = posterior median / truth − 1.)

- **Smearing the library is the systematic.** `asis` scatters about +4% on average; `observed` is consistently
  +7 to +13%, in all four runs. Adding the selection weighting on top then moves Ω_b only modestly and in both
  directions (+12 → +7, +13 → +17, +7 → +9), so **`selection.py` is not the source of the one-sided bias** and the
  previous entry's suspicion of it was wrong.
- **The effect is already in the recorded 7.5% search**, same sign and similar size: Ω_b 28.54 / 28.52 observed
  against 27.57 / 27.39 as is, i.e. +3.7 and +4.3 percentage points. It was in the numbers but had not been
  identified as a systematic.
- **The τ definition is not involved.** `sf` averages S over the t ≥ 0.25 Gyr snapshots while `asis` uses all
  t > 0; measured on the t = 7 truth trajectory the two τ agree to a median ratio of 0.9999 (5–95% 0.85–1.23),
  unbiased scatter that cannot produce a one-sided bias. Hypothesis 1 of the previous entry is therefore dead.
- **The open question for the user.** Observing the library with the same errors as the data is the principled
  recipe and the one chosen for the 7.5% run, yet it recovers Ω_b *worse* than leaving the library unsmeared, in
  4 of 4 runs at 10% and in both error models at 7.5%. Either the smearing is being applied in a way that biases
  the moments (one error draw per point, seed `OBS_SEED` = 11, rather than an average over draws, is the obvious
  suspect), or the unsmeared library is closer to the truth for a reason that does not generalise. This should be
  settled before the recipe is used on real data.

## 2026-10-31 — Test 2 complete: all six 10% posteriors, and which library weighting to use

All six searches done (257 points each, 1,542 integrations). Every run has `evals.csv`, posteriors and corner
plots under `results/<run>/` and `figures/<run>_emulator_corner_*.png`.

| run | lib | Ω_b | bias | φ_b | S | bar_frac | worst | LOO |
|---|---|---|---|---|---|---|---|---|
| t = 2, SF (truth 44.25) | asis | 43.03 | −3% | 20.2 | 0.937 | 1.311 | 2.1σ | 2.60 |
| t = 2, SF | sf | 43.73 | −1% | 17.0 | 0.986 | 0.846 | 5.3σ | 2.93 |
| t = 3, SF (truth 39.53) | asis | 42.51 | +8% | 23.9 | 0.798 | 1.707 | 4.9σ | 2.86 |
| t = 3, SF | observed | 44.11 | +12% | 20.9 | 0.830 | 1.338 | 2.1σ | 2.32 |
| t = 3, SF | sf | 42.33 | +7% | 19.8 | 0.784 | 1.375 | 1.8σ | 2.15 |
| t = 4, SF (truth 33.96) | asis | 36.37 | +7% | 22.8 | 0.895 | 1.489 | 2.2σ | 2.48 |
| t = 4, SF | sf | 38.48 | +13% | 24.2 | 0.913 | 1.359 | 2.7σ | 2.33 |
| t = 5.5, SF (truth 28.84) | asis | 28.19 | −2% | 17.5 | 0.897 | 1.186 | 3.5σ | 2.12 |
| t = 5.5, SF | observed | 32.62 | +13% | 17.9 | 0.614 | 1.959 | 14.2σ | 2.33 |
| t = 5.5, SF | sf | 33.66 | +17% | 24.9 | 0.949 | 1.405 | 3.3σ | 2.84 |
| t = 7, SF (truth 26.02) | asis | 26.85 | +3% | 20.9 | 0.889 | 1.296 | 1.4σ | 2.24 |
| t = 7, SF | observed | 27.79 | +7% | 17.8 | 0.807 | 1.161 | 2.9σ | 2.68 |
| t = 7, SF | sf | 28.23 | +9% | 21.6 | 0.823 | 1.367 | 2.3σ | 3.00 |
| t = 7, no SF | asis | 27.97 | +7% | 18.3 | 0.708 | 1.730 | 4.5σ | 3.57 |
| t = 7, no SF | observed | 29.49 | +13% | 22.5 | 0.728 | 1.868 | 5.6σ | 4.09 |

| library variant | mean \|Ω_b bias\| | mean worst offset | mean LOO |
|---|---|---|---|
| **asis** (not smeared, no S) | **5.0%** | **3.1σ** | **2.65** |
| observed (smeared, no S) | 11.2% | 6.2σ | 2.85 |
| sf (smeared + S) | 9.4% | 3.1σ | 2.65 |

- **At 10% distance error the method does not recover the truth within its errors.** The best variant, `asis`,
  still misses Ω_b by 5% on average, and the worst parameter is 3σ out in the mean. Ω_b is never recovered to
  better than 1%.
- **The selection function is not the main problem; the distance error is.** At t = 7, where both exist, adding
  the selection function changes Ω_b from +7% (no SF, `asis`) to +3% (SF, `asis`) — no worse, and bar_frac and S
  improve. The 10% smearing dominates.
- **Do not smear the library.** `observed` is the worst variant on every measure, biasing Ω_b to +7…+13% in all
  four runs where it was computed (see the entry below for the isolation of this effect). The `sf` weighting
  inherits the same smearing and most of the same bias.
- **The `sf` weighting is still not validated.** It matches `asis` on worst-offset and LOO but is worse on Ω_b
  (9.4% against 5.0%), and the sign of its advantage flips between snapshots: better at t = 3 (1.8σ against
  4.9σ), much worse at t = 2 (5.3σ against 2.1σ). Quote both variants.
- **No trend of the bias with pattern speed.** `asis` gives −3%, +8%, +7%, −2%, +3% going from Ω_b 44 down to 26,
  i.e. scatter, not a trend. Whatever drives the bias is not the bar's speed.
- **Caveat that limits several rows: bar_frac is box-limited at 2.0.** 1.868, 1.959 and 1.730 sit on the edge with
  no design points beyond it, so those are truncated, not measured; the t = 5.5 `observed` row (14.2σ, S = 0.614)
  is pathological for this reason. A re-run with a wider bar_frac axis is the obvious next step.
- **Emulator quality is the other limit:** LOO rms 2.1–4.1 against 0.45 on perfect data, so all quoted widths,
  and hence the σ offsets, are approximate.

## 2026-10-31 — Why Ω_b came out high: the dispersion ring averages. Scheme B fixes it, and the principled recipes then win

This supersedes the two entries below it, which concluded that smearing the library is harmful and that
`selection.py` does not help. Both conclusions were artefacts of the χ² scheme.

**The mechanism**, measured on the saved coefficients (no new integrations):

1. *Every library is hotter than the data.* Ring-averaged σ_R over 1 < R < 6 kpc, |z| < 1, at the truth point:

   | | σ_R − data |
   |---|---|
   | data (t = 7 SF) | 97.99 km/s |
   | library `asis` | +7.3 |
   | library `sf` | +8.1 |
   | library `observed` | **+11.6** |

   The orbits are seeded from distance-smeared stars, which inflates their energies, so even the unsmeared
   library runs ~7 km/s hot — the "hot library" problem of the perfect-data runs (~1.5 km/s), five times worse at
   10%. Re-observing the library adds another ~4.5 km/s.
2. *Ω_b is the dominant lever on library temperature.* Regressing library σ_R on the design (257 points,
   standardised), in km/s per 1 sd of each parameter:

   | variant | Ω_b | φ_b | ln S | ln bar_frac |
   |---|---|---|---|---|
   | asis | **−3.84** | +0.29 | +0.73 | +0.60 |
   | observed | **−3.46** | +0.51 | +0.45 | +0.36 |
   | sf | **−3.00** | +0.46 | +0.37 | +0.24 |

   Ω_b is about ten times more effective at cooling the library than S or bar_frac.
3. *Scheme A fits the σ_R², σ_φ² m = 0 terms*, i.e. the ring-averaged dispersions, so the fit buys agreement in
   library temperature by raising Ω_b. The hotter the library, the further it goes — which is exactly the ordering
   seen in the posteriors.

**Scheme B** (`fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024`) drops those four coefficients per block and
keeps everything else. It was introduced for the perfect-data hot library; at 10% it is decisive.

| | scheme A | scheme B |
|---|---|---|
| mean \|Ω_b bias\| | 8.1% | **4.2%** |
| mean worst offset | 3.9σ | **1.6σ** |
| mean LOO rms | 2.70 | **1.20** |

**And under scheme B the ordering of the library variants reverses** — the principled recipes win:

| | mean \|Ω_b bias\|, scheme B |
|---|---|
| `sf` (smeared + selection weighting) | **2.1%** |
| `observed` (smeared) | 2.7% |
| `asis` (not smeared, selection ignored) | 6.2% |

`sf` beats `asis` in 5 of 5 selection mocks. Full scheme-B results (truth φ_b 27, S 1, bar_frac 1):

| run | lib | Ω_b | bias | φ_b | S | bar_frac | worst |
|---|---|---|---|---|---|---|---|
| t = 2 | sf | 43.40 | −1.9% | 18.8 | 1.009 | 0.817 | 2.4σ (φ_b) |
| t = 3 | sf | 40.63 | +2.8% | 24.5 | 0.953 | 0.939 | 0.6σ |
| t = 4 | sf | 34.86 | +2.7% | 23.8 | 1.104 | 0.853 | 0.6σ |
| t = 5.5 | sf | 28.06 | −2.7% | 21.3 | 0.991 | 0.897 | 0.9σ |
| t = 7, SF | sf | 25.92 | −0.4% | 22.9 | 0.989 | 0.947 | 0.5σ |
| t = 7, no SF | observed | 26.73 | +2.7% | 27.0 | 1.005 | 1.069 | 0.5σ |

- **S and bar_frac are recovered to about 1σ** (0.95–1.10 and 0.82–0.95), against 0.71–0.95 and 1.2–1.9 under
  scheme A where bar_frac railed at the box edge. The box limit is no longer reached.
- **The t = 7 no-selection run with the correct forward model is essentially perfect** in all four parameters.
- **φ_b is the one surviving bias:** 18.8–24.5 against 27, low by ~4° and worst at t = 2 (2.4σ). This matches the
  quadrupole-node shift measured directly in the mocks (+8 to +11° at 10%), so it is the same effect seen through
  the posterior rather than a separate problem.
- **Recommendation:** use scheme B with the matched forward model (`observed` without a selection function, `sf`
  with one). Scheme A should not be used on data with distance errors at this level.
- Figures: `figures/test2_pattern_speed_recovery_scheme{A,B}.png`
  (`scripts/plot_test2_pattern_speed_recovery.py --scheme ...`).

## 2026-11-01 — Two checks: the pattern speed is 26.1 (not 25), and the Ω_b error bars are too small by ~1.5x

**1. Is the Dehnen Ω_b right?** The simulation directory has its own table,
`SchwarMAX/Bar_model_TG21/model/bar_strength_and_pattern_speed.csv`, which gives 24.64 (its Dehnen column) and
25.20 (its finite-difference column) at t = 7, against our 26.02. Checked with a third, assumption-free method
(`scripts/check_nbody_pattern_speed_independent.py`): track the m = 2 phase between snapshots in a **fixed**
frame and difference it. No velocities, no method, just where the bar is.

| t | measured m = 2 phase |
|---|---|
| 6.750 | +58.07° |
| 6.875 | +49.33° |
| 7.000 | +41.69° |
| 7.125 | +35.56° |
| 7.250 | +30.32° |

Per step the bar turns 180° + (5.2 to 8.7)°, giving Ω_b = 26.35, 26.20, 25.99, 25.86 — a decelerating bar.
Fitting all five: **Ω_b = 26.10 ± 0.05, rms residual 1.0°**, identical (to 0.02) across annuli 0.5–3, 1–4, 1–5
and 2–5 kpc and for disc-only. Per-step residuals are 1.4° at Ω_b = 26.10, but **6.6° at 25.20 and 10.6° at
24.64**, so those are excluded.

- The aliasing is not an escape: the phase is defined mod π and dt = 0.125, so aliases are π/dt = 25.13 km/s/kpc
  apart. The alternatives to 26.10 are 0.97 and 51.2, both far outside any sensible range.
- Snapshot times were verified from the files' own `Time` field: exactly 0.125 model units apart. The units are
  self-consistent too — v_c(8 kpc) = 221.1 km/s and asymmetric drift 9.4 km/s predict a mean stellar v_φ of
  211.7, and the measured value is 211.8.
- **So our two independent methods agree (26.02 Fourier, 26.10 phase) and the simulation's table is low**, by
  +0.91 on average against its Dehnen column and +0.63 against its finite-difference one, with 0.47 rms scatter.
  The mean ratio is 1.021, tantalisingly close to 1/0.9778 = 1.0227 (the model-time-to-Gyr factor), but the
  scatter is ±3% so that is suggestive, not a demonstration — the generating script is not in the directory.
- **Effect on the Test 2 results: negligible.** The truth used throughout is 26.017 against the 26.10 measured
  here, a 0.3% difference, where the biases under discussion are 2–7%.

**2. Why are the Ω_b error bars so small?** Because the tempering assumes the model fits.
`logLs = −χ²/(2√(2 N_data))` normalises χ² by √(2N), its standard deviation *when χ²/N = 1*. It is not.

| run | χ²/N at truth | χ²/N at best | Ω_b width | width × √(χ²/N) | bias | bias / inflated σ |
|---|---|---|---|---|---|---|
| t = 2 | 2.16 | 2.01 | 1.76 | 2.50 | −0.85 | 0.34 |
| t = 3 | 2.24 | 2.24 | 1.91 | 2.86 | +1.10 | 0.38 |
| t = 4 | 2.24 | 2.16 | 1.73 | 2.55 | +0.90 | 0.35 |
| t = 5.5 | 2.35 | 2.35 | 1.13 | 1.74 | −0.78 | 0.45 |
| t = 7 | 2.49 | 2.37 | 0.91 | 1.39 | −0.09 | 0.07 |

(scheme B, `sf`; width = half the 16–84% range.)

- **χ²/N ≈ 2.0–2.4 at the best design point**, so the model does not fit within the bootstrap errors. The usual
  correction is to scale the covariance by χ²_min/N, equivalently χ² → χ²/c, which flattens the curvature by c
  and widens every parameter by **√c ≈ 1.4–1.55**.
- **With that correction every run is consistent with the truth: 0.07–0.45σ.** The apparent >1σ biases are an
  artefact of treating a misspecified model's curvature as if it were statistical.
- **The two runs with the smallest error bars have the largest χ²/N** (t = 5.5 and t = 7: widths 1.13 and 0.91,
  χ²/N 2.35 and 2.49). That is the signature of misspecification rather than of better data — a worse fit makes
  the χ² surface steeper, so the uncorrected posterior is *narrower*, not wider.
- Under scheme A, χ²/N was 7–8, so √c ≈ 2.7; that is a large part of why its posteriors looked so confidently
  wrong.
- **Recommendation:** quote Ω_b errors inflated by √(χ²_min/N), or re-temper with
  `logLs = −χ²/(2 c √(2N))`. Figure: `figures/test2_pattern_speed_recovery_schemeB_inflated.png`
  (`plot_test2_pattern_speed_recovery.py --inflate`).
- **What is left is the residual χ²/N ≈ 2 itself** — the hot library, the sparse-block bootstrap covariance
  (a known 1.25 instead of 1 at |z| 1–2 kpc), and the φ_b quadrupole-node shift. Reducing it is what would
  genuinely tighten the constraint, rather than widening the error bar to cover the bias.

## 2026-11-01 — Emulator audit against the literature: the error bars are narrow because the GP error is not propagated

Follow-up to the entry below, which attributed the narrow Ω_b posteriors to χ²/N > 1 and a √(χ²_min/N) rescaling.
That factor is a **convention**, not a derivation — it is what `scipy.curve_fit` does by default
(`absolute_sigma=False`) and what the PDG uses for discrepant averages, and it assumes the excess χ² comes from
uniformly underestimated errors. If the excess is a localised model defect it is the wrong correction. It should
not have been called "the usual correction".

A second cause is better founded and is standard practice in the emulator literature:

| run | GP LOO rms | GP predictive σ at truth | 1σ of the tempered logL |
|---|---|---|---|
| t = 2 | 0.82 | 1.55 | 0.5 |
| t = 3 | 0.86 | 1.15 | 0.5 |
| t = 4 | 0.78 | 0.82 | 0.5 |
| t = 5.5 | 1.08 | 1.02 | 0.5 |
| t = 7 | 1.07 | 1.04 | 0.5 |

- **The emulator's own uncertainty is 1.6–3.1× the logL difference that defines the 1σ contour**, and
  `test1_gp_emulator.logpost` samples `emu.predict(q)[0]` — the GP **mean only**. The GP variance is computed
  (`return_std=True` is used for the slice plots) but never enters the posterior. The quoted widths are therefore
  the width of a surface known to ±1 in units where the 1σ feature is 0.5.
- **How the peer-reviewed emulators differ.** Our design (scrambled Sobol, quadratic mean + Matérn-5/2 with ARD
  length scales and a white-noise term, LOO validation) is conventional — a polynomial mean plus a GP on the
  residual is exactly the Bayes-linear emulator form of Vernon, Goldstein & Bower, and Matérn-5/2 with ARD is the
  standard emulator kernel. Two things are not conventional:
  1. **A single fixed design with no refinement.** Rogers et al. (2019) show that Bayesian optimisation of the
     training set (GP-UCB acquisition) shrinks the 95% credible volume by 90% *with fewer simulations* than a
     Latin hypercube; Vernon & Bower's history matching uses iterative waves that cut the volume by orders of
     magnitude. We evaluate 257 points once, over a box chosen from an earlier run, and never add points where
     the emulator is uncertain. Our trajectories are saved, so a second wave concentrated on the posterior ridge
     would be cheap and is the obvious improvement.
  2. **We emulate logL itself, not the observable.** Most cosmology emulators (Coyote/CosmicEmu, the Lyman-α
     emulators) emulate the *data vector* and then evaluate an exact likelihood. Emulating logL is done
     (Leclercq 2018; the Bayesian-optimisation-for-evidence line of work) but is harder, because logL spans a
     huge dynamic range and GP accuracy degrades away from the peak — which is what the LOO numbers above show.
- **Recommendation, in order of value:** (i) propagate the GP variance into the posterior, e.g. draw
  logL ~ N(mean, σ) inside `logpost` or marginalise over it; (ii) run a second design wave on the posterior
  ridge; (iii) only then argue about χ²/N rescaling.

**Pattern speed, continued.** Asked to use the simulation directory's tabulated Ω_b as the truth, which the
plots `figures/test2_pattern_speed_recovery_schemeB_truth_tabulated{,_fd}.png` now do
(`plot_test2_pattern_speed_recovery.py --truth tabulated`). On those values the mean |Ω_b bias| of the `sf`
weighting rises from 2.1% to 4.8% (−1.1, +7.0, +9.0, +1.8, +5.2% at t = 2, 3, 4, 5.5, 7), i.e. the method looks
worse, not better. Our own measurement was checked further: the m = 2 phase agrees with the inertia-tensor major
axis to 0.2°, a direct image cross-correlation of consecutive snapshots gives Ω_b = 26.05, and the snapshot time
unit was confirmed by integrating 3000 near-circular outer-disc stars in the axisymmetrised potential — dt =
0.1250 model units fits best (5.1° median error) against 0.1278 (7.1°) for the "Time field is in Gyr"
hypothesis. Four independent routes give 26.0–26.1.

## 2026-11-01 — Truth switched to the tabulated Ω_b; emulating the coefficients fails, and why; library noise found

**Truth.** Adopted the simulation directory's **finite-difference** column as the truth (user's decision;
`plot_test2_pattern_speed_recovery.py --truth tabulated_fd`, now the default). Scheme B, `sf` weighting:

| run | truth (fd) | recovered | bias | σ |
|---|---|---|---|---|
| t = 2 | 43.90 | 43.40 | −1.1% | −0.28 |
| t = 3 | 38.67 | 40.63 | +5.1% | +1.02 |
| t = 4 | 33.28 | 34.86 | +4.7% | +0.91 |
| t = 5.5 | 28.21 | 28.06 | −0.5% | −0.13 |
| t = 7 | 25.20 | 25.92 | +2.9% | +0.80 |

Mean |bias| 2.9% (against 2.1% on our own measurement, 4.8% on the tabulated Dehnen column).
Our independent value stands at 26.0–26.1 for t = 7 by four routes, but is not used.
Figure `figures/test2_pattern_speed_recovery_schemeB_fdtruth.png`.

**Emulating the Fourier coefficients instead of logL — tried, and it is worse.**
`scripts/emulate_fourier_coefficients.py` emulates the library coefficient vector (663 numbers for scheme B),
whitened by the data covariance so that PCA truncation error maps directly onto χ², then evaluates χ² exactly.
This is what the cosmology emulators do (emulate the data vector, keep an exact likelihood), so it was expected
to win.

| n PC | variance explained | held-out logLs error rms |
|---|---|---|
| 16 | 75.4% | 2.75 |
| 48 | 84.4% | 2.74 |
| 128 | 95.2% | 2.72 |
| direct logL emulator | — | **1.07** |

- **The error is flat in the number of components**, so it is not PCA truncation: the coefficient variation
  itself is not a smooth function of θ. The unpredictable part is per-point **library sampling noise**, which
  enters χ² = h|e|² quadratically (and with a positive bias) and so cannot be emulated away, whereas the direct
  logL emulator absorbs it into its WhiteKernel.
- Conclusion: **keep emulating logL** for this problem. The usual argument for emulating the data vector assumes
  a noise-free simulator; ours is not.

**The library noise is real, and the likelihood ignores it.** Bootstrapping the 50k orbits 100 times at the
t = 7 truth point (selection weighting, scheme B) gives the library's own coefficient covariance C_lib:

- C_lib / C_data, median over blocks: **0.091** (16–84%: 0.045–0.202).
- χ²/N at the truth: **2.489** with C_data alone, **1.931** with C_data + C_lib.

So neglected library noise accounts for **37% of the excess** of χ²/N above 1. It is a genuine likelihood
misspecification, it is cheap to fix (one orbit bootstrap per run), and fixing it changes the posterior itself,
not just the error bar — it down-weights the blocks where the library is noisiest. The remaining χ²/N ≈ 1.93 is
still unexplained (hot library, sparse-block bootstrap covariance, the φ_b quadrupole shift).

**On reframing as simulation-based inference.** Not the right tool here: the likelihood is explicit (Gaussian in
the coefficients with a bootstrap covariance), the simulator is deterministic in θ given fixed initial
conditions, and NPE/NLE/NRE want 10³–10⁵ simulations against our 257 at 150 s each. The useful part of the SBI
instinct — "characterise the simulator's own stochasticity" — is precisely the C_lib measurement above, which
costs one bootstrap rather than 10⁴ integrations.

## 2026-11-01 — Emulator error is NOT why the Ω_b errors look small (negative result; retracts the entry below)

The entry below claimed the narrow Ω_b posteriors came from not propagating the GP uncertainty, on the grounds
that the GP's predictive σ ≈ 1.0 exceeds the 0.5 that defines the 1σ contour. That reasoning was wrong and the
measurement says so. Recorded so nobody repeats it.

**What the GP's σ actually is.** The fitted kernel for `test2_t7_sf_10pct` (scheme B, `sf`) is

    103² × Matérn(length_scale=[1.97, 6.89, 1.53, 1.73], ν=2.5) + WhiteKernel(0.731)

so σ ≈ 1.05 splits into **white noise 0.855** (scatter in the logL *values*, i.e. library sampling noise — the
same thing measured as C_lib) and **interpolation 0.608** (genuine "not sampled here"). `return_std` is the
predictive σ for a *new observation*; it is not the uncertainty of the fitted surface, which is far smaller
because 256 points average the noise down.

**Three ways to propagate it, all tested, all negligible:**

| treatment | Ω_b width |
|---|---|
| GP mean only (current) | 0.93 |
| + ½σ²_total (log-normal mean) | 0.94 |
| + ½σ²_interp only | 0.94 |
| 24 refits on logL jittered by σ_white, pooled | 0.94 (**1.01×**) |

The realisation test also gives a median scatter of only **±0.07** against a width of 0.93. Jittering 256
points by 0.855 and refitting moves the surface by ≈ 0.855/√256 ≈ 0.05, which is exactly what is seen.

**Why the naive argument failed.** logL spans ≈ 100 across the design box (kernel amplitude 103) while the
emulator uncertainty is ≈ 0.7. Near the peak the design is dense, so adding ½σ² lifts only regions that are
already excluded, and the posterior barely moves.

**So the error bars are not obviously too small.** The tempering `−χ²/(2√(2N))` already inflates parameter
widths by √(√(2N)) ≈ 6 over the naive Δχ² = 1 criterion. The problem is not the width but the **bias**: the
systematic offset in Ω_b (2–5%) is comparable to the tempered width (3–4%), which is a model-misspecification
problem, not an error-bar one. Of the two candidate corrections raised earlier, √(χ²_min/N) is a fitting
convention rather than a derivation, and emulator-error propagation is now measured to be a 1% effect. Neither
should be used to explain away the bias.

**Still worth doing** (for correctness, not for the error bars): add C_lib to the covariance — it removes 37%
of the χ² excess and changes the posterior, not just its width.

## 2026-11-01 — Is sqrt(2N) the right temperature? No: it is ~4x too wide, and it hides the bias

`scripts/calibrate_posterior_temperature.py` (uses saved coefficients only, no integrations).

**The stated justification for sqrt(2N) does not hold.** sqrt(2N) is the standard deviation of the *absolute*
chi2 over data realisations. A posterior width depends instead on the scatter of the *difference* between two
models on the *same* data, where the data fluctuation is common to both and cancels. For predictions m1, m2:

    Delta^2 = (m1-m2)^T C^-1 (m1-m2),   E[dchi2] = +/- Delta^2,   std[dchi2] = 2 Delta

so the noise scale is **2*Delta, set by the model separation, not by N**.

**Numerically it happens to land close.** For the 20 design points nearest the truth (t = 7, scheme B, `sf`),
Delta = 17.9–22.7, so 2*Delta = 35.9–45.4 against sqrt(2N) = 36.4. That is a coincidence of this design's
density, not a principle — halve the design spacing and 2*Delta halves while sqrt(2N) does not move.

**What the statistical error actually is.** Fisher matrix from the saved library coefficients (local linear fit
of the whitened predictions over the 60 nearest design points, marginal over the other three parameters):

| run | Fisher sigma(Omega_b) | tempered posterior width | ratio | observed bias |
|---|---|---|---|---|
| t = 7, SF | **0.234** (0.90%) | 0.93 | **4.0x** | +0.72 |
| t = 5.5, SF | **0.273** (0.95%) | 1.13 | **4.1x** | −0.15 |

**Coverage test.** 400 synthetic data sets drawn from the truth model, posterior built from the 256 Sobol
points by importance weighting:

| temperature | coverage of the 68% interval | median width |
|---|---|---|
| 8 | 41% | 0.07 |
| 16 | 87% | 0.28 |
| **36.4 = sqrt(2N)** | **100%** | 0.85 |
| 64 | 100% | 1.40 |

- **The current tempering is far too conservative**, not too tight: 100% coverage where 68% is wanted.
- *Caveat:* with only 256 discrete points a narrow posterior is under-resolved, so the widths below T ≈ 16 are
  underestimates and the low-T coverage numbers are not trustworthy. The T = 36.4 row is reliable (its width,
  0.85, matches the 0.93 that emcee on the GP gives).
- **The consequence matters: the tempering is absorbing the systematic into the error bar.** The t = 7 bias of
  +0.72 is **3.1 sigma** against the Fisher error but only 0.8 sigma against the tempered width. Shrinking the
  temperature to its statistically correct value would not fix the method - it would reveal the bias as highly
  significant.

**How to get a faithful posterior (recommendation).** Separate the two errors instead of conflating them:
1. **Fix the covariance, then use T = 1.** Add C_lib (worth 37% of the chi2 excess) and address the sparse-block
   bootstrap miscalibration. The width is then the statistical one, about 0.9% in Omega_b.
2. **Measure the systematic separately** by running the method over several mock realisations and quoting the
   scatter/offset as a systematic error, rather than inflating the statistical error until it covers.

**One real correlation issue the user raised.** The covariance is block-diagonal. For the *data* that is
correct: each star belongs to exactly one (ring, layer) block, so the bootstrap noise is independent between
blocks. For the *library* it is not: every orbit visits many blocks, so C_lib has cross-block terms that the
per-block bootstrap measured here does not capture. Any C_lib added to the likelihood should be estimated with
those cross-block terms, or the correlated part of the library noise will still be unaccounted for.

## 2026-11-01 — Error budget for Ω_b, S and bar_frac; and the bar angle from the v_R quadrupole instead

**Error budget** (`scripts/quantify_test2_errors.py`, scheme B, `sf`, truth = tabulated finite-difference
column, mean over the five snapshots). Statistical = Fisher 1σ of the likelihood itself (correct covariance,
T = 1); systematic = mean(recovered − truth) across the five runs, ± scatter/√5.

| parameter | statistical | tempered (quoted) | ratio | systematic | sys/stat |
|---|---|---|---|---|---|
| Ω_b | 0.317 (0.9%) | 1.489 | 4.7× | **+0.72 ± 0.47** (+2.1%) | 2.3 |
| φ_b | 1.205° | 5.838° | 4.8× | **−4.76 ± 1.01°** | **3.9** |
| S | 3.4% | 16.0% | 4.7× | **+0.8 ± 2.5%** | **0.2** |
| bar_frac | 6.4% | 32.6% | 5.1× | **−11.7 ± 2.8%** | 1.8 |

- The tempering inflates every parameter by a consistent 4.7–5.1×, confirming it is a global temperature and
  not a per-parameter effect.
- **S is the one parameter recovered without bias** (+0.8 ± 2.5%, sys/stat 0.2).
- **bar_frac is biased low by 12%** and Ω_b high by 2%; both are 2σ-level systematics.
- **φ_b is the worst by far** — but *not* through degeneracy: its statistical precision is excellent (1.2°).
  It is a −4.8° bias, 3.9× the statistical error. Written to
  `results/test2_error_budget_sf_schemeB.json`.

**The bar angle does not need the orbit library.** A bar drives a quadrupole in ⟨v_R⟩ whose node lies along its
major axis, so φ_b is directly visible — the distance error just rotates the node, and if σ_DM is known that
rotation can be calibrated out. `scripts/plot_bar_angle_vs_distance_error.py`, t = 7, 10 independent error
draws per level, node fitted per ring with the per-star least squares (so the selection does not bias it):

| fractional distance error | node [deg] | scatter | shift from the true −27° |
|---|---|---|---|
| 0% | −23.90 | — | +3.10 |
| 2.5% | −24.68 | 2.28 | +2.32 |
| 5% | −22.11 | 2.53 | +4.89 |
| 7.5% | −20.48 | 3.97 | +6.52 |
| 10% | −17.67 | 2.81 | +9.33 |
| 12.5% | −14.26 | 1.93 | +12.74 |
| 15% | −12.05 | 1.69 | +14.95 |
| 20% | −8.69 | 2.21 | +18.31 |

- **The relation is linear and steep:** shift = 1.4° + 0.81° per 1% distance error (fit over 0–12.5%), with a
  1.07° rms residual about the line.
- **The selection function barely matters** — the with/without curves agree to about 2°, well inside the
  realisation scatter, which is the same conclusion the posteriors reached.
- **Measured against true positions the node does not move** (−23.9° at every error level, by construction),
  confirming the whole shift comes from the distance smearing.
- **Precision of the calibrated angle: about 2.3°** (the single-realisation scatter), against the library
  method's −4.8° bias. So this is the better route to φ_b, with two caveats:
  1. there is an **intrinsic +1.4 to +3.1° offset at zero distance error** — the quadrupole node is not exactly
     the bar's density axis, so that offset has to be calibrated from a simulation and is not model-free;
  2. σ_DM must be known well: a 1% error in the assumed distance error moves the recovered angle by 0.81°.
- Figure `figures/test2_bar_angle_vs_distance_error_nbody.png`,
  table `results/nbody_bar_angle_vs_distance_error.csv`.

## 2026-11-01 — The bar-angle calibration is NOT universal: it depends on the true bar angle

`scripts/plot_bar_angle_calibration_universality.py`. The galaxy is rotated about z so the bar sits at a
chosen angle to the Sun–GC line, re-observed from the fixed Sun, smeared, and the v_R quadrupole node
remeasured. The same 50k stars are used at every angle (the R < 8 kpc pool is rotation-invariant), so only the
geometry relative to the Sun changes — which is exactly what the line-of-sight distance error couples to.

Node − true bar angle [deg], 5 error draws per cell:

| distance error | 15° | 20° | 25° | 27° | 30° | 35° |
|---|---|---|---|---|---|---|
| 0% | 3.10 | 3.10 | 3.10 | 3.10 | 3.10 | 3.10 |
| 5% | 2.96 | 3.71 | 4.25 | 5.16 | 5.67 | 4.79 |
| 10% | 6.60 | 7.97 | 10.12 | 10.22 | 11.61 | 12.58 |
| 15% | 8.72 | 11.06 | 14.08 | 15.03 | 17.34 | 18.35 |
| 20% | 9.27 | 13.58 | 16.81 | 17.17 | 19.11 | 23.20 |

- **The bias roughly doubles from 15° to 35°** (6.6° → 12.6° at a 10% distance error). The spread across bar
  angles is 4.3° median and 13.9° maximum, against a single-realisation scatter of only 1.8°, so the dependence
  is real, not noise. A single calibration curve therefore cannot be applied without already knowing the angle.
- At 0% the shift is identical (3.10°) at every angle, as it must be for the same stars rigidly rotated — a
  useful check that the rotation and the measurement are consistent.
- **It is still invertible, just less sensitive.** The measured node remains monotonic in the true angle, but
  the gradient flattens as the error grows:

  | distance error | d(node)/d(φ_true) | φ error per 1° of node error |
  |---|---|---|
  | 0% | −1.00 | 1.00 |
  | 5% | −0.88 | 1.13 |
  | 10% | −0.69 | 1.45 |
  | 15% | −0.49 | 2.05 |
  | 20% | −0.34 | 2.97 |

- **At a 10% distance error, φ_b is recoverable to about 2.7°** (node scatter 1.8° × 1.45), against the
  library's −4.8° ± 1.0° bias. So the quadrupole route still wins, but it needs the full two-dimensional
  calibration surface node(φ_true, σ_DM) rather than a single curve, and it degrades quickly beyond ~15%.
- Figure `figures/test2_bar_angle_calibration_universality_nbody.png`,
  table `results/nbody_bar_angle_calibration_universality.csv`.

**Correction to the previous entry.** With 10 realisations now run for the selection-function curve too, the
no-selection and with-selection calibrations agree within 1–2° at every error level, inside their ±1.7–3.0°
scatter, so "the selection function barely matters" is now supported rather than asserted. The claimed
"intrinsic +1.4 to +3.1° offset at zero error" was overstated: it is a property of the particular 50k draw
(scatter ±2.3°), and on all 1,493,610 parent particles the zero-error offset is only **+0.8°**.

## 2026-11-01 — The inversion curve: observed quadrupole node against true bar angle at 10% distance error

`scripts/plot_bar_angle_inversion.py`, t = 7, bar rotated to 15–35 deg in 2.5 deg steps, 15 distance-error
realisations per angle. The selection function is applied as a property of the *survey* — the same S(l, b, d)
tables at every bar angle, since a footprint does not know where the bar is — giving ~50,100 selected stars
from the full 2 M each time.

| true angle | node, no selection | node, with selection |
|---|---|---|
| 15° | −9.02 ± 1.91 | −8.88 ± 2.18 |
| 20° | −13.06 ± 2.75 | −12.02 ± 2.35 |
| 25° | −16.03 ± 2.66 | −15.77 ± 1.92 |
| 27.5° | −18.50 ± 3.12 | −16.20 ± 1.38 |
| 30° | −20.15 ± 3.30 | −18.81 ± 1.54 |
| 35° | −21.89 ± 3.85 | −22.16 ± 2.03 |

| | fitted relation | node scatter | φ_b recoverable to |
|---|---|---|---|
| no selection | node = −0.664 φ_true + 0.42 | 2.75° | **4.1°** |
| with selection | node = −0.669 φ_true + 1.29 | 2.14° | **3.2°** |

- **The relation is flatter than 1:1** — slope −0.66 to −0.67, not −1. The node under-reports the bar angle
  (a true 35° reads about −22°), and because the mapping is compressed a given node error becomes a 1.5×
  larger angle error. The node is therefore not a direct readout at 10% error; it has to be inverted through
  this curve.
- **The selection function makes no significant difference** to either the slope (−0.664 vs −0.669) or the
  scatter, consistent with everything else we have found about it.
- **φ_b is recoverable to 3–4°**, against the orbit library's −4.8° ± 1.0° bias. Better, but not dramatically,
  and the calibration is specific to this galaxy and this error level.

**Correction.** An earlier 5-realisation run of the same script gave scatters of 3.44° (no selection) and
1.29° (with), implying 5.6° and 1.9°, and suggested the selection case was much the better one. With 15
realisations the two agree (2.75° and 2.14°). A standard deviation from 5 samples carries about 30%
uncertainty, so that apparent factor of 2.7 was noise and should not have been quoted as a difference.

## 2026-11-01 — Does the bar-angle calibration transfer between snapshots? Mostly, with a −1.9° bias

`scripts/plot_bar_angle_transfer_test.py`. The relation fitted at t = 7 is applied to the measured quadrupole
nodes of t = 4 — a genuinely different bar (A₂ 0.472 vs 0.397, R_bar 4.25 vs 5.00 kpc, Ω_b 33.3 vs 25.2), both
with 10% distance errors and the DR3-Mira-like selection, 15 realisations per angle.

| | fitted relation |
|---|---|
| t = 7 (the calibration) | node = **−0.669** φ_true **+1.29** |
| t = 4 (its own relation) | node = **−0.639** φ_true **+1.79** |

Applying the t = 7 calibration to the t = 4 nodes:

| true φ_b | node | inferred φ_b | error |
|---|---|---|---|
| 15.0 | −7.77 | 13.54 | −1.46 |
| 20.0 | −10.92 | 18.24 | −1.76 |
| 25.0 | −14.22 | 23.17 | −1.83 |
| 30.0 | −17.45 | 28.00 | −2.00 |
| 35.0 | −20.60 | 32.71 | −2.29 |

- **Mean error −1.86°, rms 1.89°**, against a per-measurement scatter of 2.20°. So a single measurement would
  be consistent with the truth at about 1σ — the calibration does transfer in that practical sense.
- **But the offset is a genuine systematic, not noise.** All nine residuals have the same sign and their
  scatter about the mean is only 0.34°, so it will not average down with more stars or more realisations. The
  slopes differ by 4.5% (−0.669 vs −0.639) and the intercepts by 0.5°, and together these produce the −1.9°
  shift across 15–35°.
- **Budget for the transferred method:** about 1.9° systematic plus 2.2° statistical, so **≈2.9° in
  quadrature**, against the orbit library's −4.8° ± 1.0° bias. Still the better route, but the margin is
  narrower than the 3.2° within-snapshot figure implied.
- **Implication for real data.** A calibration taken from a simulation and applied to the Milky Way would
  carry a systematic of at least this order, and plausibly more: t = 4 and t = 7 are the same galaxy 2.9 Gyr
  apart, which is far closer to each other than either is to the MW. The quadrupole route is not calibration-free.
- Figure `figures/test2_bar_angle_transfer_nbody_to_nbody_t4.png`;
  tables `results/nbody{,_t4}_bar_angle_inversion_10pct.csv`.

## 2026-11-01 — The transfer offset IS the selection function interacting with a different bar

Follow-up to the transfer test. The −1.86° offset when the t = 7 bar-angle calibration is applied to t = 4 was
attributed to "the bars being different". Re-running the same transfer with the selection function switched
off separates the two causes:

| transfer t = 7 → t = 4 | mean error | rms | slopes (t7 / t4) | intercepts |
|---|---|---|---|---|
| **with** selection | **−1.86°** | 1.89° | −0.669 / −0.639 | +1.29 / +1.79 |
| **without** selection | **−0.35°** | 0.47° | −0.664 / −0.645 | +0.42 / +0.16 |

- **Without the selection the calibration transfers almost perfectly** (−0.35°, residuals −0.07 to −1.11°),
  even though the two bars are genuinely different (A₂ 0.472 vs 0.397, R_bar 4.25 vs 5.00 kpc, Ω_b 33.3 vs
  25.2). The intrinsic difference between the bars costs only ~0.35°.
- **The selection contributes the other −1.5°**, and roughly doubles both the slope difference (2.9% → 4.5%)
  and the intercept difference (0.26° → 0.50°). The same fixed S(l, b, d) selects a different *part* of each
  bar, because the two have different density distributions, so it weights the kinematic maps differently and
  shifts the node by a different amount.
- **This reconciles two statements that looked contradictory.** At a *fixed* galaxy the selection is almost
  irrelevant — it adds +0.54° to the bias and leaves the calibration slope unchanged to 1%. *Across* galaxies
  it is the dominant term in the transfer error. So the selection is harmless when calibrating and applying on
  the same galaxy, and harmful exactly in the configuration real data forces: calibrate on a simulation, apply
  to the Milky Way, whose bar differs from any snapshot by far more than t = 4 differs from t = 7.
- Practical consequence: the calibration systematic cannot be reduced by choosing a snapshot with the right
  pattern speed; it is set by how well the simulation reproduces the bar's *density distribution* as seen
  through the survey's selection.

## 2026-11-01 — Testing Zhang et al. (2024): <|v_R|/v_tot> beats the v_R quadrupole, but is not galaxy-free

Zhang et al. (2024) propose the m = 2 phase of <|v_R|/v_tot> as a bar-orientation proxy, more robust to
distance error than a <v_R> map; the GaiaNIR follow-up (arXiv:2604.18329) quotes 1–7° accuracy.
`scripts/compare_bar_angle_estimators.py` and `scripts/test_estimator_transfer.py` test it here.

**1. Robustness to distance error at a fixed galaxy — the claim holds, strongly.** t = 7, native bar angle,
with the selection function, 15 realisations:

| | median drift over 5–20% | median scatter |
|---|---|---|
| ⟨v_R⟩ quadrupole node | **9.06°** | 3.08° |
| ⟨\|v_R\|/v_tot⟩ m = 2 phase | **0.55°** | 1.79° |

Better on accuracy *and* precision. The raw phase (true axis −27.0°) reads −27.71, −27.76, −27.80, −26.89 at
0, 5, 10, 15% error — within 0.9° of the truth with **no calibration at all** — and only breaks down at 20%
(−23.69).

**Why.** ⟨v_R⟩ is *signed*: its inflow and outflow lobes partially cancel when distance smearing mixes
azimuths, and because the smearing is along the line of sight rather than isotropic the cancellation is
asymmetric, which rotates the node. ⟨\|v_R\|/v_tot⟩ is unsigned — it tracks σ_R/v_tot, which peaks along the
bar major axis — so mixing dilutes the amplitude without turning the phase.

**2. But the direct read-out is galaxy-dependent.** Taking the measured phase as the bar angle, 10% error plus
selection, over φ_true = 15–35°:

| | mean error | rms | max |
|---|---|---|---|
| t = 7, ⟨v_R⟩ node | +8.41° | 8.59° | 11.05° |
| **t = 7, ⟨\|v_R\|/v_tot⟩** | **+0.44°** | **0.62°** | **0.96°** |
| t = 4, ⟨v_R⟩ node | +10.97° | 11.28° | 14.63° |
| **t = 4, ⟨\|v_R\|/v_tot⟩** | **+4.44°** | 4.59° | 6.22° |

At t = 7 the estimator is essentially exact without calibration; at t = 4 it over-reports by 3–6°. The fitted
relations show why — the slope is near −1 but not equal to it, and it differs between snapshots:

| estimator | t = 7 slope | t = 4 slope |
|---|---|---|
| ⟨v_R⟩ node | −0.731 | −0.594 |
| ⟨\|v_R\|/v_tot⟩ | **−0.947** | **−0.825** |

**3. Transfer between snapshots: no better than the quadrupole.** Applying the t = 7 relation to t = 4 gives a
mean error of −3.49° (⟨v_R⟩) and −4.23° (⟨\|v_R\|/v_tot⟩). The slope *ratio* favours
⟨\|v_R\|/v_tot⟩ (0.87 vs 0.81, i.e. less compression) but the intercept difference works against it, and the
two end up comparable. So the estimator removes the distance-error sensitivity but **not** the galaxy
dependence — which, per the entry above, is the selection function interacting with each bar's density
distribution.

**Caveats.** The 9 angles in each scan share the same 10 distance-error seeds, so they are correlated and the
formal slope errors (±0.01–0.02) understate the true uncertainty; an earlier independent run gave −0.669 for
the t = 7 ⟨v_R⟩ slope against −0.731 here, a ~1.8σ difference by the formal errors. The transfer numbers in
particular carry roughly ±1.5° of realisation noise, so −3.49 vs −4.23 should not be read as a real ordering.
What is robust across runs is the factor ~16–22 difference in distance-error drift and the slope being much
closer to −1 for ⟨\|v_R\|/v_tot⟩.

Figures `figures/test2_bar_angle_estimators_nbody{,_sf}.png`,
`figures/test2_estimator_transfer_nbody_to_nbody_t4.png`.

## 2026-10-31 — Estimating the selection function from the data instead of knowing it (BACKFILLED log entry)

`scripts/selection.py` (`empirical_tables`, `ratio_tables`), `scripts/rescore_empirical_selection.py`.
The `sf` weighting used the S tables stored in the mock, i.e. perfect knowledge of the selection, which no
survey has. Two estimators built from the data alone, both rescoring the finished searches from saved
trajectories (no integrations):

- **`sfratio`** — the direct ratio ρ_obs/ρ_lib per (HEALPix nside 16 × 1 kpc distance) cell (user's proposal).
- **`sfemp`** — a separable S = s_sky × s_dist × s_R fitted to the data's three margins by weighted IPF.

Settings came from a scan at the truth: nside 32 and a **single** IPF pass for `sfemp`; iterating to
self-consistency through the 1/τ feedback makes it worse (correlation with the true S falls 0.73 → 0.55).

| how S is obtained | mean \|Ω_b bias\| | mean worst offset | mean LOO |
|---|---|---|---|
| none (`asis`) | 7.3% | 2.6σ | 1.34 |
| **ρ_obs/ρ_lib** (`sfratio`) | **3.6%** | 1.5σ | 1.05 |
| separable fit (`sfemp`) | 4.0% | 1.5σ | 0.94 |
| true S (`sf`) | 2.1% | 1.0σ | 0.92 |

- **An estimated S recovers about 70% of the gain from knowing S exactly**, and the direct ratio beats the
  separable fit in 4 of 5 snapshots.
- **The better estimate of S is the worse estimator of the parameters.** The separable fit reproduces the true
  S much better (correlation 0.73 against 0.48) yet recovers Ω_b slightly worse. The likelihood does not need
  S to be right, only right enough to fix the weighting within each block — which is why the noisy cell-by-cell
  ratio is good enough.
- Posteriors and a README in `figures/test2_schemeB_selection/` (20 corner plots, four treatments × five
  snapshots).

## 2026-11-01 — Combining ring phases by amplitude, not by median: the v_R bias becomes universal (BACKFILLED)

The per-ring m = 2 phases were being combined with a plain **median over a fixed 1–4 kpc window**, which gives
a ring with no bar signal the same vote as one with a strong signal. Replacing that with the **sum of the
complex m = 2 coefficients** (each ring weights itself by its own amplitude) changes the conclusions.

Drift of the measured angle from 0 to 10% distance error, selection on, 8 realisations:

| scheme | t = 2 | t = 3 | t = 4 | t = 5.5 | t = 7 | spread |
|---|---|---|---|---|---|---|
| ⟨\|v_R\|/v_tot⟩, median 1–4 | +4.87 | +7.40 | +4.17 | +2.97 | −0.37 | 7.77 |
| ⟨\|v_R\|/v_tot⟩, vecsum >30% peak | +3.85 | +5.98 | +3.95 | +3.09 | +2.06 | **3.92** |
| ⟨v_R⟩ node, median 1–4 | +9.95 | +8.88 | +9.18 | +9.05 | +11.08 | 2.20 |
| **⟨v_R⟩ node, vecsum 1–4** | +10.77 | +10.45 | +10.36 | +9.58 | +10.61 | **1.19** |

- Weighting halves the snapshot spread for ⟨\|v_R\|/v_tot⟩ — the user's suggestion was right that signal-free
  rings were voting.
- **It also removes t = 7's apparent immunity** (−0.37 → +2.06). The earlier "⟨\|v_R\|/v_tot⟩ is unbiased"
  result was an artefact of the median combination at one unrepresentative snapshot.
- **The decisive number is the scatter across snapshots, because that is what survives calibration:**
  ⟨v_R⟩ node **+10.35° ± 0.46°**, ⟨\|v_R\|/v_tot⟩ +4.01° ± 2.52°. The ±0.46° is comparable to the individual
  measurement errors, so the ⟨v_R⟩ drift is consistent with being *perfectly* universal across a factor-1.8
  range of pattern speed. **A large known bias beats a small unknown one**, so the ⟨v_R⟩ quadrupole is the
  better estimator after all — reversing the recommendation of the entry below.

**Verification** (`scripts/bar_angle_vecsum_verification.py`, 12 realisations, independent seeds per point):

| estimator | direction | mean error | rms |
|---|---|---|---|
| **⟨v_R⟩ node** | t = 7 → t = 4 | **−0.89°** | 1.05 |
| **⟨v_R⟩ node** | t = 4 → t = 7 | **+0.93°** | 1.25 |
| ⟨\|v_R\|/v_tot⟩ | t = 7 → t = 4 | −3.71° | 3.79 |
| ⟨\|v_R\|/v_tot⟩ | t = 4 → t = 7 | +4.16° | 4.26 |

The errors are antisymmetric under swapping source and target, as a real transfer bias must be — a useful
internal check. **With amplitude weighting the ⟨v_R⟩ calibration transfers between snapshots to ~1°**, against
the 1.9–3.5° found with the median combination.
Figures `figures/test2_bar_angle_vecsum_{drift,transfer}.png`.

**Methodological fix that made this measurable.** The earlier scripts used the *same* distance-error seeds at
every bar angle, correlating the points. Measured effect on the fitted slope: spread 0.125 across four seed
choices, against a formal fit error of 0.019 — 6.6×. With independent seeds per angle the spread halves to
0.064. This is why two transfer runs of the same quantity disagreed (−1.86 vs −3.49) and why slopes quoted to
three decimals earlier in this log should be read as ±0.05.

## 2026-11-01 — Why the <|v_R|/v_tot> drift differs between snapshots (BACKFILLED; inconclusive on cause)

With the median combination, the ⟨\|v_R\|/v_tot⟩ drift at 10% ranged from −1.05° (t = 7) to +5.98° (t = 3).
Measured bar properties against that drift, all five snapshots:

| snapshot | A₂ | R_bar | kinematic twist (1–4 kpc) | drift 0→10% |
|---|---|---|---|---|
| t = 2 | 0.472 | 4.25 | 17.2° | +3.86 ± 0.46 |
| t = 3 | 0.446 | 4.00 | 10.1° | +5.98 ± 0.79 |
| t = 4 | 0.472 | 4.25 | 9.1° | +2.70 ± 0.34 |
| t = 5.5 | 0.448 | 4.50 | 4.7° | +1.22 ± 0.64 |
| t = 7 | 0.397 | 5.00 | 4.8° | −1.05 ± 0.43 |

| correlation with drift | |
|---|---|
| bar length R_bar | **−0.96** |
| kinematic twist | +0.65 |
| bar strength A₂ | +0.66 |

- "Kinematic twist" is the peak-to-peak variation of the per-ring ⟨\|v_R\|/v_tot⟩ phase over 1–4 kpc at zero
  distance error. t = 7's bar is almost untwisted (−27.3 to −26.1°, 2.1° with finer binning); t = 4's runs from
  −30.4° at 1 kpc to −20.3° at 4 kpc. Smearing mixes radii, so a twisted bar has its phase rotated while a flat
  one does not — the natural mechanism.
- **But the cause is not established.** These five snapshots are one secular sequence, so R_bar, A₂, twist and
  Ω_b all co-vary; with n = 5 and mutually confounded predictors, −0.96 against +0.65 does not single out bar
  length. The obvious mechanism — the fixed 1–4 kpc window sitting at different fractions of the bar — was
  tested and **ruled out**: rescaling to 0.2–0.8 R_bar changed t = 7 from −1.17 to −1.17 and t = 4 from +2.97
  to +2.73.
- Superseded in practice by the amplitude-weighted combination (entry above), which shrinks the spread and
  makes the ⟨v_R⟩ node the preferred estimator regardless.

## 2026-11-02/03 — Test 2 at a 12.5% distance error, and the selection function you can actually have

The 10% campaign repeated at **12.5%** (σ_DM = 0.125·5/ln10 = 0.2714 mag). Nothing about the method
changed: same 257-point focused Sobol design, same DR3-Mira-like selection, same χ² schemes. Six
searches, 1,542 integrations: four on the cluster (t = 2, 3, 4, 5.5), two on the workstation
(t = 7 with and without the selection function). Scripts `run_test2_12p5pct_{local,cluster}.sh`;
`run_test2_posteriors.sh` is now parameterised by `SIGMA`/`TAG` and its `ONLY` is an ordered queue,
so the 10% defaults reproduce the old behaviour exactly.

**All biases below are against the tabulated finite-difference column** (the convention adopted
2026-11-01), scheme B, five snapshots.

| how S is obtained | available for real data? | 10% | 12.5% |
|---|---|---|---|
| none (`asis`) | — control | 5.3% | 9.3% |
| **ρ_obs/ρ_lib (`sfratio`)** | **yes** | **2.0%** | **3.1%** |
| separable fit (`sfemp`) | yes | 2.2% | 3.1% |
| true S (`sf`) | no — upper bound | 2.9% | 4.1% |

Per snapshot at 12.5% (Ω_b bias): t = 2 −4.1%, t = 3 +2.5%, t = 4 +2.3%, t = 5.5 −4.7%, t = 7 −2.0%
for `sfratio`; the 68% half-width grows 1.55 → 1.74 km/s/kpc (×1.13) from 10% to 12.5%.

**The method degrades gracefully**: a 25% larger distance error costs about one point of bias and
13% of posterior width.

**An estimated S beats the exact S — and the earlier claim to the contrary was a truth-convention
artefact.** The 2026-10-31 entry and `figures/test2_schemeB_selection/README.md` report
`sf` 2.1% against `sfratio` 3.6%, i.e. "an estimated S recovers about 70% of the gain from knowing
S exactly". Those numbers are scored against **our own measured Ω_b**, not the fd column. The same
posteriors scored against fd give `sfratio` 2.0% against `sf` 2.9% — the ordering **reverses**.
Both the README and that entry's table should be read as superseded. *(Retraction.)*

The mechanism (user, 2026-11-03): what the likelihood needs is for the library's spatial
distribution to match the data's, not for S to be recovered correctly. ρ_obs/ρ_lib forces exactly
that, cell by cell. The exact S cannot, because the library is not the simulation — it is orbits
integrated in a *trial* potential, time-sampled from t ≥ 0.25 Gyr and folded, so weighting it by the
true S still leaves a density mismatch. This also explains the "better S estimate, worse parameters"
puzzle recorded on 2026-10-31 (`sfemp` reproduces the true S at correlation 0.73 against `sfratio`'s
0.48, yet recovers Ω_b no better): fidelity to S is the wrong target.

**`sf` carries a positive Ω_b bias that the estimated S removes**: at 12.5%, t = 3 +6.0%, t = 4
+9.5%, t = 7 +2.6% with the true S, against +2.5%, +2.3%, −2.0% with the ratio. The t = 4 outlier
reported mid-campaign (+4.7% → +9.5% from 10% to 12.5%) is largely a known-S artefact — with the
ratio it is +0.2% → +2.3%.

**Trajectory time cadence: 50 Myr costs nothing measurable** (`scripts/test_time_cadence.py`,
rescoring finished searches on `traj[:, ::k]`). Storing every other snapshot halves the trajectory
files, which are 98% of a run's bytes.

| t = 2, scheme B, `sf` | Ω_b | bias |
|---|---|---|
| 25 Myr (30 library samples) | 42.96 [41.26–44.91] | −2.1% |
| 50 Myr (15) | 43.14 [41.32–45.19] | −1.7% |
| 100 Myr (8) | 42.87 [40.89–45.15] | −2.3% |

Medians move by < 0.2 km/s/kpc, a tenth of the error bar; the emulator's LOO rms is unchanged. What
*does* grow is χ²/N at the truth — +3.6% at t = 7 but **+12.4% at t = 2**, the fast bar, where a
coarser sampling thins the phase coverage more. The posterior barely notices because the extra
library noise is close to common across design points. `test2b_qmc_distance_error_run.py --cadence-step K`
thins before scoring *and* before saving, so the stored file and any later rescoring share one grid.

**The fd truth column is itself noisy**, which caps how well a per-snapshot bias can be measured:
over the 49 snapshots in 2 ≤ t ≤ 8 the residual about a smooth fit is **0.66 km/s/kpc (2.1%)**, and
adjacent snapshots 122 Myr apart differ by up to 4.15 — far more than a bar's pattern speed can
physically change. Against a mean |bias| of 3.1% (≈1.0 km/s/kpc), roughly two-thirds of what we
quote may be the method and one-third the truth column. Individual outliers deserve correspondingly
less weight.

**In progress.** Seven fill-in snapshots (t = 2.375, 2.625, 3.375, 3.625, 4.5, 5.0, 6.25), two per
interval between the four fastest runs and one between t = 5.5 and 7, launched on the cluster
2026-11-03 12:11 at the 50 Myr cadence. Twelve points from 25.2 to 43.9, largest gap 5.4 → 3.0, the
range deliberately unchanged. Purpose: separate a systematic trend in Ω_b from snapshot scatter,
which n = 5 cannot do. **Caveat that does not go away with more points** — these are one galaxy in a
secular sequence, so Ω_b, A₂ and R_bar co-vary; twelve snapshots is not twelve independent galaxies,
and this gives a systematic *along this sequence*, not a Milky Way systematic.

Figures: `figures/key_results_12p5pct/` (summary, corner plots per snapshot, README),
`figures/test2_error_level_comparison{,_sfratio}.png`, `figures/test2_snapshot_coverage.png`,
`figures/mock_validation/test2_mock_grid_12.5pct.png`. `figures/` itself was re-filed — see
`figures/README.md` and `scripts/organise_figures.py`.

## Figure index (Test 1b, N-body, perfect mock)

Figures are now filed by purpose - see `figures/README.md`. Nothing was renamed, so every name
below is still found by `find figures -name <name>`; most of the Test 1b set is in
`figures/archive/test1/`, the emulator slices and leave-one-out checks in
`figures/archive/emulator_diagnostics/`, and the corner plots in `figures/archive/posteriors/`.

| file | what it shows |
|---|---|
| `test1_nbody_square_maps_bar_length_vs_data.png` | face-on 0.5 kpc square maps of ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ: data, ground truth, S 1.25, S 0.8, data errors, residuals |
| `test1_nbody_square_maps_bar_strength_vs_data.png` | the same for bar_frac 1.4 and 0.7 |
| `test1_nbody_bar_length_strength_radial_profiles.png` | m = 2 amplitudes and azimuthal means vs R for the S/f models (before the 1/N fix) |
| `test1_nbody_bar_length_strength_polar_maps_*.png`, `..._polar_pull_maps.png`, `..._m2_profiles.png` | polar-cell maps and m = 2 profiles of the S/f models |
| `test1_nbody_square_maps_wrong_bar_angle_20deg_vs_27deg.png` | square maps: data and library for trial bar angles 20° (wrong) and 27° (true), and their residuals |
| `test1_nbody_ring_avg_dispersion_residuals_vs_R.png` | ring-averaged σ_R, σ_φ, σ_z: models − data vs R, with the data error band (the "hot library") |
| `test1_nbody_fourier_compression_vs_mmax.png` | Δχ² and tempered ΔlogL of the S/f models vs the highest Fourier mode kept |
| `test1_nbody_emulator_corner_logLs_meansdisp.png` | posterior: current sector likelihood, means + dispersions (tempered) |
| `test1_nbody_emulator_corner_logLs_means.png` | posterior: current sector likelihood, means only (tempered; made before the 1/N fix, which does not affect means) |
| `test1_nbody_emulator_corner_logLs_fourier_vR_m24_vphi_sigR_sigphi_m024.png` | posterior: proposed Fourier scheme |
| `test1_nbody_emulator_corner_logLs_fourier_vR_sigR_sigphi_m24_vphi_m024.png` | posterior: proposed scheme with the σ ring averages dropped |
| `test1_nbody_emulator_corner_logLs_fourier_means_m024_sigmas_m24.png` | posterior: means m ≤ 4 + all σ m = 2–4 (incl. v_z) |
| `test1_nbody_emulator_corner_logLs_fourier_{mmax4,m2to4}_meansdisp.png` | posterior: all six moments, m ≤ 4 or m = 2–4 |
| `test1_nbody_emulator_corner_logLs_fourier_symm_vR_sin24_vphi_sigR_sigphi_cos024.png` | posterior: proposed scheme with symmetry-allowed terms only (no sin for ⟨v_φ⟩/σ, no cos for ⟨v_R⟩) |
| `test1_nbody_logl_corner_logLs_*.png` | log L of the 65 evaluated points (no emulator) |
| `test1_nbody_qmc_focused_emulator_corner_logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024.png` | focused run (172 points), posterior: proposed scheme, Fourier terms fitted to the stars, joint covariance |
| `test1_nbody_qmc_focused_emulator_corner_logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024.png` | focused run, posterior: the same without the σ ring averages |
| `test1_nbody_qmc_focused_emulator_corner_logLs_fourier_R0to6_vR_m24_vphi_sigR_sigphi_m024.png` | focused run, posterior: proposed scheme from the sector fit (control for the fitting method) |
| `test1_nbody_qmc_focused_emulator_corner_logLs_sectors_meansdisp.png` | focused run, posterior: current polar-cell likelihood, means + dispersions |
| `test1_nbody_qmc_focused_logl_corner_<column>.png` | focused run: evaluated log L (upper) and GP profile (lower) for the columns above and `logLs_fourier_stars_R0to8_…` |
| `test1_nbody_qmc_focused_emulator_{slices,loo}_<column>.png` | focused run: 1D emulator slices through truth and maximum; leave-one-out check |
| `test1_nbody_fourier_star_fit_vs_true_maps_absz_{0to1,1to2}.png` | per-star Fourier fit (m ≤ 4) from the 50k mock vs the true N-body maps (1.49 M particles): truth, reconstruction, difference, pulls, truncation |
| `test1_nbody_qmc_focused_bar_length_strength_degeneracy_<column>.png` | why S and bar_frac are weakly constrained: (S, bar_frac) posterior with curves of fixed bar potential; fractional width of the bar m = 2 potential vs R against S, bar_frac, Ω_b; Q_T(R) posterior vs truth |
| `test1_nbody_qmc_focused_fourier_coefficients_data_vs_models_absz_{0to1,1to2}.png` | Fourier coefficients (m = 0, cos/sin 2φ, cos/sin 4φ) of ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ vs R: data with errors against the truth model, two ridge models, a rejected Ω_b model and the all-particle N-body |
| `test1_nbody_distance_error_fourier_coefficients_perfect_vs_smeared_absz_{0to1,1to2}.png` | Test 2b: Fourier coefficients of the perfect and the distance-smeared (σ_DM = 0.265) mock against the truth model and the truth model observed with the same errors |
| `test2b_smeared_initial_conditions_fourier_coefficients_absz_{0to1,1to2}.png` | Test 2b: smeared data against libraries in the true potential from smeared initial conditions (as is; observed again) and from true initial conditions (observed) |
| `test2b_posterior_distance_draws_fourier_coefficients_absz_0to1.png` | Test 2b: smeared data against libraries from posterior distance draws (density prior; full phase-space prior), smeared stars and true positions, all observed with the errors |
| `test2b_qmc_distance_error_7p5_coadded_error_emulator_corner_logLs_{observed,asis}.png` | Test 2b, 7.5% distance error, full focused grid: posterior with the co-added error (library observed again / as is); `..._logl_corner_...` evaluated points; without `_coadded_error` and with the long column name: bootstrap-only error |
| `test2_nbody_mock_selection_vs_dr3_miras_sigmaDM0.163.png` | Test 2 mock: sky map, distance and R distributions and face-on map of the selection-function mock against the DR3 Miras |
| `test1_nbody_density_stationarity.png` | density change of the libraries; counts statistic retracted (needs the selection function) |

## Next steps

Agreed approach (2026-10-27): start from the simplest data and add realism one step at a time.

**Likelihood.** Fit 3D mean and dispersion maps of v_R, v_φ, v_z (later also cross terms).
- Cells are (R, φ mod π, |z|) or 0.5 kpc boxes, with bootstrap errors on the data.
- A per-star local Gaussian with cell-constant μ, C is algebraically the same thing. The map form
  makes the error model explicit (bootstrap errors capture non-Gaussianity) and lets us choose
  which moments and cells enter.

1. **Test 1b (perfect N-body mock):** the global search launched 2026-10-27 14:06 saves the library
   box sums and polar cell moments for every point.
   - Fit the 3D moment maps (means; means + dispersions; with and without R < 1 kpc), plus the
     local Gaussian and the KDE for comparison.
   - Fit a GP emulator and draw the corner plot for each, and check recovery of
     (Ω_b, φ_b, S, bar_frac).
   - Expected bias: bar_frac about +3%, from the spherical halo.
2. **Test 2a (selection function):** a Mira-like N-body mock with a Gaia-like footprint, extinction
   and magnitude limit.
   - Library orbits weighted by 1/τ_i with τ_i = ⟨S(x(t))⟩ along the orbit.
   - Model cell moments computed with library samples weighted by S(x_j).
   - Unfolded cells for the data (S is not symmetric); the library can stay folded.
3. **Test 2b (distance and RV errors):** PL distance errors of about 10–15% plus v_los jitter.
   - Compare (i) no correction, (ii) forward-modelling the errors into the library and comparing
     moments of (μ_l, μ_b, v_los) in (l, b, distance-modulus) cells, as in Portail+17 and
     Clarke & Gerhard 2022, and (iii) seeding several orbits per star from its distance posterior
     (to avoid counting the errors twice when the library starts from noisy initial conditions).
4. **Later:**
   - age-dependent selection (hotter vs colder Mira populations) with Auriga;
   - model mismatch (an N-body bar outside the AGAMA family);
   - DR4 sample size and footprint.
