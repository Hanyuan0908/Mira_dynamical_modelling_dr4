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

## Figure index (Test 1b, N-body, perfect mock)

All in `figures/`.

| file | what it shows |
|---|---|
| `test1_nbody_square_maps_bar_length_vs_data.png` | face-on 0.5 kpc square maps of ⟨v_R⟩, ⟨v_φ⟩, σ_R, σ_φ: data, ground truth, S 1.25, S 0.8, data errors, residuals |
| `test1_nbody_square_maps_bar_strength_vs_data.png` | the same for bar_frac 1.4 and 0.7 |
| `test1_nbody_bar_length_strength_radial_profiles.png` | m = 2 amplitudes and azimuthal means vs R for the S/f models (before the 1/N fix) |
| `test1_nbody_bar_length_strength_polar_maps_*.png`, `..._polar_pull_maps.png`, `..._m2_profiles.png` | polar-cell maps and m = 2 profiles of the S/f models |
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
