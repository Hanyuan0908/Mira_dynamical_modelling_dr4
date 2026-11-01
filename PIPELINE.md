# Pipeline guide: mock data → orbit integration → likelihood → posterior

This is a working guide to the repository, written for running the pipeline elsewhere (e.g. by an agent on a
computing cluster). Results and their history are in [PROGRESS.md](PROGRESS.md). The overview and the Au18 set-up
are in [README.md](README.md).

## 1. What the method does

The method is a steady-state test of a barred potential, using a 6D stellar sample.

1. Take the observed stars (positions and velocities).
2. Rotate them by a trial bar angle φ_b, so that the assumed bar lies on the x-axis.
3. Integrate every star for 1 Gyr in a trial potential rotating at a trial pattern speed Ω_b. The trial potential has
   the real galaxy's axisymmetric part, plus its non-axisymmetric part rescaled in strength (`bar_frac`) and
   length (`S`).
4. Time-average the orbits (t ≥ 0.25 Gyr) into a "library".
5. Compare the library's kinematic maps with the data's.

Only in the right potential does the sample stay stationary in the bar frame, so the library reproduces the data.
The four parameters are (|Ω_b|, φ_b, S, bar_frac). They are explored with a scrambled Sobol design; a Gaussian-process
emulator of the log-likelihood then gives the posterior.

**Current likelihood** ("proposed scheme"):
- **Fourier fit.** In each block of 0.3 kpc rings (R < 6 kpc) and |z| layers (0–1, 1–2, 2–3 kpc), fit a Fourier
  series directly to the individual stars:
  - ⟨v_R⟩ with m = 2, 4;
  - ⟨v_φ⟩, σ_R² and σ_φ² with m = 0, 2, 4.
- **Data errors.** The joint covariance of each block's coefficients, from 300 bootstrap resamplings of the stars.
- **χ² and tempering.** χ² over all blocks, tempered as logL = −χ² / (2√(2 N_data)).

**Ground truth of the test galaxy:** |Ω_b| = 26.02 km/s/kpc, φ_b = 27°, S = 1, bar_frac = 1.

**Status (2026-10-29):**
- Test 1 (perfect data) recovers the truth.
- Test 2b (distance errors): the method needs distances good to about 2%. At 7.5% the posterior shifts by about +2.5
  km/s/kpc in Ω_b and +40% in bar_frac.
- Test 2a (selection function): the mocks exist — 10% distance errors with the DR3-Mira-like selection at model
  times 2, 3, 4, 5.5 and 7 (Ω_b 44.3 down to 26.0 km/s/kpc), plus a no-selection benchmark at t = 7. The
  **library-side selection weighting is still to be written** (section 5.2), as are the potentials for the four
  new snapshots.

## 2. Environment

- Conda environment `astro312`. Versions used:

  | package | version |
  |---|---|
  | Python | 3.12.13 |
  | AGAMA | 1.0.156 |
  | numpy | 2.4.3 |
  | scipy | 1.17.1 |
  | scikit-learn | 1.8.0 |
  | emcee | 3.1.6 |
  | corner | 2.2.3 |
  | h5py | 3.16.0 |
  | pandas | 2.3.3 |
  | matplotlib | 3.10.8 |

- `scikit-learn` is needed only by `test1_gp_emulator.py`, `healpy` only by `make_mock_selection.py`. Both are
  optional for the integration runs.
- Older stacks work: numpy 1.26 / scipy 1.15 / AGAMA 1.0 (May 2024) / sklearn 1.8 / healpy 1.17 runs the pipeline.
  Pin `healpy` to a build matching your numpy's C-API (1.20 needs numpy 2.x).
- AGAMA integrates orbits with OpenMP. `common.py` defaults `OMP_NUM_THREADS` to `os.cpu_count()`; set it
  explicitly for a shared machine or a batch job and that wins. Everything else is single-threaded numpy.
  Measured scaling on an 8-core/16-thread i7-10700KF (10k orbits): 4 threads 178 orbits/s, 8 → 356, 12 → 391,
  16 → 426, so hyperthreading is worth about 20%. chi2 is bit-identical at 8 and 16 threads.
- Run scripts from `scripts/`: they import each other as plain modules (`import common`, `import observe`, …).
- Units throughout: Msun, kpc, km/s. AGAMA's time unit is 1 kpc/(km/s) = 0.9778 Gyr.
- Reference timing (Apple M5 Pro, 15 cores): one trial potential with 50k orbits, 1 Gyr, 41 snapshots, accuracy
  10⁻⁵ takes **80–95 s**; 500k orbits take 650–750 s. On an i7-10700KF (8 cores / 16 threads) the same 50k
  integration takes **~150 s**, and a whole design point (integration + data fit + library fit + χ²) 2.6 min,
  so a 257-point search is about 11 h.
  - Concurrent heavy jobs slow it 2–3×.
  - Peak memory per 50k-orbit integration has not been measured. The trajectory alone is 47 MB (float32) and about
    100 MB in AGAMA's float64 output.

## 3. Files that are not in git (copy them, or rebuild them)

`data/` and `mock/` are git-ignored. For the N-body pipeline you need:

| file | size | made by | needed for |
|---|---|---|---|
| `mock/Nbody/t_t0_7` | 126 MB | the simulation (NEMO format) | rebuilding everything below |
| `mock/nbody_stars_mwframe.hdf5` | 86 MB | `prepare_nbody_mwframe.py` | mocks; the 1.49 M-particle ground truth |
| `mock/nbody_dm_mwframe.hdf5` | 43 MB | `prepare_nbody_mwframe.py` | `build_nbody_potential.py` |
| `mock/nbody_potential/*.ini` | 2 MB | `build_nbody_potential.py`, `potential_family.py` | every trial potential |
| `mock/nbody_mock_perfect_R8_N50k.hdf5` | 4 MB | `make_mock_perfect.py --galaxy nbody` | the data (Test 1) and the base of all mocks |
| `mock/nbody_mock_distance_error_sigmaDM*.hdf5` | 8 MB each | `make_mock_distance_error.py` | Test 2b data |
| `data/gaia_dr3_lpv_mira.csv` | 87 MB | external | only the default σ_DM (median DM_err) |
| `results/test1_nbody_qmc_focused/design.csv` | small, in git | `test1_nbody_qmc_focused_run.py` | the 257 trial points (Sobol) |

`mock/halo_18/` (4 GB, Auriga Au18) is not needed for the N-body pipeline.

Big outputs go to `common.RUNS_DIR` (`ROOT/runs`, or `$MIRA_DR4_RUNS`); see section 7. Allow
**~12.6 GB per 257-point run** (one 47 MB trajectory per point).

## 4. Frames and conventions (get these right)

- **MW-like frame:** astropy's default Galactocentric frame (`common.py`).
  - Sun at x = −R0, with R0 = 8.122 kpc, z_sun = 0.0208 kpc and v_sun = (12.9, 245.6, 7.78) km/s.
  - The disc rotates clockwise seen from +z (L_z < 0).
  - The bar lies at −27° from the Sun–Galactic-centre line.
- **Bar frame:** `steady_state_likelihood.rotate_z(stars_mw, +φ_b)` puts the assumed bar on the x-axis.
  - Orbits are integrated in the bar frame with `agama.orbit(..., Omega=-|Ω_b|)`; negative means clockwise.
- **Heliocentric observables:** `observe.to_galactic` / `observe.to_galactocentric` (AGAMA): l, b [deg], d [kpc],
  pm_l*, pm_b [mas/yr], v_los [km/s].
- **Distance errors:** Gaussian in distance modulus, DM_obs = DM_true + σ_DM N(0,1)
  (`observe.scatter_distance_modulus`); the fractional error is ln(10)/5 · σ_DM.
  - Angles, proper motions and v_los are kept exact.
  - Positions and tangential velocities are recomputed from d_obs.
- **Trial potential:** `steady_state_likelihood.bar_aligned_potential(bar_frac, S, "nbody")`.
  - Fixed parts: the spherical DM halo plus the axisymmetric baryons.
  - Non-axisymmetric part: `bar_frac · [Φ_bary(x/S) − Φ_bary,axi(x/S)]`, merged into one mmax-8 CylSpline.
  - AGAMA `scale=[bar_frac·S, S]` multiplies by bar_frac and stretches by S (checked to 6×10⁻⁶).
- **Library:**
  - Snapshots with t ≥ 0.25 Gyr (`T_MIN_GYR`); 41 snapshots over 1 Gyr, every 25 Myr.
  - Each orbit is weighted by 1/τ, where τ is the fraction of its time inside the selection volume. That volume is
    currently R < 8 kpc (`R_SEL`).
- **Folding:** |z| and v_z·sign(z). Only even m is used, so φ needs no folding.

## 5. Pipeline stages

### 5.0 One-off galaxy set-up (already done; rebuild only if the inputs change)

```
python prepare_nbody_mwframe.py        # mock/Nbody/t_t0_7 -> mock/nbody_{stars,dm}_mwframe.hdf5 (centre, spin, bar at -27 deg)
python build_nbody_potential.py        # -> mock/nbody_potential/{dm_spherical,bulge,disc,...}.ini  (AGAMA expansions)
python potential_family.py --galaxy nbody   # -> mock/nbody_potential/Nbody_{full,axi,baryon_full,baryon_axi}.ini (default --galaxy is au18)
python measure_au18_pattern_speed.py --galaxy nbody   # -> results/nbody_ground_truth.json (Omega_p = 26.02)
```

`measure_au18_pattern_speed.py` imports external code from `external/patternSpeed-main` (`$PATTERN_SPEED_DIR`; see section 7).

### 5.1 Mock data

```
python make_mock_perfect.py --galaxy nbody                 # 50k random star particles with R < 8 kpc, seed 2026
python make_mock_distance_error.py                         # sigma_DM = DR3 Mira median (0.2654 mag, 12.2%), seed 2029
python make_mock_distance_error.py --frac 0.10             # 10% (sigma_DM = 0.2171 mag); --sigma-dm also works
```

**Snapshots as test galaxies.** The N-body run is one secularly evolving galaxy, so each snapshot is a different
pattern speed with the same physics (`measure_nbody_pattern_speed_evolution.py`). Any script taking `--galaxy`
accepts `nbody_t<T>` for the snapshot at model time T, and `prepare_nbody_mwframe.py --time T` builds it:

```
python prepare_nbody_mwframe.py --time 5.5                 # -> mock/nbody_t5.5_{stars,dm}_mwframe.hdf5
python make_mock_selection.py --galaxy nbody_t5.5 --frac 0.10
```

`nbody` is model t = 7 and keeps the original un-tagged file names; `nbody_t<T>` tags every file with its time.
Model time is not Gyr: 1 unit = 0.9778 Gyr, so t = 7 is 6.84 Gyr.

- **`make_mock_perfect.py`:**
  - Draws N stars without replacement from the particles with R < R_MAX.
  - Writes their true phase space (x..vz) and the error-free observables (l, b, d, pml, pmb, vlos).
  - Records `star_index` (index into the full particle file) and `pool_size`.
- **`make_mock_distance_error.py`:**
  - Reads the perfect mock and adds the distance-modulus error.
  - Writes the observed x..vz and d under the usual names, so every reader works unchanged; the true values go to
    `<name>_true`, plus `DM_true`, `DM_obs`, `DM_err`.
  - The same seed gives the same normal deviates for every σ_DM. Levels are therefore nested, e.g. the 7.5% mock is
    the 12.2% mock scaled.
  - Caveat: the stars are selected on their **true** R < 8 kpc (not realistic). 1.7–2.7% end up at observed R > 8.

### 5.2 Selection function (Test 2) — mock implemented, library weighting not yet

```
python make_mock_selection.py --frac 0.10                      # 10% (sigma_DM 0.2171), 50k expected stars
python make_mock_selection.py --galaxy nbody_t4 --frac 0.10    # another snapshot
```

- **Mock (`make_mock_selection.py`).** Recipe agreed 2026-10-29:
  - Observe all 2 M N-body star particles with the distance-modulus error and keep those with observed R < 8 kpc.
  - Select them with S(l, b, d) = S_sky[HEALPix nside 32] × S_dist[0.25 kpc in d] × S_R[0.25 kpc in Galactocentric R].
  - S is fitted by iterative proportional fitting so the mock reproduces the DR3 Miras' (DM_err < 0.4, R < 8 kpc)
    sky map, distance histogram and R histogram. Each particle is kept with probability min(S, 1).
  - The S tables (`S_sky`, `dist_edges`/`S_dist`, `R_edges`/`S_R`) are stored in the mock file.
  - Same output fields as the other mocks.
- **Still to do: the library.** Earlier plan, written before the mock existed:

Plan agreed earlier (PROGRESS.md, "Next steps"):
- **Data.** Draw the mock through a Gaia-like selection S(l, b, m) (footprint, extinction, magnitude limit),
  applied to the *observed* quantities.
- **Library.**
  - Weight each orbit by 1/τ_i, with τ_i = ⟨S(x(t))⟩ along the orbit (time average).
  - Weight each library sample by S(x_j).
  - This generalises the current R < 8 kpc indicator (`R_SEL`, used in
    `test1_nbody_qmc_focused_run.library_samples` and `test1_nbody_qmc_focused_fourier_stars.library_fit`).
- **Folding.** S breaks the folding symmetries, so the data need unfolded cells (the library can stay folded).
- **Where it plugs in:**
  - a new `make_mock_selection.py` (or an option of the mock scripts);
  - a selection callable passed to the library weighting in place of `R < ssl.R_SEL`.

### 5.3 Trial design (Sobol)

- `test1_nbody_qmc.py`: the truth plus 64 scrambled Sobol points (seed 2027) in a wide box:
  - |Ω_b| 18–38, φ_b 14–44°, S 0.75–1.45, bar_frac 0.55–1.6;
  - S and bar_frac log-uniform.
- `test1_nbody_qmc_focused_run.py`: those 65 plus 192 new Sobol points (seed 2028) in a box focused on the
  perfect-data posterior:
  - |Ω_b| 20–32, φ_b 15–39°, S 0.6–1.65, bar_frac 0.45–2.0;
  - 257 points in total, in `results/test1_nbody_qmc_focused/design.csv` (columns idx, kind, batch, Omega_abs,
    phi_b_deg, S, bar_frac, traj_file).
  - Caveat: with 7.5% distance errors the posterior reaches the box edges at high Ω_b and bar_frac. For noisy data,
    widen the box.

### 5.4 Orbit integration (the expensive step) and the likelihood

The core is three calls (see `test2b_qmc_distance_error_run.py`):

```python
pot = ssl.bar_aligned_potential(bar_frac, S, "nbody")
_, traj = ssl.integrate(pot, ssl.rotate_z(stars_mw, np.radians(phi_b)), -abs(Omega_abs))   # (N, 41, 6) float32, bar frame
cd, cov, n, ok = fs.data_fit(stars_mw, phi_b, boot_w)          # data coefficients + joint bootstrap covariance
cl = fs.library_fit(traj, times, ok)                           # library coefficients (t >= 0.25 Gyr, 1/tau weights)
chi2, ndata = fs.scheme_chi2(cd, cl, cov, valid, fs.SCHEMES[...])
```

`ssl.ACCURACY` must be set to 1e-5 (the global searches use `old.ACCURACY`).

Runners (each loops over all design points serially; resumable, skipping points whose output file exists):

| script | data | per point | outputs |
|---|---|---|---|
| `test1_nbody_qmc.py` | perfect mock | integrate + KDE / local-Gaussian / moment likelihoods (heavier) | `~/mira_dr4_runs/test1_nbody/{NNN.npz, traj_NNN.npy}`, `results/test1_nbody/qmc_evals.csv` |
| `test1_nbody_qmc_focused_run.py` | perfect mock | integrate new points; library sums on a fine grid | `~/mira_dr4_runs/test1_nbody_qmc_focused/{lib_NNN.npz, traj_NNN.npy}` |
| `test1_nbody_qmc_focused_fourier_stars.py` | perfect mock | per-star Fourier fit of data (300 bootstraps) and library, χ² for 4 schemes | `fourier_stars_NNN.npz`, merged into `results/test1_nbody_qmc_focused/evals.csv` |
| `test2b_qmc_distance_error_run.py` | a distance-error mock (`MOCK` constant, currently the 7.5% file) | integrate the observed stars; library as is and observed again with the same error; χ² | `~/mira_dr4_runs/test2b_qmc_7p5/{traj_NNN.npy, fourier_NNN.npz}`, `results/test2b_qmc_7p5/evals.csv` |

- **Per-point cost:** the integration (80–95 s) plus the data bootstrap at that point's angle (about 15 s).
- **The data bootstrap could be done once.** A different φ_b only rotates each mode's (cos, sin) pair by m Δφ_b, so
  c' = T c and C' = T C Tᵀ exactly.
  - Implemented in `test2b_qmc_7p5_rescore_coadded.rotation`.
  - Checked against per-angle fits to 2×10⁻¹⁰.

Post-processing:
- `test2b_qmc_7p5_rescore_coadded.py`: rescores a finished search with the co-added error (user's recipe for real
  data). The error is the covariance over 300 re-draws of every observed star's distance, plus the bootstrap
  covariance of the unsmeared data, rotated to each angle. Writes `results/test2b_qmc_7p5/evals_coadded.csv`.
- `run_test1_nbody_qmc_focused_postprocess.sh`: the Test 1 post-processing chain.

### 5.5 Posterior and figures

```
python test1_gp_emulator.py --run <run> --column <logLs column> --loo-top 32
python plot_test1_logl_corner.py --run <run> --column <logLs column> --levels=-100,-50,-20,-10,-4.5,-2,-0.5 --lines=-4.5,-2
```

- **Runs** are registered in `test1_gp_emulator.RUNS`: `nbody_qmc_focused`, `test2b_qmc_7p5`, `test2b_qmc_7p5_coadded`.
  Each entry gives the CSV, box, truth and output prefix; add an entry for each new search.
- **Emulator:**
  - quadratic mean, plus a GP with a Matern-5/2 kernel (one length scale per parameter) and white noise, trained on
    the `kind == "qmc"` rows;
  - posterior from emcee, seeded (32 walkers, 3000 steps, 1000 burn-in);
  - outputs `results/<run>/emulator_results_<column>.json` and `figures/<prefix>_emulator_{corner,slices,loo}_<column>.png`.
- **Main columns:**
  - perfect data: `logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024`;
  - 7.5% search: `logLs_observed` / `logLs_asis` (co-added error) and `..._observed` / `..._asis` (bootstrap error).

## 6. Key modules

| module | role |
|---|---|
| `common.py` | paths (`ROOT`, `MOCK_DIR`, `FIG_DIR`, `MIRA_DR3_FILE`), solar parameters, bar angle, plot style |
| `galaxies.py` | registry of the test galaxies (`au18`, `nbody`): star files, potential directories, truth files |
| `nbody.py` | reader for the NEMO snapshot `mock/Nbody/t_t0_7` (units: G = 1, mass unit 232,508 Msun) |
| `observe.py` | Galactocentric ↔ Galactic conversions, distance-error scattering |
| `potential_family.py` | builds the trial-potential component files; `files("nbody")` |
| `steady_state_likelihood.py` (`ssl`) | `bar_aligned_potential`, `rotate_z`, `integrate`, `cell_moments`, constants (`T_GYR`, `N_SNAP`, `T_MIN_GYR`, `R_SEL`, `ACCURACY`) |
| `test1_nbody_qmc_focused_fourier_stars.py` (`fs`) | the current likelihood: grid (`RING_EDGES`, `LAYER_EDGES`), `prepare`, two-stage per-star least-squares `fit`, `data_fit` (bootstrap), `library_fit`, `scheme_chi2`, `SCHEMES` |
| `test1_gp_emulator.py` | emulator, run registry, posterior |
| `plot_distance_error_effect_on_fourier_coefficients.py` | `observe_library_samples` / `forward_smeared_library_fit`: observe a library with distance errors |

Checks worth re-running after changes:
- `check_fourier_star_fit_vs_true_maps.py`: the Fourier fit and its errors against the 1.49 M-particle truth.
- `check_bootstrap_errors_with_distance_errors.py`: bootstrap errors with distance errors.

## 7. Running on another machine

The scripts read and write through repo-relative paths, so a fresh checkout needs only the data in place
(section 3) and two optional environment variables.

1. **Inputs.** `data/` and `mock/` are git-ignored. Either copy the files in, or symlink them:

   ```
   ln -s /path/to/gaia_dr3_lpv_mira.csv        data/gaia_dr3_lpv_mira.csv
   ln -s /path/to/Bar_model_TG21/model         mock/Nbody
   ln -s /path/to/patternSpeed-main            external/patternSpeed-main
   ```

   `external/` holds third-party code (the Dehnen, Semczuk & Schönrich pattern-speed package) and is
   git-ignored too. Override its location with `$PATTERN_SPEED_DIR`.

2. **Bulky outputs** (one 47 MB trajectory per design point; **~12.6 GB per 257-point run**) go to
   `common.RUNS_DIR`, which is `ROOT/runs` by default. Point it at scratch space with either a symlink
   (`ln -s /scratch/you/mira_dr4_runs runs`) or the environment variable:

   ```
   export MIRA_DR4_RUNS=/scratch/you/mira_dr4_runs
   ```

   Used by `test1_nbody_qmc.py`, `test1_nbody_qmc_focused_run.py`, `test2b_qmc_distance_error_run.py`,
   `test2b_smeared_initial_conditions.py` and `test2b_distance_error_tolerance.py`.

3. **Trajectory paths in `design.csv`** are stored relative to `RUNS_DIR` (`test1_nbody/traj_007.npy`).
   `test1_nbody_qmc_focused_run.traj_path` resolves them, and still accepts the absolute paths written by
   older runs when those files exist, so an existing design keeps working.

4. **Snapshot times** come from `ssl.snapshot_times()` (`np.linspace(0, T_GYR, N_SNAP)`), not from a previous
   run's output file, so nothing has to exist before the first point is integrated.

5. **AGAMA version.** `agama.orbit`'s `verbose` argument was added after 1.0;
   `potential_family.ORBIT_KW` probes for it, so the scripts run on either build.

6. **Threads and seeds.**
   - Set `OMP_NUM_THREADS` per task. AGAMA integrates orbits with OpenMP over all cores; everything else
     is single-threaded numpy.
   - Seeds are fixed in the scripts: mock 2026, distance error 2029, Sobol 2027/2028, bootstrap 0, library
     observing 11. Re-running the fits on saved trajectories reproduces the saved coefficients exactly.
   - Whether AGAMA's integrations are bit-identical across machines has not been tested.

7. **Parallelism.** The runners process the whole design serially and append to one `evals.csv`. For job arrays,
   give each task a subset of points (add an index-range option) and write only per-point files (`fourier_NNN.npz`).
   Then merge into `evals.csv` in one final step. The per-point `.npz` already contain every χ²/N_data needed.

8. **Mock choice.** `test2b_qmc_distance_error_run.py` has the mock file and output directory as constants (`MOCK`,
   `OUT_DIR`, `RES_DIR`). Make them arguments to run other error levels or selection-function mocks.

## 8. Known issues to keep in mind

- **Joint bootstrap covariance in sparse blocks.** It overweights blocks with ≲ 500 stars (|z| ≥ 1 kpc): χ²/N 1.25
  instead of 1 at |z| 1–2. Options: bigger blocks there, or a model-based covariance.
- **Hot library.** The truth library is ~1.5 km/s hotter than the data in the ring-averaged dispersions; cause
  unknown. It is the main model misfit at the truth on perfect data (χ²/N ≈ 0.16–0.18 of the 1.50).
- **Distance errors break the "orbits from the observed stars" idea.**
  - Tolerance is about 2% in distance.
  - Re-weighting candidate distances per star overfits any potential (χ²/N ~ 10⁻⁴ everywhere).
  - Posterior draws with the true N-body as a prior still leave χ²/N ≈ 4.
  - Whether the 7.5% pattern-speed shift is systematic needs more mock realisations: different stars, error draws
    and galaxies.
