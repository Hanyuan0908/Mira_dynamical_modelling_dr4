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
- Test 2a (selection function) is next and is **not implemented yet** (section 5.2).

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

- AGAMA integrates orbits with OpenMP over all cores; set `OMP_NUM_THREADS` to the cores of the job. Everything else
  is single-threaded numpy.
- Run scripts from `scripts/`: they import each other as plain modules (`import common`, `import observe`, …).
- Units throughout: Msun, kpc, km/s. AGAMA's time unit is 1 kpc/(km/s) = 0.9778 Gyr.
- Reference timing (Apple M5 Pro, 15 cores): one trial potential with 50k orbits, 1 Gyr, 41 snapshots, accuracy
  10⁻⁵ takes **80–95 s**; 500k orbits take 650–750 s.
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

Big outputs are written outside the repo, under `~/mira_dr4_runs/<run>/` (`Path.home()` in the scripts). Allow
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

`measure_au18_pattern_speed.py` imports external code from a hard-coded path (`PATTERN_SPEED_DIR`; see section 7).

### 5.1 Mock data

```
python make_mock_perfect.py --galaxy nbody                 # 50k random star particles with R < 8 kpc, seed 2026
python make_mock_distance_error.py                         # sigma_DM = DR3 Mira median (0.2654 mag, 12.2%), seed 2029
python make_mock_distance_error.py --sigma-dm 0.16286      # 7.5% (sigma_DM = 0.075 / (ln 10 / 5))
```

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

### 5.2 Selection function (Test 2a) — NOT IMPLEMENTED yet

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

## 7. Running on a cluster: what to change first

The scripts were written for one workstation. Before launching:

1. **Paths.**
   - `~/mira_dr4_runs` is hard-coded through `Path.home()` in `test1_nbody_qmc.py` (`ARR_DIR`),
     `test1_nbody_qmc_focused_run.py`, `test2b_qmc_distance_error_run.py`, `test2b_smeared_initial_conditions.py`
     and `test2b_distance_error_tolerance.py`. Point it at scratch space (e.g. through an environment variable).
   - `results/test1_nbody_qmc_focused/design.csv` stores **absolute** trajectory paths of this machine in
     `traj_file`, read by `test1_nbody_qmc_focused_run.py` and `..._fourier_stars.py`.
   - Absolute paths: `PATTERN_SPEED_DIR` in `measure_au18_pattern_speed.py` and `PY=` in
     `run_test1_nbody_qmc_focused_postprocess.sh`.
2. **Snapshot times.** Eleven scripts read them from an old run file (`np.load(old.ARR_DIR / "000.npz")["times"]`),
   including both runners, the likelihood module and every `test2b_*` script.
   On a fresh machine use `ssl.integrate`'s returned times, or `np.linspace(0, ssl.T_GYR, ssl.N_SNAP)`. They are the
   same.
3. **Parallelism.** The runners process the whole design serially and append to one `evals.csv`. For job arrays,
   give each task a subset of points (add an index-range option) and write only per-point files (`fourier_NNN.npz`).
   Then merge into `evals.csv` in one final step. The per-point `.npz` already contain every χ²/N_data needed.
4. **Mock choice.** `test2b_qmc_distance_error_run.py` has the mock file and output directory as constants (`MOCK`,
   `OUT_DIR`, `RES_DIR`). Make them arguments to run other error levels or selection-function mocks.
5. **Threads and seeds.**
   - Set `OMP_NUM_THREADS` per task.
   - Seeds are fixed in the scripts: mock 2026, distance error 2029, Sobol 2027/2028, bootstrap 0, library
     observing 11. Re-running the fits on saved trajectories reproduces the saved coefficients exactly.
   - Whether AGAMA's integrations are bit-identical across machines has not been tested.
6. **Not needed on Linux:** `caffeinate` (macOS; used only in ad-hoc commands, not in the scripts).

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
