# paper_figures/

The results worth publishing, copied from across `figures/`. Everything here is a copy — the
originals stay where their scripts write them, so nothing breaks if this folder is reorganised.

Rebuild with `python scripts/collect_paper_figures.py`; add one with
`--add <figure> [--as <name>] [--section <dir>]`.

**Standard throughout** (adopted 2026-11-03/04): orbit trajectories saved at a **50 Myr** cadence;
the library observed with **K = 10** pooled distance-error draws per sample; the selection function
**estimated from the data** as S = rho_obs/rho_lib; chi2 **scheme B**; truth = the simulation
table's **finite-difference** column. No figure here uses the mock's own S tables — no survey has
them. See `emulator_validation/` for why K = 10 and what the emulator's errors actually are.

---

## test2_12p5pct/ — Test 2: 12.5% distance error + a DR3-Mira-like selection function

**Twelve snapshots** of one N-body bar, Omega_b 25-44 km/s/kpc, 257-point focused Sobol design each.

| file | what it shows |
|---|---|
| `01_pattern_speed_recovery_12_snapshots.png` | **the main result**: recovered against true Omega_b, and the fractional bias |
| `02_pattern_speed_10pct_vs_12p5pct.png` | how the recovery changes from a 10% to a 12.5% distance error (five matched snapshots) |
| `03_mock_vs_dr3_miras.png` | the t = 7 mock against the DR3 Miras: sky, distance, R, face-on |
| `04_all_mocks_vs_dr3_miras.png` | all twelve mocks face-on, with their agreement with the DR3 margins |
| `05_snapshot_coverage.png` | where the twelve sit on the run's pattern-speed track |
| `corner_sfratio/` | posterior per snapshot in (Omega_b, phi_b, S, bar_frac), named by snapshot |

### The twelve points

| snapshot | Omega_b true | A2 | recovered | bias | LOO | truth at |
|---|---|---|---|---|---|---|
| t2 | 43.90 | 0.469 | 42.11 [40.03–44.46] | -4.1% | 1.22 | +0.97σ |
| t2.375 | 40.87 | 0.416 | 41.70 [39.41–44.31] | +2.0% | 0.68 | -0.09σ |
| t2.625 | 40.21 | 0.442 | 40.34 [38.19–42.76] | +0.3% | 1.05 | -0.10σ |
| t3 | 38.67 | 0.452 | 39.11 [37.08–41.37] | +1.1% | 1.04 | +0.19σ |
| t3.375 | 36.83 | 0.459 | 37.76 [35.95–39.85] | +2.5% | 0.91 | -0.08σ |
| t3.625 | 34.79 | 0.450 | 35.97 [33.98–38.28] | +3.4% | 1.01 | +0.02σ |
| t4 | 33.28 | 0.479 | 35.18 [32.53–37.41] | +5.7% | 0.87 | -0.50σ |
| t4.5 | 32.07 | 0.470 | 32.98 [30.66–35.10] | +2.8% | 1.30 | -0.26σ |
| t5 | 30.60 | 0.461 | 29.23 [27.69–31.79] | -4.5% | 1.10 | +0.70σ |
| t5.5 | 28.21 | 0.450 | 26.77 [25.73–28.01] | -5.1% | 1.20 | +1.82σ |
| t6.25 | 26.58 | 0.433 | 26.34 [25.16–27.75] | -0.9% | 1.21 | +0.67σ |
| t7 | 25.20 | 0.407 | 24.83 [24.00–25.73] | -1.5% | 1.06 | +1.37σ |

**mean |bias| 2.8%, mean bias +0.2%, rms 3.3%.**

The method is **unbiased on average** but the scatter is **structured, not random**: above
Omega_b ~ 31 the mean bias is **+1.7%** (n = 8), below it **-3.0%** (n = 4), with a sharp
transition between Omega_b = 32.07 (+2.8%) and 30.60 (-4.5%). The original five snapshots read
-4.1, +1.1, +5.7, -5.1, -1.5 and looked like noise; the structure only appears once the gaps are
filled. Where exactly the step sits is not pinned — only that it lies between those two points.

Caveats that belong in the text:
- the twelve snapshots are **one galaxy in a secular sequence**, so Omega_b, A2 and R_bar co-vary;
  this is a systematic *along this sequence*, not a Milky Way systematic;
- the finite-difference truth column carries **0.66 km/s/kpc (2.1%)** of its own scatter, so
  individual points have that much slop (the 7.3-point step is far outside it);
- phi_b is recovered poorly by the orbit library throughout; the bar angle is better taken from
  the v_R quadrupole (PROGRESS 2026-11-01).

## emulator_validation/ — is the GP emulator trustworthy?

| file | what it shows |
|---|---|
| `01_leave_one_out.png` | 1,536 held-out predictions: error against position in the design, and whether the emulator's own error bars are honest |
| `02_pooled_draws_K10.png` | why K = 10: the library noise falls as 1/sqrt(K), chi2 falls by a parameter-dependent amount, and the emulator improves 2.3x |
| `loo_per_run/` | the per-run leave-one-out scatter, named by snapshot |

Headline: near the likelihood peak — the top 32 points per run, carrying **94% of the posterior
weight** — the held-out rms is **1.34** in tempered logL and the pull sd is **0.76**, i.e. the
emulator's error bars are *conservative* exactly where it matters. The large errors visible in the
left panel are in the tails, which carry 0.0% of the weight.

## test1_perfect_mock/ — Test 1b: the control, no errors and no selection

| file | what it shows |
|---|---|
| `01_posterior_corner_schemeA.png` | posterior with all Fourier terms |
| `02_posterior_corner_schemeB.png` | the same under scheme B, the scheme the Test 2 campaigns use |
| `03_fourier_coefficients_data_vs_models_absz_0to1.png` | the actual fit: m = 0, cos/sin 2phi, cos/sin 4phi of the four moments against R |
| `03_..._absz_1to2.png` | the same for 1 < |z| < 2 kpc |
