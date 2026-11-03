# Key results - Test 2 at a 12.5% distance error

The 10% campaign (2026-10-29 to 11-01) repeated at **12.5%** (sigma_DM = 0.125 * 5 / ln 10 = 0.2714 mag),
run 2026-11-02/03. Nothing about the method changed: the same 257-point focused Sobol design, the same
DR3-Mira-like selection function, the same chi2 schemes and the same `sf` library weighting. Only the
distance error of the mock is different.

Six searches, 1,542 orbit integrations. Four ran on the cluster (t = 2, 3, 4, 5.5) and two here
(t = 7 with and without the selection function). Full history in `PROGRESS.md`.

**Scheme B** (the dispersion ring averages dropped) and the `sf` weighting throughout - the combination
the 10% campaign settled on. Scheme A trades chi2 against library temperature and reads high; the
launcher emulates scheme A by default, so the scheme B posteriors here were produced separately from
the same `evals.csv` (no re-integration).

Truth is the simulation directory's **finite-difference** column, the convention adopted 2026-11-01.

---

## 1. Pattern speed: 10% against 12.5%  (`01_pattern_speed_10pct_vs_12p5pct.png`)

| snapshot | truth | 10% | bias | 12.5% | bias | change |
|---|---|---|---|---|---|---|
| t = 2 | 43.90 | 43.40 | -1.1% | 42.96 | -2.1% | -1.0% |
| t = 3 | 38.67 | 40.63 | +5.1% | 40.99 | +6.0% | +0.9% |
| t = 4 | 33.28 | 34.86 | +4.7% | 36.43 | +9.5% | +4.7% |
| t = 5.5 | 28.21 | 28.06 | -0.5% | 28.34 | +0.5% | +1.0% |
| t = 7 | 25.20 | 25.92 | +2.9% | 25.85 | +2.6% | -0.3% |
| **mean \|bias\|** | | | **2.9%** | | **4.1%** | |

- **The method degrades gracefully.** Mean |bias| 2.9% -> 4.1%, and the 68% interval widens by
  only 16% (mean half-width 1.49 -> 1.72 km/s/kpc). A 25% larger distance error costs
  about 1 point of bias.
- **Three of the five snapshots barely move.** t = 7 improves slightly (+2.9% -> +2.6%), t = 5.5 and t = 2
  shift by about 1 point. The mean is carried by one snapshot.
- **t = 4 is the outlier:** +4.7% -> +9.5%, which is 4.7 of the 6.0 points of total bias change; the other
  four average +0.4%. Checked and **not** a design-box artefact - its posterior (Omega_b 34.09-38.46,
  phi_b 17.5-30.9) sits well inside the box (23.71-49.18, 14.20-43.84). t = 4 was also the snapshot whose
  bar-angle calibration transferred worst in the 2026-11-01 work (slope -0.639 against t = 7's -0.669), so
  the same sensitivity may be behind both. **Not explained.**
- The quoted widths are the tempered ones, which the 10% campaign showed over-cover (100% against 68%);
  they are not statistical 1 sigma. See PROGRESS 2026-11-01.

## 2. Recovery at 12.5% alone  (`02_pattern_speed_recovery_12p5pct.png`)

The standard figure in the 10% campaign's format, all library variants, for direct comparison with
`figures/test2_pattern_speed_recovery_schemeB_fdtruth.png`.

## 3. Corner plots  (`corner_by_snapshot/`)

GP emulator posteriors in (Omega_b, phi_b, S, bar_frac); orange = truth.

| file | run |
|---|---|
| `t2_sf_schemeB_corner.png` | t = 2, selection function, scheme B, `sf` |
| `t3_sf_schemeB_corner.png` | t = 3, selection function, scheme B, `sf` |
| `t4_sf_schemeB_corner.png` | t = 4, selection function, scheme B, `sf` |
| `t5.5_sf_schemeB_corner.png` | t = 5.5, selection function, scheme B, `sf` |
| `t7_sf_schemeB_corner.png` | t = 7, selection function, scheme B, `sf` |
| `t7_NOSF_benchmark_schemeA_{observed,asis}_corner.png` | t = 7 no-selection benchmark (scheme A: it has no `sf` variant) |

The no-selection benchmark is the control that isolates what the selection function does - it draws
uniformly from the centrally concentrated R < 8 kpc particles (median R 1.88 kpc against 3.33), so it is
not a sixth snapshot.

## Reproducing

```
cd scripts
sh run_test2_12p5pct_local.sh        # t = 7 with and without the selection function
sh run_test2_12p5pct_cluster.sh      # t = 2, 3, 4, 5.5
python plot_test2_error_level_comparison.py                        # figure 1
python plot_test2_pattern_speed_recovery.py --runs <the six 12.5% runs> --tag _12p5pct_schemeB
```

Scheme B posteriors (the launcher does scheme A):

```
python test1_gp_emulator.py --run <run> \
  --column logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_sf --loo-top 32
```

Trajectories are archived on the portable disk; `runs/<run>/` here keeps only the per-point `.npz`,
which is all the emulator needs.
