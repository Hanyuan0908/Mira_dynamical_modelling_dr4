# Scheme B posteriors — how the selection function is obtained

Test 2, 10% distance error, five snapshots of the N-body bar, each with a DR3-Mira-like selection
function. χ² is **scheme B** (`fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024`), which drops
the σ_R², σ_φ² ring averages.

The question these answer: **can S be estimated from the data instead of being known?**

| tag | how S is obtained | available for real data? |
|---|---|---|
| `S-none` | no selection weighting at all | — (control) |
| `S-ratio` | ρ_obs/ρ_lib per (HEALPix nside 16 × 1 kpc distance) cell | **yes** |
| `S-separable` | S = s_sky × s_dist × s_R fitted to the data's three margins | **yes** |
| `S-true` | the exact tables that generated the mock | no — upper bound |

> **Superseded (2026-11-03).** The table below is scored against *our own measured* Omega_b. The
> project adopted the tabulated **finite-difference** column as the truth on 2026-11-01, and the same
> posteriors scored that way reverse the ordering: `S-ratio` 2.0%, `S-separable` 2.2%, `S-none` 5.3%,
> `S-true` **2.9%** - the estimated S *beats* the exact S, rather than recovering 70% of its gain.
> See PROGRESS 2026-11-02/03.

## Summary

| how S is obtained | mean \|Ω_b bias\| | mean worst offset | mean LOO |
|---|---|---|---|
| `S-none` | 7.3% | 2.6σ | 1.34 |
| `S-ratio` | 3.6% | 1.4σ | 1.05 |
| `S-separable` | 4.0% | 1.4σ | 0.94 |
| `S-true` | 2.1% | 1.0σ | 0.92 |

**An estimated S recovers about 70% of the gain from knowing S exactly**, and the direct density
ratio is marginally better than the separable fit (4 of 5 snapshots). Neither needs anything the
data do not already contain.

## Per run

| file | truth Ω_b | recovered | bias | φ_b | S | bar_frac | worst | LOO |
|---|---|---|---|---|---|---|---|---|
| `t2__S-none.png` | 44.25 | 40.93 | -7.5% | 21.7 | 0.930 | 1.388 | 1.8σ | 1.30 |
| `t2__S-ratio.png` | 44.25 | 42.55 | -3.8% | 18.6 | 1.033 | 0.811 | 2.5σ | 0.96 |
| `t2__S-separable.png` | 44.25 | 42.14 | -4.8% | 18.8 | 1.014 | 0.815 | 2.3σ | 0.83 |
| `t2__S-true.png` | 44.25 | 43.40 | -1.9% | 18.8 | 1.009 | 0.817 | 2.4σ | 0.82 |
| `t3__S-none.png` | 39.53 | 38.14 | -3.5% | 27.6 | 0.849 | 1.695 | 3.5σ | 1.76 |
| `t3__S-ratio.png` | 39.53 | 39.03 | -1.3% | 24.5 | 0.914 | 0.926 | 0.6σ | 1.15 |
| `t3__S-separable.png` | 39.53 | 38.90 | -1.6% | 25.0 | 0.914 | 0.965 | 0.5σ | 0.98 |
| `t3__S-true.png` | 39.53 | 40.63 | +2.8% | 24.5 | 0.953 | 0.939 | 0.6σ | 0.86 |
| `t4__S-none.png` | 33.96 | 30.24 | -10.9% | 19.3 | 0.896 | 1.216 | 2.8σ | 1.19 |
| `t4__S-ratio.png` | 33.96 | 33.34 | -1.8% | 22.9 | 1.039 | 0.858 | 0.7σ | 0.86 |
| `t4__S-separable.png` | 33.96 | 33.15 | -2.4% | 22.6 | 1.063 | 0.850 | 0.8σ | 0.85 |
| `t4__S-true.png` | 33.96 | 34.86 | +2.7% | 23.8 | 1.104 | 0.853 | 0.6σ | 0.78 |
| `t5.5__S-none.png` | 28.84 | 26.32 | -8.7% | 22.4 | 0.911 | 1.347 | 3.1σ | 1.37 |
| `t5.5__S-ratio.png` | 28.84 | 26.88 | -6.8% | 20.4 | 0.961 | 0.861 | 2.1σ | 1.17 |
| `t5.5__S-separable.png` | 28.84 | 26.87 | -6.8% | 19.8 | 0.959 | 0.864 | 2.2σ | 0.98 |
| `t5.5__S-true.png` | 28.84 | 28.06 | -2.7% | 21.3 | 0.991 | 0.897 | 0.9σ | 1.08 |
| `t7__S-none.png` | 26.02 | 24.54 | -5.7% | 26.3 | 0.913 | 1.486 | 1.8σ | 1.09 |
| `t7__S-ratio.png` | 26.02 | 24.86 | -4.4% | 22.6 | 0.990 | 0.902 | 1.4σ | 1.13 |
| `t7__S-separable.png` | 26.02 | 24.84 | -4.5% | 22.8 | 1.011 | 0.887 | 1.5σ | 1.08 |
| `t7__S-true.png` | 26.02 | 25.92 | -0.4% | 22.9 | 0.989 | 0.947 | 0.5σ | 1.07 |

Truth for the rest: φ_b = 27°, S = 1, bar_frac = 1. Orange marks the truth in each panel.

φ_b stays low by ~4° throughout — the one bias scheme B does not remove. It matches the +8 to +11°
rotation of the v_R quadrupole node that a 10% distance error produces (PROGRESS, 2026-10-29), so it
is that systematic seen through the posterior rather than a separate problem.

## Regenerating

```
cd scripts
# estimate S and rescore a finished search (uses the saved trajectories, no new integrations):
python rescore_empirical_selection.py --run <run> --estimator ratio --nside 16 --dbin 1.0
python rescore_empirical_selection.py --run <run> --estimator ipf   --nside 32
# then the posterior:
python test1_gp_emulator.py --run <run> --column logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_<variant> --loo-top 32
```

Variants: `asis`, `sfratio`, `sfemp`, `sf`. Summary plot:
`figures/test2_pattern_speed_recovery_schemeB_withS.png`.

`S-none` and `S-true` are copies of the same files in `figures/test2_schemeB/`, kept here so the
four treatments can be compared side by side.
