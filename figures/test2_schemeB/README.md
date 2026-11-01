# Scheme B corner plots — Test 2, 10% distance error

Posterior corner plots (GP emulator, 256 Sobol points) for the **scheme B** χ²,
`fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024`:

> 0.3 kpc rings to R < 6 kpc × |z| layers 0–1, 1–2, 2–3 kpc. Per block it fits
> ⟨v_R⟩ at m = 2,4; ⟨v_φ⟩ at m = 0,2,4; σ_R² and σ_φ² at **m = 2,4 only** —
> i.e. the dispersion *ring averages* are dropped. Those terms are where the
> over-hot orbit library bites, and keeping them makes the fit buy agreement in
> library temperature with a spuriously fast bar (PROGRESS, 2026-10-31).

Orange lines / squares mark the truth. File names are `<snapshot>_<mock>__lib-<weighting>.png`.

## Library weightings

| tag | library smeared? | selection weighting? | right for |
|---|---|---|---|
| `lib-asis` | no | no | nothing — kept as a control |
| `lib-observed` | yes | no | mocks **without** a selection function |
| `lib-sf` | yes | yes, S(x)/τ | mocks **with** a selection function |

## What is here

| file | truth Ω_b | recovered | bias | φ_b | S | bar_frac | worst | LOO | matched model |
|---|---|---|---|---|---|---|---|---|---|
| `t2_selection__lib-asis.png` | 44.25 | 40.93 | -7.5% | 21.7 | 0.930 | 1.388 | 1.8σ | 1.30 |  |
| `t2_selection__lib-sf.png` | 44.25 | 43.40 | -1.9% | 18.8 | 1.009 | 0.817 | 2.4σ | 0.82 | yes |
| `t3_selection__lib-asis.png` | 39.53 | 38.14 | -3.5% | 27.6 | 0.849 | 1.695 | 3.5σ | 1.76 |  |
| `t3_selection__lib-sf.png` | 39.53 | 40.63 | +2.8% | 24.5 | 0.953 | 0.939 | 0.6σ | 0.86 | yes |
| `t4_selection__lib-asis.png` | 33.96 | 30.24 | -10.9% | 19.3 | 0.896 | 1.216 | 2.8σ | 1.19 |  |
| `t4_selection__lib-sf.png` | 33.96 | 34.86 | +2.7% | 23.8 | 1.104 | 0.853 | 0.6σ | 0.78 | yes |
| `t5.5_selection__lib-asis.png` | 28.84 | 26.32 | -8.7% | 22.4 | 0.911 | 1.347 | 3.1σ | 1.37 |  |
| `t5.5_selection__lib-sf.png` | 28.84 | 28.06 | -2.7% | 21.3 | 0.991 | 0.897 | 0.9σ | 1.08 | yes |
| `t7_selection__lib-asis.png` | 26.02 | 24.54 | -5.7% | 26.3 | 0.913 | 1.486 | 1.8σ | 1.09 |  |
| `t7_noselection__lib-asis.png` | 26.02 | 25.80 | -0.8% | 25.6 | 0.908 | 1.473 | 1.1σ | 1.43 |  |
| `t7_noselection__lib-observed.png` | 26.02 | 26.73 | +2.7% | 27.0 | 1.005 | 1.069 | 0.5σ | 1.60 | yes |
| `t7_selection__lib-sf.png` | 26.02 | 25.92 | -0.4% | 22.9 | 0.989 | 0.947 | 0.5σ | 1.07 | yes |

Truth for the other three parameters is φ_b = 27°, S = 1, bar_frac = 1.

## Regenerating

```
cd scripts
python test1_gp_emulator.py --run <run> \
    --column logLs_fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024_<variant> --loo-top 32
```

That writes into `figures/` with the long auto-generated name; move it here and rename.
`common.savefig` now accepts a sub-directory in the name, so new scripts can write
straight into `figures/test2_schemeB/...`.

Scheme A (the main scheme, dispersion ring averages kept) corner plots remain in `figures/`
under their original names, for comparison.
