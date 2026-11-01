# Key results — steady-state bar test with Mira-like data

Two results from the Test 2 campaign (2026-10-29 to 11-01). Full history in `PROGRESS.md`.

Conditions throughout: isolated N-body bar, 50k stars, **10% distance error**, DR3-Mira-like
selection function, chi2 **scheme B** (dispersion ring averages dropped), library weighted by
S(x)/tau. Truth pattern speeds are the simulation directory's finite-difference column.

---

## 1. Pattern speed: recovered vs true  (`01_pattern_speed_recovery.png`)

The method recovers Omega_b across a factor-1.8 range of pattern speed, from a bar at 43.9 down
to 25.2 km/s/kpc.

| snapshot | true Omega_b | recovered | bias |
|---|---|---|---|
| t = 2 | 43.90 | 43.40 | -1.1% |
| t = 3 | 38.67 | 40.63 | +5.1% |
| t = 4 | 33.28 | 34.86 | +4.7% |
| t = 5.5 | 28.21 | 28.06 | -0.5% |
| t = 7 | 25.20 | 25.92 | +2.9% |

**Mean |bias| {np.mean(rows):.1f}%.** Error budget (mean over the five snapshots):

| parameter | statistical (Fisher) | systematic | sys/stat |
|---|---|---|---|
| Omega_b | 0.9% | +2.1% | 2.3 |
| phi_b | 1.2 deg | -4.8 deg | 3.9 |
| S (bar length) | 3.4% | **+0.8%** | **0.2** |
| bar_frac (strength) | 6.4% | -11.7% | 1.8 |

- **Bar length is recovered unbiased** - the one parameter the method measures cleanly.
- Strength is 12% low, pattern speed 2% high, both ~2 sigma systematics.
- The quoted (tempered) widths are ~4.7x the statistical error, so the tempering absorbs the
  systematic rather than exposing it; see PROGRESS 2026-11-01.

---

## 2. Bar angle from the v_R quadrupole  (`02_bar_angle_inversion_10pct.png`)

The orbit library recovers phi_b poorly (-4.8 deg bias, 3.9x its statistical error), so the bar
angle is better taken from the quadrupole in <v_R>, whose node lies along the bar axis. The
distance error rotates that node, so the relation must be calibrated:

| | fitted relation (t = 7) | node scatter | phi_b recoverable to |
|---|---|---|---|
| no selection | node = -0.664 phi_true + 0.42 | 2.75 deg | 4.1 deg |
| with selection | node = -0.669 phi_true + 1.29 | 2.14 deg | **3.2 deg** |

- **The relation is flatter than 1:1** (slope -0.67, not -1): the node under-reports the angle,
  and a given node error becomes a ~1.5x larger angle error.
- **The selection function does not matter** - same slope and scatter with and without it.
- **The calibration is not universal in the bar angle**: the bias roughly doubles from 15 to 35 deg
  (PROGRESS 2026-11-01), so the full surface node(phi_true, sigma_DM) is needed, not one curve.
- Snapshot dependence: see section 3.

---

## 3. Does the calibration transfer?  (`03_bar_angle_transfer_t7_to_t4.png`)

The t = 7 relation applied to a different bar (t = 4: A2 0.472 vs 0.397, R_bar 4.25 vs 5.00 kpc,
Omega_b 33.3 vs 25.2):

| | fitted relation |
|---|---|
| t = 7 (the calibration) | node = -0.669 phi + 1.29 |
| t = 4 (its own) | node = -0.639 phi + 1.79 |

**Mean error -1.86 deg, rms 1.89 deg**, against a 2.20 deg per-measurement scatter - so a single
measurement sits within ~1 sigma of the truth and the calibration does transfer in practice.

But the offset is systematic, not noise: all nine residuals share the same sign and scatter about
their mean by only 0.34 deg, so it will not average down. **Budget: ~1.9 deg systematic + 2.2 deg
statistical = ~2.9 deg**, versus the orbit library's -4.8 deg bias.

**Caveat.** t = 4 and t = 7 are the same galaxy 2.9 Gyr apart - far more similar to each other than
either is to the Milky Way. A simulation-derived calibration applied to real data would carry a
systematic of at least this size, probably larger. The method is not calibration-free.

---

## Reproducing

```
cd scripts
python plot_test2_pattern_speed_recovery.py                 # figure 1
python plot_bar_angle_inversion.py --frac 0.10 --n-real 15  # figure 2
python plot_bar_angle_inversion.py --galaxy nbody_t4 --frac 0.10 --n-real 15
python plot_bar_angle_transfer_test.py --ref nbody --test nbody_t4   # figure 3
python quantify_test2_errors.py                             # the error budget
```
