"""Can the saved orbit library support a frequency-ratio x1 cut? No - it is the cadence, not the baseline.

The dynamical bar length (Petersen et al. 2021; Zhang et al. 2024) needs orbits classified by their
frequencies in the bar frame. Zhang et al. (2024, MNRAS 534, 1175; arXiv:2406.06678) call an orbit
bar-supporting if Omega_R / Omega_x = 2 +- 0.2 with r_apo < 7 kpc, which is the cut used here. The
obvious question is whether the project's saved library - traj_NNN.npy, 1 Gyr at 41 snapshots, bar
frame, already integrated - can be re-used for that instead of integrating anything new.

The test integrates one set of orbits ONCE at 2 Gyr / 2001 samples (dt = 1 Myr), then re-measures
the frequencies on thinnings of that identical trajectory. Nothing but the sampling changes, so any
disagreement is the sampling alone. Two thinnings separate the two candidate causes:
  * 1 Gyr at 1 Myr   - the library's BASELINE, a fine cadence  (frequency resolution df = 1/Gyr)
  * 2 Gyr at 25 Myr  - the library's CADENCE, a full baseline  (df = 0.5/Gyr, better than reference)

Frequencies follow /data/hz420/OSARGs/x1_orbits_obs.py, which is the Zhang et al. (2024) code: the
highest peak of the FFT periodogram of the complex signal (x + i v_x) and (R + i v_R).

Usage: python diagnose_x1_frequency_sampling.py [--galaxies nbody au20 au26] [--n-star 600]
Table:  results/x1_frequency_sampling.csv
"""
import argparse

import numpy as np
import pandas as pd
import h5py
import agama
from scipy.signal import find_peaks

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, rel
import steady_state_likelihood as ssl
import potential_family as pf

RATIO, TOL, R_APO_MAX = 2.0, 0.2, 7.0          # Zhang et al. (2024) x1 cut
T_REF_GYR, N_REF = 2.0, 2001                   # dt = 1 Myr


def freq(t, y):
    """Highest peak of the FFT periodogram (OSARGs `extract_orbital_period`, n_peaks = 1)."""
    n = len(t)
    f = np.fft.fftfreq(n, t[1] - t[0])
    p = (1.0 / n) * abs(np.fft.fft(y)) ** 2
    m = f > 0
    f, p = f[m], p[m]
    pk, _ = find_peaks(p)
    return float(f[pk][np.argmax(p[pk])]) if len(pk) else np.nan


def measure(tr, t):
    """(Omega_R/Omega_x, Omega_x, Omega_R, r_apo) for a (N, nt, 6) bar-frame trajectory."""
    x, y, vx, vy = tr[:, :, 0], tr[:, :, 1], tr[:, :, 3], tr[:, :, 4]
    R = np.hypot(x, y)
    vR = (vx * x + vy * y) / np.where(R > 0, R, 1)
    fx = np.array([freq(t, x[i] + 1j * vx[i]) for i in range(len(tr))])
    fR = np.array([freq(t, R[i] + 1j * vR[i]) for i in range(len(tr))])
    rapo = np.sqrt(x ** 2 + y ** 2 + tr[:, :, 2] ** 2).max(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return fR / fx, fx, fR, rapo


def is_x1(ratio, rapo):
    return (np.abs(ratio - RATIO) <= TOL) & (rapo < R_APO_MAX)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--galaxies", nargs="+", default=["nbody", "au20", "au26"])
    ap.add_argument("--n-star", type=int, default=600)
    args = ap.parse_args()
    rng = np.random.default_rng(0)
    rows = []

    for gal in args.galaxies:
        Om = abs(pf.omega_true(gal))
        with h5py.File(MOCK_DIR / f"{gal}_mock_selection_sigmaDM0.271_N50k.hdf5", "r") as h:
            idx = np.sort(rng.choice(h["x_true"].shape[0], args.n_star, replace=False))
            w = np.column_stack([h[k][:][idx] for k in
                                 ("x_true", "y_true", "z_true", "vx_true", "vy_true", "vz_true")])
        pot = ssl.bar_aligned_potential(1.0, 1.0, gal)
        orbs = agama.orbit(potential=pot, ic=ssl.rotate_z(w, np.deg2rad(BAR_ANGLE_DEG)),
                           time=T_REF_GYR / ssl.TIME_UNIT_GYR, trajsize=N_REF, Omega=-Om,
                           accuracy=ssl.ACCURACY, **ssl.ORBIT_KW)
        tr = np.stack(orbs[:, 1])
        t = np.linspace(0, T_REF_GYR, N_REF)
        half = N_REF // 2 + 1                                 # the first 1 Gyr
        step25, step50 = (N_REF - 1) // (T_REF_GYR * 40), (N_REF - 1) // (T_REF_GYR * 20)
        s25, s50 = int(step25), int(step50)
        sets = {
            "reference      2 Gyr /  1 Myr": (tr, t),
            "baseline only  1 Gyr /  1 Myr": (tr[:, :half], t[:half]),
            "cadence only   2 Gyr / 25 Myr": (tr[:, ::s25], t[::s25]),
            "library saved  1 Gyr / 25 Myr": (tr[:, :half:s25], t[:half:s25]),
            "library scored 1 Gyr / 50 Myr": (tr[:, :half:s50], t[:half:s50]),
        }
        ref = None
        print(f"\n########## {gal}  Omega_b = {Om:.2f} km/s/kpc")
        for lab, (a, tt) in sets.items():
            ratio, fx, fR, rapo = measure(a, tt)
            sel = is_x1(ratio, rapo)
            dt, T = tt[1] - tt[0], tt[-1] - tt[0]
            if ref is None:
                ref = (ratio, sel, rapo, fx, fR)
                rec, cont, dmed = len(ratio), 0, 0.0
                print(f"  x1 frequencies: Omega_x median {np.nanmedian(fx[sel]):5.1f}/Gyr, "
                      f"Omega_R median {np.nanmedian(fR[sel]):5.1f}/Gyr "
                      f"(99th {np.nanpercentile(fR[sel], 99):.0f})")
            else:
                g = np.isfinite(ratio) & np.isfinite(ref[0])
                dmed = float(np.nanmedian(np.abs(ratio[g] - ref[0][g])))
                rec = int((sel & ref[1]).sum())
                cont = int((sel & ~ref[1]).sum())
            alias = float(np.nanmean(ref[4][ref[1]] > 1 / (2 * dt)))
            print(f"  {lab}:  f_Nyq {1/(2*dt):6.1f}/Gyr  df {1/T:4.2f}/Gyr | x1 {int(sel.sum()):4d}"
                  f" | recovered {rec:4d}/{int(ref[1].sum())}  false {cont:3d}"
                  f" | median |d(ratio)| {dmed:4.2f} | Omega_R aliased {100*alias:3.0f}%")
            rows.append(dict(galaxy=gal, Omega_b=Om, sampling=lab.strip(), dt_Myr=1000 * dt,
                             T_Gyr=T, f_nyq=1 / (2 * dt), df=1 / T, n_x1=int(sel.sum()),
                             n_recovered=rec, n_false=cont, med_abs_dratio=dmed,
                             frac_aliased=alias,
                             R_x1_max=float(np.hypot(tr[sel, :, 0], tr[sel, :, 1]).max())
                             if sel.any() else np.nan))
    D = pd.DataFrame(rows)
    D.to_csv(RESULTS_DIR / "x1_frequency_sampling.csv", index=False)
    print(f"\nwrote {rel(RESULTS_DIR / 'x1_frequency_sampling.csv')}")


if __name__ == "__main__":
    main()
