"""Why does the distance-error phase drift track corotation? Test the sigma/R mechanism.

R_CR correlating with the drift has no mechanism of its own - a bar does not resist distance
smearing by rotating slowly. The candidate mechanism is geometric:

  the distance error displaces each star along its line of sight, which within a few kpc of the
  GC is close to the Sun-GC axis, so a displacement of sigma ~ 0.125 * 8 kpc ~ 1 kpc scrambles
  AZIMUTH by roughly sigma / R.

If that is what is happening, the drift is set by the radius at which the kinematic signal lives,
and R_CR only predicts it because a bar's kinematic pattern reaches out towards corotation. The
test is then direct: measure, at zero distance error, where each galaxy's m = 2 signal actually
sits, and see whether it beats R_CR.

Three measurements per galaxy, all in the window the calibration uses, all at zero distance error:
  R_eff     amplitude-weighted mean radius of the per-ring m = 2 coefficients - "where the signal is"
  R_extent  width of the selected ring range in kpc
  twist     spread of the per-ring m = 2 phases across the window (deg), amplitude-weighted:
            an alternative mechanism, since mixing radii rotates the sum only if the phase varies
            with radius

Scored against the drift measured by diagnose_bar_angle_outliers.py.

Usage: python diagnose_drift_mechanism.py [--workers 4]
Table: results/bar_angle_drift_mechanism.csv
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

from common import RESULTS_DIR, rel
from bar_angle_vecsum_verification import ring_coeffs, Rc
from plot_smeared_velocity_maps import sample
from scan_bar_angle_windows import NBODY, AURIGA


def geometry(job):
    """Where the zero-error m = 2 signal sits, and how much its phase twists across the window."""
    gname, seed, frac_thr = job
    w = sample(gname, 0.0, seed)
    a2, b2 = ring_coeffs(w, "absvR")
    z = a2 + 1j * b2
    ok = np.isfinite(z)
    amp = np.where(ok, np.abs(z), np.nan)
    m = ok & (amp > frac_thr * np.nanmax(amp))          # the self:0.4 window
    if not m.any():
        return gname, dict(R_eff=np.nan, R_extent=np.nan, twist=np.nan, A_mean=np.nan)
    wgt = amp[m]
    R_eff = float(np.sum(wgt * Rc[m]) / np.sum(wgt))
    # per-ring m = 2 phase in degrees (axis angle, so defined mod 180); take the amplitude-weighted
    # spread about the summed phase, unwrapped onto (-90, 90]
    ph = 0.5 * np.degrees(np.angle(z[m]))
    ph0 = 0.5 * np.degrees(np.angle(np.sum(z[m])))
    d = (ph - ph0 + 90) % 180 - 90
    twist = float(np.sqrt(np.sum(wgt * d ** 2) / np.sum(wgt)))
    return gname, dict(R_eff=R_eff, R_extent=float(Rc[m].max() - Rc[m].min()),
                       twist=twist, A_mean=float(np.mean(wgt)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--thr", type=float, default=0.4, help="self:<thr> window")
    ap.add_argument("--snaps", nargs="+", default=NBODY + AURIGA)
    args = ap.parse_args()

    jobs = [(g, 11, args.thr) for g in args.snaps]
    print(f"zero-error m=2 geometry for {len(jobs)} galaxies", flush=True)
    rows = []
    with ProcessPoolExecutor(args.workers) as ex:
        for g, d in ex.map(geometry, jobs):
            rows.append(dict(galaxy=g, **d))
            print(f"  {g:14s} R_eff {d['R_eff']:5.2f} kpc  extent {d['R_extent']:4.2f}  "
                  f"twist {d['twist']:5.2f} deg", flush=True)
    G = pd.DataFrame(rows)

    D = pd.read_csv(RESULTS_DIR / "bar_angle_outlier_diagnosis.csv")
    M = D.merge(G, on="galaxy")
    S = pd.read_csv(RESULTS_DIR / "auriga_bar_survey.csv").set_index("halo")
    Q = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    au = M[M.suite == "auriga"].copy()
    h = [int(g[2:]) for g in au.galaxy]
    au["R_CR"] = Q.loc[h, "R_CR"].values
    au["R_bar"] = S.loc[h, "R_bar"].values
    au["A2_max"] = S.loc[h, "A2_max"].values
    au["inv_R_eff"] = 1.0 / au.R_eff
    M.to_csv(RESULTS_DIR / "bar_angle_drift_mechanism.csv", index=False)

    print(f"\npredictors of the 12.5% phase drift, {len(au)} Auriga galaxies\n")
    print(f"  {'quantity':24s} {'r':>7s} {'p':>8s} {'rho':>7s}")
    print("  " + "-" * 50)
    res = []
    for c in ["inv_R_eff", "R_eff", "R_CR", "twist", "R_extent", "R_bar", "A2_max", "A_mean"]:
        r = stats.pearsonr(au[c], au.absvR_drift)
        sp = stats.spearmanr(au[c], au.absvR_drift)
        res.append((abs(r.statistic), c, r, sp))
    for _, c, r, sp in sorted(res, reverse=True):
        print(f"  {c:24s} {r.statistic:+7.2f} {r.pvalue:8.4f} {sp.statistic:+7.2f}")

    # Does R_CR still say anything once R_eff is known? Partial correlation.
    def partial(x, y, z):
        rx = x - np.polyval(np.polyfit(z, x, 1), z)
        ry = y - np.polyval(np.polyfit(z, y, 1), z)
        return stats.pearsonr(rx, ry)
    p1 = partial(au.R_CR.values, au.absvR_drift.values, au.R_eff.values)
    p2 = partial(au.R_eff.values, au.absvR_drift.values, au.R_CR.values)
    print(f"\n  partial corr(R_CR,  drift | R_eff) = {p1.statistic:+.2f} (p = {p1.pvalue:.3f})")
    print(f"  partial corr(R_eff, drift | R_CR ) = {p2.statistic:+.2f} (p = {p2.pvalue:.3f})")
    rr = stats.pearsonr(au.R_CR, au.R_eff)
    print(f"  corr(R_CR, R_eff) = {rr.statistic:+.2f} (p = {rr.pvalue:.4f})")

    # The sigma/R prediction is quantitative, not just monotonic: drift ~ k * sigma / R_eff.
    sig = 0.125 * 8.0
    pred = np.degrees(sig / au.R_eff)
    k = np.sum(pred * au.absvR_drift) / np.sum(pred ** 2)
    print(f"\n  sigma/R_eff prediction (sigma = {sig:.1f} kpc), fitted scale k = {k:.2f}:")
    print(f"  {'galaxy':8s} {'R_eff':>6s} {'deg(sigma/R)':>13s} {'predicted':>10s} {'measured':>9s}")
    for _, x in au.sort_values("R_eff").iterrows():
        pv = np.degrees(sig / x.R_eff)
        print(f"  {x.galaxy:8s} {x.R_eff:6.2f} {pv:13.2f} {k * pv:10.2f} {x.absvR_drift:9.2f}")
    resid = au.absvR_drift - k * pred
    print(f"  rms residual {np.sqrt(np.mean(resid ** 2)):.2f} deg "
          f"against a {au.absvR_drift.std():.2f} deg spread in the drift")
    print(f"\nwrote {rel(RESULTS_DIR / 'bar_angle_drift_mechanism.csv')}")


if __name__ == "__main__":
    main()
