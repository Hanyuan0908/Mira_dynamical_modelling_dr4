"""Which kinematic window makes the bar-angle calibration depend least on the galaxy?

The window is the free knob: the measured m = 2 phase is the amplitude-weighted sum over whichever
rings the rule keeps, and the rule can be anything. What we want from it is NOT the narrowest band
but the weakest dependence on bar properties, because that dependence is the thing that stops the
calibration transferring between galaxies (PROGRESS 2026-11-04/05, still an open question).

Three numbers per rule, all computed from results/bar_angle_ring_cache.npz, so scanning costs
nothing - the per-ring coefficients are already measured for 25 galaxies x 6 angles x 8
realisations:

  band / slope   the Auriga 16-84% width and the median's slope -> the formal +-phi_b
  blind rms      calibrate on the 13 Auriga, infer the 12 N-body snapshots through that median,
                 compare with their known truth. Includes the suite-to-suite offset, which the
                 band does not.
  rho(R_CR)      Spearman of each Auriga galaxy's mean offset against its corotation radius. This
                 is the bar-property dependence itself; |rho| near 0 is the target.

Usage: python scan_kinematic_window_transfer.py [--top 6]
Table: results/bar_angle_window_transfer_scan.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from common import RESULTS_DIR, rel
from scan_bar_angle_windows import phi_from, CACHE, NBODY, AURIGA

EST = ("vR", "absvR")


def measured(Z, snaps, angles, ei, rule, est):
    """Mean measured angle over realisations, shape (n_galaxy, n_angle)."""
    return np.array([[np.nanmean([phi_from(Z[g, a, s, ei], rule, Z[g, a, s, 0], est)
                                  for s in range(Z.shape[2])])
                      for a in range(len(angles))] for g in range(len(snaps))])


def invert(ang, curve, m):
    if curve[0] > curve[-1]:
        return np.interp(m, curve[::-1], ang[::-1])
    return np.interp(m, curve, ang)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=6)
    args = ap.parse_args()

    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], [str(x) for x in d["snaps"]], list(d["angles"])
    ang = np.array(angles)
    iau = [i for i, g in enumerate(snaps) if g.startswith("au")]
    inb = [i for i, g in enumerate(snaps) if not g.startswith("au")]
    Q = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    rcr = np.array([float(Q.loc[int(snaps[i][2:]), "R_CR"]) for i in iau])

    RULES = ([f"self:{f:g}" for f in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
             + [f"vR:{f:g}" for f in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)]
             + [f"peakfrac:{a:g}-{b:g}" for a, b in ((0.3, 1.2), (0.4, 1.5), (0.5, 2.0),
                                                     (0.2, 1.0), (0.6, 2.5))]
             + ["fixed:1-4", "fixed:1-3", "fixed:0.5-3", "fixed:2-6"])

    rows = []
    for ei, est in enumerate(EST):
        for rule in RULES:
            P = measured(Z, snaps, angles, ei, rule, est)
            if not np.isfinite(P).all():
                continue
            A = P[iau]
            med = np.nanmedian(A, 0)
            width = float(np.mean(np.nanpercentile(A, 84, 0) - np.nanpercentile(A, 16, 0)))
            slope = float(np.polyfit(ang, med, 1)[0])
            if slope <= 0.05:
                continue
            err = np.array([invert(ang, med, P[i]) - ang for i in inb]).ravel()
            off = np.nanmean(A - A.mean(0), axis=1)          # per-Auriga-galaxy offset
            rho = float(stats.spearmanr(rcr, off).statistic)
            rows.append(dict(est=est, rule=rule, band=width, slope=slope,
                             phi_err=0.5 * width / slope, blind_mean=float(err.mean()),
                             blind_rms=float(np.sqrt((err ** 2).mean())),
                             blind_max=float(np.abs(err).max()), rho_RCR=rho))

    T = pd.DataFrame(rows)
    T.to_csv(RESULTS_DIR / "bar_angle_window_transfer_scan.csv", index=False)
    for est in EST:
        S = T[T.est == est].sort_values("blind_rms")
        print(f"\n=== {est}: ranked by blind-check rms on the 12 N-body snapshots ===")
        print(f"  {'rule':18s} {'band':>6s} {'slope':>6s} {'+-phi':>7s} | "
              f"{'blind mean':>10s} {'rms':>6s} {'max':>6s} | {'rho(R_CR)':>9s}")
        print("  " + "-" * 76)
        for _, r in S.head(args.top).iterrows():
            print(f"  {r.rule:18s} {r.band:6.2f} {r.slope:6.2f} {r.phi_err:7.2f} | "
                  f"{r.blind_mean:+10.2f} {r.blind_rms:6.2f} {r.blind_max:6.2f} | {r.rho_RCR:+9.2f}")
        w = S.reindex(S.rho_RCR.abs().sort_values().index).head(3)
        print(f"  weakest bar-property dependence |rho(R_CR)|:")
        for _, r in w.iterrows():
            print(f"    {r.rule:16s} rho {r.rho_RCR:+.2f}  blind rms {r.blind_rms:5.2f}  "
                  f"+-phi {r.phi_err:5.2f}")
    print(f"\nwrote {rel(RESULTS_DIR / 'bar_angle_window_transfer_scan.csv')}")


if __name__ == "__main__":
    main()
