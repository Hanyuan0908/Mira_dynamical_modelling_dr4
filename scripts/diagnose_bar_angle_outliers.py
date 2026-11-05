"""Why are Au5 and Au20 outliers? Split their angle error into the two things that can cause it.

A galaxy's measured bar angle can be wrong for two independent reasons, and they call for different
fixes:

  intrinsic  the kinematic m = 2 quadrupole is not aligned with the bar's density major axis even
             with perfect distances. Nothing about the observation causes it; it is the galaxy.
  drift      the distance error mixes azimuths asymmetrically about the Sun-GC line and rotates
             the measured pattern. This is what the calibration is built to undo.

The calibration fits one (slope, intercept) to all galaxies, so it can only absorb the part of each
that is common to them. A galaxy with an unusual intrinsic misalignment, or an unusual drift, lands
off the relation - and the two look identical in the transfer plot while having nothing in common.

Everything is measured in the galaxy's native frame (bar at BAR_ANGLE_DEG), with the selection
function applied, averaged over n_real distance-error realisations, in the window the calibration
actually uses.

Usage: python diagnose_bar_angle_outliers.py [--n-real 6] [--window kin0.4] [--workers 4]
Table: results/bar_angle_outlier_diagnosis.csv
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

from common import RESULTS_DIR, BAR_ANGLE_DEG, rel
from bar_angle_vecsum_verification import phase
from plot_smeared_velocity_maps import sample
from scan_bar_angle_windows import NBODY, AURIGA, phi_from, CACHE


def one(job):
    gname, frac, seed, window = job
    w = sample(gname, frac, seed)
    win = window if window.startswith("kin") else tuple(float(v) for v in window.split("-"))
    # phase() returns the bar axis in the input frame; +BAR_ANGLE_DEG turns it into
    # (measured - true), so 0 means the estimator is unbiased. Without it every intrinsic
    # value reads as about -27 and looks like a catastrophe rather than a degree or two.
    return gname, frac, {e: phase(w, e, win) + BAR_ANGLE_DEG for e in ("vR", "absvR")}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-real", type=int, default=6)
    ap.add_argument("--window", default="kin0.4")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--snaps", nargs="+", default=NBODY + AURIGA)
    args = ap.parse_args()

    # frac = 0 is deterministic apart from the selection draw, so it needs fewer realisations
    jobs = ([(g, 0.0, 11 + s, args.window) for g in args.snaps for s in range(2)]
            + [(g, 0.125, 11 + s, args.window) for g in args.snaps for s in range(args.n_real)])
    print(f"{len(jobs)} phase measurements over {len(args.snaps)} galaxies", flush=True)
    acc = {}
    with ProcessPoolExecutor(args.workers) as ex:
        for k, (g, fr, ph) in enumerate(ex.map(one, jobs), 1):
            for e, v in ph.items():
                acc.setdefault((g, fr, e), []).append(v)
            if k % 20 == 0:
                print(f"  {k}/{len(jobs)}", flush=True)

    # the transfer offset each galaxy shows in the calibration, from the ring cache
    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], [str(x) for x in d["snaps"]], list(d["angles"])
    P = np.array([[np.nanmean([phi_from(Z[gi, ai, s, 1], "self:0.4", Z[gi, ai, s, 0], "absvR")
                               for s in range(Z.shape[2])]) for ai in range(len(angles))]
                  for gi in range(len(snaps))])
    offset = dict(zip(snaps, np.nanmean(P - P.mean(0), axis=1)))

    rows = []
    for g in args.snaps:
        r = dict(galaxy=g, suite="auriga" if g.startswith("au") else "nbody",
                 transfer_offset=offset.get(g, np.nan))
        for e in ("vR", "absvR"):
            p0 = float(np.mean(acc[(g, 0.0, e)]))
            p1 = float(np.mean(acc[(g, 0.125, e)]))
            r[f"{e}_intrinsic"] = p0
            r[f"{e}_smeared"] = p1
            r[f"{e}_drift"] = p1 - p0
            r[f"{e}_drift_sem"] = float(np.std(acc[(g, 0.125, e)], ddof=1)
                                        / np.sqrt(len(acc[(g, 0.125, e)])))
        rows.append(r)
    D = pd.DataFrame(rows)
    out = RESULTS_DIR / "bar_angle_outlier_diagnosis.csv"
    D.to_csv(out, index=False)

    print(f"\nangle budget in the galaxy's native frame, window {args.window} "
          f"(truth = 0 by construction)\n")
    print(f"  {'galaxy':10s} {'suite':7s} | {'absvR intrinsic':>15s} {'drift':>8s} "
          f"{'= smeared':>10s} | {'transfer offset':>15s}")
    print("  " + "-" * 78)
    for _, r in D.sort_values("absvR_intrinsic").iterrows():
        print(f"  {r.galaxy:10s} {r.suite:7s} | {r.absvR_intrinsic:+15.2f} "
              f"{r.absvR_drift:+8.2f} {r.absvR_smeared:+10.2f} | {r.transfer_offset:+15.2f}")

    au = D[D.suite == "auriga"]
    print()
    for col in ("absvR_intrinsic", "absvR_drift", "absvR_smeared"):
        rr = stats.pearsonr(au[col], au.transfer_offset)
        ss = stats.spearmanr(au[col], au.transfer_offset)
        print(f"  corr({col:17s}, transfer offset) over 13 Auriga: "
              f"r = {rr.statistic:+.2f} (p = {rr.pvalue:.4f}), rho = {ss.statistic:+.2f}")
    print(f"\n  spread over the 13 Auriga: intrinsic "
          f"{au.absvR_intrinsic.min():+.1f} to {au.absvR_intrinsic.max():+.1f} deg, "
          f"drift {au.absvR_drift.min():+.1f} to {au.absvR_drift.max():+.1f} deg")
    print(f"\nwrote {rel(out)}")


if __name__ == "__main__":
    main()
