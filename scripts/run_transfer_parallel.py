"""Measure the bar-angle transfer table in parallel over galaxies.

test_bar_angle_transfer_12p5.measure() is pure and independent per (galaxy, angle), so the
25 galaxies x 6 angles are farmed out to a process pool. Seeds are indexed by position in the
FULL SNAPS list, exactly as the serial script does, so the table comes out identical either way.
Serially the 25 galaxies take about six hours; this is what makes the Auriga extension practical.

Usage: python run_transfer_parallel.py [--workers 5] [--n-real 8]
Then:  python test_bar_angle_transfer_12p5.py --replot     (leave-one-out table + figure)
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from common import RESULTS_DIR, rel
from nbody import SNAPSHOT_DIR
import test_bar_angle_transfer_12p5 as T


def one(job):
    gi, gn, ai, a, n_real, kin = job
    r = T.measure(gn, a, n_real, 20_000 + 7919 * ai + 104_729 * gi, kin)
    rec = dict(galaxy=gn, t=T.t_of(gn), phi_true=a)
    for (e, w), (m, s) in r.items():
        rec[f"{e}_{w}"] = m
        rec[f"{e}_{w}_err"] = s
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--n-real", type=int, default=8)
    ap.add_argument("--kin", type=float, default=0.7)
    ap.add_argument("--snaps", nargs="+", default=T.SNAPS)
    args = ap.parse_args()

    tb = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")
    au_tb = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")

    jobs = [(T.SNAPS.index(gn), gn, ai, a, args.n_real, args.kin)
            for gn in args.snaps for ai, a in enumerate(T.ANGLES)]
    print(f"{len(jobs)} measurements over {len(args.snaps)} galaxies "
          f"on {args.workers} workers", flush=True)
    rows = []
    with ProcessPoolExecutor(args.workers) as ex:
        for k, rec in enumerate(ex.map(one, jobs), 1):
            rows.append(rec)
            print(f"  [{k:3d}/{len(jobs)}] {rec['galaxy']:14s} phi {rec['phi_true']:5.1f}  "
                  f"absvR/kin {rec['absvR_kin']:6.2f}", flush=True)

    D = pd.DataFrame(rows)
    om = []
    for gn in D.galaxy:
        if T.is_auriga(gn):
            h = int(str(gn)[2:])
            om.append(float(au_tb.loc[h, "Omega_p"]) if h in au_tb.index else np.nan)
        else:
            om.append(float(tb[np.isclose(tb.snap_time, T.t_of(gn))]["pattern_speed_fd"].iloc[0]))
    D.insert(3, "Omega_b", om)

    out = RESULTS_DIR / "bar_angle_transfer_12p5.csv"
    if out.exists():                      # keep galaxies measured on an earlier run
        old = pd.read_csv(out)
        D = pd.concat([old[~old.galaxy.isin(D.galaxy.unique())], D], ignore_index=True)
    D = D.sort_values(["galaxy", "phi_true"]).reset_index(drop=True)
    D.to_csv(out, index=False)
    print(f"\nwrote {rel(out)}: {len(D)} rows, {D.galaxy.nunique()} galaxies")


if __name__ == "__main__":
    main()
