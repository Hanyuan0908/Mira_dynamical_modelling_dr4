"""Do the Milky Way stars complete fewer orbits / bar passages in the library window than the mock stars? (user,
2026-11-24: could the weak pattern-speed signal be an integration-time effect?)

The libraries use the snapshots at 0.25-1 Gyr of a 1 Gyr integration for every potential. If the Milky Way rotation curve
gives longer orbital or bar-relative periods, its stars sample fewer bar encounters in that window. For N_SUB random stars
per case, integrated for 1 Gyr at 1 Myr in the bar frame at the reference bar (S = bar_frac = 1), over 0.25-1 Gyr:
  n_peri   radial oscillations (pericentre passages: local minima of R)
  n_rev    revolutions relative to the bar (|unwrapped bar-frame azimuth change| / 2 pi); each passes the bar twice.
           Stars trapped by the bar (librating, not circulating) have n_rev ~ 0.
reported overall and in bins of the star's mean R, with the fraction of stars with n_rev < 0.5 ("not circulating").

Usage: python count_orbits_in_library_window.py [--n-sub 3000]
Output: results/orbit_counts_in_library_window.csv and a printed table
"""
import argparse

import numpy as np
import pandas as pd
import agama

from common import RESULTS_DIR
import steady_state_likelihood as ssl
import test1_nbody_qmc as old
from potential_family import TIME_UNIT_GYR, ORBIT_KW
from library_noise_floor_test import load_stars

CASES = [("OSARG stars, Sormani+22", "osarg", "sormani22", 32.5), ("steady-state mock, Hunter+24", "h24mock", "mw", 32.5),
         ("t7 stars, t7 family", "t7cov", "nbody", 26.0), ("OSARG stars, t7 family", "osarg", "nbody", 26.0)]
RBINS = [(0, 2), (2, 4), (4, 6), (6, 8)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-sub", type=int, default=3000)
    args = ap.parse_args()
    rows = []
    for lab, stars, family, om in CASES:
        w, _, _, phi = load_stars(stars)
        sel = np.random.default_rng(5).choice(len(w), min(args.n_sub, len(w)), replace=False)
        pot = ssl.bar_aligned_potential(1.0, 1.0, family)
        orbs = agama.orbit(potential=pot, ic=ssl.rotate_z(w[sel], np.radians(phi)), time=1.0 / TIME_UNIT_GYR, trajsize=1001,
                           Omega=-abs(om), accuracy=old.ACCURACY, **ORBIT_KW)
        t = orbs[0, 0] * TIME_UNIT_GYR
        X = np.stack(orbs[:, 1])                                # (n, 1001, 6) bar frame
        k = t >= 0.25 - 1e-9
        R = np.hypot(X[:, k, 0], X[:, k, 1])
        ph = np.unwrap(np.arctan2(X[:, k, 1], X[:, k, 0]), axis=1)
        n_rev = np.abs(ph[:, -1] - ph[:, 0]) / (2 * np.pi)
        n_peri = np.sum((R[:, 1:-1] < R[:, :-2]) & (R[:, 1:-1] < R[:, 2:]), axis=1)
        Rm = R.mean(axis=1)
        for lo, hi in RBINS + [(0, 8)]:
            m = (Rm >= lo) & (Rm < hi)
            if m.sum() < 20:
                continue
            q = lambda a: np.percentile(a[m], [16, 50, 84])
            r_rev, r_per = q(n_rev), q(n_peri)
            rows.append(dict(case=lab, R_lo=lo, R_hi=hi, n=int(m.sum()), rev_p16=r_rev[0], rev_p50=r_rev[1], rev_p84=r_rev[2],
                             bar_passages_p50=2 * r_rev[1], peri_p16=r_per[0], peri_p50=r_per[1], peri_p84=r_per[2],
                             frac_not_circulating=float(np.mean(n_rev[m] < 0.5))))
            print(f"{lab:30s} R {lo}-{hi} kpc (n {m.sum():4d}): bar revolutions {r_rev[1]:.1f} ({r_rev[0]:.1f}-{r_rev[2]:.1f}), "
                  f"i.e. ~{2 * r_rev[1]:.0f} bar passages; radial oscillations {r_per[1]:.0f} ({r_per[0]:.0f}-{r_per[2]:.0f}); "
                  f"not circulating {np.mean(n_rev[m] < 0.5):.2f}", flush=True)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "orbit_counts_in_library_window.csv", index=False)


if __name__ == "__main__":
    main()
