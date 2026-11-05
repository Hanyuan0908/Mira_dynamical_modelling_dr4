"""Measure the bar in every Auriga z = 0 halo, to find galaxies that widen the calibration suite.

The bar-angle calibration band is currently built from twelve snapshots of ONE isolated N-body
galaxy plus Au18, so it samples a narrow range of bar strength (A2 = 0.32-0.35 as observed, 0.37-0.50
intrinsic). Since the residual scatter in the calibration tracks bar strength, a suite spanning a
wider range of A2 is what the band needs.

This measures, for each halo's z = 0 snapshot, the same quantities prepare_au18_mwframe measures for
Au18 - the m = 2 amplitude profile of the stars, its peak, the bar half-length and the bar phase -
so the numbers are directly comparable with Au18 and with the N-body sequence. Fragkoudi et al.
(2025, MNRAS, arXiv:2406.09453) classify Auriga galaxies as barred at A2 > 0.25; that threshold is
reported here but the measurement is our own.

Usage: python survey_auriga_bars.py [--halos 2 3 ...] [--root PATH]
Table: results/auriga_bar_survey.csv
"""
import argparse
import traceback

import numpy as np
import pandas as pd

import auriga as au
from common import RESULTS_DIR, rel
from prepare_au18_mwframe import shrinking_sphere, measure_bar, RMAX, R_DISC

ROOT = "/home/hz420/austreams/Auriga_simulation"
A2_BARRED = 0.25          # Fragkoudi et al. 2025 threshold


def one(halo, root):
    """Same recipe as prepare_au18_mwframe.main(): centre, de-spin, measure the m = 2 profile."""
    from pathlib import Path
    hd = Path(root) / f"halo_{halo}"
    hdr = au.read_header(hd / "snapdir_127")
    h = au.main_halo(hd / "groups_127", hdr)
    # read_particles adds Coordinates itself and returns them relative to the given centre
    st = au.read_particles(hd / "snapdir_127", 4,
                           ["Velocities", "Masses", "GFM_StellarFormationTime"],
                           h["pos"], RMAX, hdr)
    keep = st["GFM_StellarFormationTime"] > 0          # drop wind particles
    st = {k: v[keep] for k, v in st.items()}
    x = st["Coordinates"].astype(np.float64)
    v = st["Velocities"].astype(np.float64)
    m = st["Masses"].astype(np.float64)
    x -= shrinking_sphere(x, m)
    s = np.linalg.norm(x, axis=1) < R_DISC
    v -= np.average(v[s], axis=0, weights=m[s])
    L = np.sum(m[s, None] * np.cross(x[s], v[s]), axis=0)
    xr = x @ au.rotation_to(L, np.array([0.0, 0.0, -1.0])).T
    bar = measure_bar(xr[:, 0], xr[:, 1], xr[:, 2], m)
    disc = np.hypot(xr[:, 0], xr[:, 1]) < R_DISC
    return dict(halo=halo, n_star=int(len(m)),
                M_star_disc=float(m[disc].sum() * 1e10),
                M200=float(h["M200"]),
                A2_max=float(bar["A2_max"]), R_peak=float(bar["R_peak"]),
                R_bar=float(bar["R_bar"]), barred=bool(bar["A2_max"] > A2_BARRED))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--halos", nargs="+", type=int, default=None)
    ap.add_argument("--root", default=ROOT)
    args = ap.parse_args()
    from pathlib import Path
    halos = args.halos or sorted(int(p.name.split("_")[1]) for p in Path(args.root).glob("halo_*")
                                 if (p / "snapdir_127").is_dir())
    rows = []
    for h in halos:
        try:
            r = one(h, args.root)
            rows.append(r)
            print(f"  Au{h:<3d} A2 {r['A2_max']:.3f}  R_peak {r['R_peak']:5.2f}  "
                  f"R_bar {r['R_bar']:5.2f} kpc  {'BARRED' if r['barred'] else 'unbarred':8s} "
                  f"({r['n_star']:,} stars)", flush=True)
        except Exception as e:
            print(f"  Au{h:<3d} FAILED: {e.__class__.__name__}: {e}", flush=True)
    D = pd.DataFrame(rows).sort_values("A2_max", ascending=False)
    out = RESULTS_DIR / "auriga_bar_survey.csv"
    D.to_csv(out, index=False)
    print(f"\nwrote {rel(out)}")
    b = D[D.barred]
    print(f"  {len(b)} of {len(D)} barred at A2 > {A2_BARRED}; "
          f"A2 spans {b.A2_max.min():.3f}-{b.A2_max.max():.3f}")


if __name__ == "__main__":
    main()
