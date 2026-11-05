"""Is each Auriga galaxy isolated at z = 0, or interacting with a companion?

An interacting galaxy violates the steady-state assumption the orbit-library modelling rests on, and
a companion also drives its own m = 2 response that would contaminate a bar-angle measurement.
Fragkoudi et al. (2025) flag interacting Auriga galaxies by eye, from disturbances in the stellar
distribution within a 500 Myr lookback. This is a quantitative stand-in using the Subfind catalogue
at z = 0: for each main halo, find the most massive companion subhalo within R200 and report its
stellar-mass ratio and separation.

A galaxy is flagged as disturbed when a companion exceeds MASS_RATIO of the host's stellar mass
inside SEP_R200 of R200 - i.e. massive enough and close enough to matter dynamically. The thresholds
are deliberately generous: this is a screen for exclusion, not a merger census.

Usage: python check_auriga_interactions.py [--halos 2 9 ...]
Table: results/auriga_interactions.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import auriga as au
from common import RESULTS_DIR, rel
from survey_auriga_bars import ROOT

MASS_RATIO = 0.05        # companion stellar mass / host stellar mass
SEP_R200 = 1.0           # companion separation, in units of the host's R200


def companions(halo, root):
    import h5py
    hd = Path(root) / f"halo_{halo}"
    hdr = au.read_header(hd / "snapdir_127")
    h = au.main_halo(hd / "groups_127", hdr)
    u = au._units(hdr) if hasattr(au, "_units") else None
    pos, mstar = [], []
    for fn in sorted((hd / "groups_127").glob("fof_subhalo_tab_*")):
        with h5py.File(fn, "r") as f:
            if "Subhalo" not in f or "SubhaloPos" not in f["Subhalo"]:
                continue
            pos.append(f["Subhalo"]["SubhaloPos"][:])
            mstar.append(f["Subhalo"]["SubhaloMassType"][:, 4])
    pos = np.concatenate(pos) * (hdr["HubbleParam"] ** -1) * hdr["Time"] * 1e3   # kpc
    mstar = np.concatenate(mstar) * (hdr["HubbleParam"] ** -1) * 1e10            # Msun
    d = np.linalg.norm(pos - h["pos"], axis=1)
    host = int(np.argmin(d))
    m_host = mstar[host]
    other = np.ones(len(d), bool)
    other[host] = False
    near = other & (d < SEP_R200 * h["R200"]) & (mstar > 0)
    if not near.any():
        return dict(halo=halo, M_star_host=float(m_host), n_comp=0, best_ratio=0.0,
                    best_sep=np.nan, R200=float(h["R200"]), disturbed=False)
    ratio = mstar[near] / m_host
    k = int(np.argmax(ratio))
    return dict(halo=halo, M_star_host=float(m_host), n_comp=int((ratio > MASS_RATIO).sum()),
                best_ratio=float(ratio[k]), best_sep=float(d[near][k]), R200=float(h["R200"]),
                disturbed=bool(ratio[k] > MASS_RATIO))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--halos", nargs="+", type=int, default=None)
    ap.add_argument("--root", default=ROOT)
    args = ap.parse_args()
    if args.halos:
        halos = args.halos
    else:
        s = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv")
        halos = sorted(s.halo)
    rows = []
    for h in halos:
        try:
            r = companions(h, args.root)
            rows.append(r)
            print(f"  Au{h:<3d} M* {r['M_star_host']:.2e}  most massive companion within R200: "
                  f"{r['best_ratio']:.3f} of host at {r['best_sep']:.0f} kpc  "
                  f"({r['n_comp']} above {MASS_RATIO:.0%})  "
                  f"{'DISTURBED' if r['disturbed'] else 'isolated'}", flush=True)
        except Exception as e:
            print(f"  Au{h:<3d} FAILED: {e.__class__.__name__}: {e}", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(RESULTS_DIR / "auriga_interactions.csv", index=False)
    print(f"\nwrote {rel(RESULTS_DIR / 'auriga_interactions.csv')}")
    print(f"  {int(D.disturbed.sum())} of {len(D)} have a companion above "
          f"{MASS_RATIO:.0%} within R200")


if __name__ == "__main__":
    main()
