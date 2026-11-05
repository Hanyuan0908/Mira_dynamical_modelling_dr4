"""Does the galaxy-to-galaxy spread in phase drift survive with the selection function switched off?

The drift varies from -1 to +11 deg across the 13 Auriga galaxies and tracks corotation radius. Two
very different things could produce that:

  dynamics   the distance error scatters stars across the galaxy, and what they carry when they
             land in the measurement window depends on where they came from - in particular on
             whether they came from inside or outside corotation, where the bar-driven m = 2
             velocity response has opposite sign. R_CR would then matter on its own.
  selection  the DR3-Mira-like selection picks a different part of each galaxy, because the
             galaxies have different density distributions and extinction sightlines. The drift
             would then be a property of the selection-galaxy pairing, not of the bar, and R_CR
             would only be along for the ride.

Switching the selection off and keeping the distance error separates them. If the spread collapses,
the selection function is doing the work; if it survives, the effect is in the smearing itself.

Both arms use identical galaxies, seeds and window, so the comparison is paired.

Usage: python diagnose_drift_without_selection.py [--workers 4] [--n-real 6] [--nmax 300000]
Table: results/bar_angle_drift_no_selection.csv
"""
import argparse
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import h5py
from scipy import stats

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, rel
from galaxies import galaxy
from bar_angle_vecsum_verification import phase
from scan_bar_angle_windows import AURIGA
import observe, selection

FRAC = 0.125
WIN = "kin0.4"


def one(job):
    """Phase with and without the selection function, at zero and 12.5% distance error."""
    gname, seed, nmax = job
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{FRAC * 5 / np.log(10):.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*W.T)
    # a fixed random subsample for the no-selection arm: the full file is ~2M stars and the
    # drift is a systematic, so a few 1e5 is ample and keeps the two arms comparable in cost
    sub = np.random.default_rng(999).permutation(len(o["d"]))[:nmax]

    out = {}
    for fr, tag in ((0.0, "none"), (FRAC, "smeared")):
        d = (observe.scatter_distance_modulus(o["d"], FRAC * 5 / np.log(10),
                                              np.random.default_rng(seed)) if fr > 0 else o["d"])
        w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"],
                                                      o["vlos"]))
        S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
        keep = np.random.default_rng(seed + 500_000).random(len(S)) < np.minimum(S, 1.0)
        out[f"sel_{tag}"] = phase(w[keep], "absvR", WIN) + BAR_ANGLE_DEG
        out[f"nosel_{tag}"] = phase(w[sub], "absvR", WIN) + BAR_ANGLE_DEG
        out[f"n_sel_{tag}"] = int(keep.sum())
    out["n_nosel"] = len(sub)
    return gname, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--n-real", type=int, default=6)
    ap.add_argument("--nmax", type=int, default=300_000)
    ap.add_argument("--snaps", nargs="+", default=AURIGA)
    args = ap.parse_args()

    jobs = [(g, 11 + s, args.nmax) for g in args.snaps for s in range(args.n_real)]
    print(f"{len(jobs)} runs over {len(args.snaps)} galaxies", flush=True)
    acc = {}
    with ProcessPoolExecutor(args.workers) as ex:
        for k, (g, ph) in enumerate(ex.map(one, jobs), 1):
            acc.setdefault(g, []).append(ph)
            if k % 10 == 0:
                print(f"  {k}/{len(jobs)}", flush=True)

    rows = []
    for g, lst in acc.items():
        m = {k: float(np.mean([x[k] for x in lst])) for k in lst[0]}
        rows.append(dict(galaxy=g,
                         sel_intrinsic=m["sel_none"], sel_smeared=m["sel_smeared"],
                         sel_drift=m["sel_smeared"] - m["sel_none"],
                         nosel_intrinsic=m["nosel_none"], nosel_smeared=m["nosel_smeared"],
                         nosel_drift=m["nosel_smeared"] - m["nosel_none"],
                         n_sel=m["n_sel_smeared"], n_nosel=m["n_nosel"]))
    D = pd.DataFrame(rows)
    Q = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    h = [int(g[2:]) for g in D.galaxy]
    D["R_CR"] = Q.loc[h, "R_CR"].values
    D = D.sort_values("R_CR").reset_index(drop=True)
    D.to_csv(RESULTS_DIR / "bar_angle_drift_no_selection.csv", index=False)

    print(f"\n12.5% distance-error drift, with and without the selection function "
          f"(deg, window {WIN})\n")
    print(f"  {'galaxy':8s} {'R_CR':>6s} | {'WITH selection':>21s} | {'NO selection':>21s}")
    print(f"  {'':8s} {'':6s} | {'intr':>6s} {'drift':>7s} {'N':>6s} | "
          f"{'intr':>6s} {'drift':>7s} {'N':>6s}")
    print("  " + "-" * 68)
    for _, r in D.iterrows():
        print(f"  {r.galaxy:8s} {r.R_CR:6.2f} | {r.sel_intrinsic:+6.2f} {r.sel_drift:+7.2f} "
              f"{r.n_sel:6.0f} | {r.nosel_intrinsic:+6.2f} {r.nosel_drift:+7.2f} {r.n_nosel:6.0f}")
    print("  " + "-" * 68)
    print(f"  {'spread (sd)':15s} | {D.sel_intrinsic.std():+6.2f} {D.sel_drift.std():+7.2f} "
          f"{'':6s} | {D.nosel_intrinsic.std():+6.2f} {D.nosel_drift.std():+7.2f}")
    print(f"  {'range':15s} | {'':6s} {D.sel_drift.max()-D.sel_drift.min():7.2f} "
          f"{'':6s} | {'':6s} {D.nosel_drift.max()-D.nosel_drift.min():7.2f}")
    print(f"  {'mean':15s} | {'':6s} {D.sel_drift.mean():+7.2f} "
          f"{'':6s} | {'':6s} {D.nosel_drift.mean():+7.2f}")

    print()
    for c in ("sel_drift", "nosel_drift"):
        r = stats.pearsonr(D.R_CR, D[c])
        s = stats.spearmanr(D.R_CR, D[c])
        print(f"  corr(R_CR, {c:12s}) = {r.statistic:+.2f} (p = {r.pvalue:.4f}), "
              f"rho = {s.statistic:+.2f}")
    r = stats.pearsonr(D.sel_drift, D.nosel_drift)
    print(f"  corr(sel_drift, nosel_drift)   = {r.statistic:+.2f} (p = {r.pvalue:.4f})")
    print(f"\nwrote {rel(RESULTS_DIR / 'bar_angle_drift_no_selection.csv')}")


if __name__ == "__main__":
    main()
