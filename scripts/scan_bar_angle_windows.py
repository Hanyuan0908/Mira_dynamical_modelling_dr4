"""Which radial window makes the bar-angle calibration transfer best?

The measured m = 2 phase depends on which rings are summed. A window fixed in kpc samples a
different fraction of the bar at each epoch, because the bar grows (density m = 2 peaks at 2.4 kpc
at t = 2 and 3.4 kpc at t = 7), so it leaves a time-ordered residual in the calibration. A window
defined from the data's own kinematic amplitude should remove that - but <|v_R|/v_tot| has a much
flatter radial profile than <v_R>, so a threshold localises it poorly and the residual may survive.

This caches the per-ring m = 2 coefficients once, for every (galaxy, bar angle, realisation), then
scans window rules over the cache at no further cost. Rules tested:

  fixed:L-H       rings between L and H kpc
  self:F          rings whose OWN m = 2 amplitude exceeds F of its peak (self-consistent)
  vR:F            rings selected by the <v_R> amplitude profile, for BOTH estimators. <v_R> is a
                  streaming signal that switches off beyond the bar, so it localises the bar
                  sharply; <|v_R|/v_tot> tracks sigma_R/v_tot, which the disc has everywhere.
  peakfrac:L-H    rings between L and H times the radius where the OWN amplitude peaks

Scored by the leave-one-out transfer rms: calibrate on each galaxy in turn, infer the bar angle of
all the others, and take the rms error - the systematic a user incurs applying a simulation
calibration to a galaxy they have not calibrated on.

Usage: python scan_bar_angle_windows.py [--rebuild] [--n-real 8]
Cache: results/bar_angle_ring_cache.npz;  table: results/bar_angle_window_scan.csv
"""
import argparse

import numpy as np
import pandas as pd
import h5py

from common import MOCK_DIR, RESULTS_DIR, BAR_ANGLE_DEG, rel
from galaxies import galaxy
import observe, selection, steady_state_likelihood as ssl
from bar_angle_vecsum_verification import ring_coeffs, Rc

NBODY = ["nbody", "nbody_t2", "nbody_t2.375", "nbody_t2.625", "nbody_t3", "nbody_t3.375",
         "nbody_t3.625", "nbody_t4", "nbody_t4.5", "nbody_t5", "nbody_t5.5", "nbody_t6.25"]
# Barred Auriga galaxies with a well-defined bar. Au2 and Au13 carry a companion above 5% of the
# host stellar mass inside R200 (check_auriga_interactions) and are excluded: an interaction breaks
# the steady state and drives its own m = 2.
AURIGA = ["au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24",
          "au26", "au27", "au28"]
SNAPS = NBODY + AURIGA
ANGLES = [15.0, 20.0, 25.0, 27.0, 30.0, 35.0]
FRAC = 0.125
EST = ["vR", "absvR"]
CACHE = RESULTS_DIR / "bar_angle_ring_cache.npz"


def build(snaps, angles, n_real):
    """Per-ring m = 2 coefficients for every (galaxy, angle, realisation, estimator)."""
    out = np.full((len(snaps), len(angles), n_real, len(EST), len(Rc)), np.nan, complex)
    for gi, gn in enumerate(snaps):
        g = galaxy(gn)
        tab = selection.load_tables(
            MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{FRAC * 5 / np.log(10):.3f}_N50k.hdf5")
        with h5py.File(g["stars"], "r") as f:
            W = np.column_stack([f[k][:].astype(np.float64)
                                 for k in ("x", "y", "z", "vx", "vy", "vz")])
        for ai, a in enumerate(angles):
            o = observe.to_galactic(*ssl.rotate_z(W, np.radians(BAR_ANGLE_DEG - a)).T)
            seed0 = 20_000 + 7919 * ai + 104_729 * gi
            for s in range(n_real):
                d = observe.scatter_distance_modulus(o["d"], FRAC * 5 / np.log(10),
                                                     np.random.default_rng(seed0 + s))
                w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"],
                                                              o["pmb"], o["vlos"]))
                S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
                sel = w[np.random.default_rng(seed0 + 500_000 + s).random(len(S)) < np.minimum(S, 1)]
                for ei, e in enumerate(EST):
                    a2, b2 = ring_coeffs(sel, e)
                    out[gi, ai, s, ei] = a2 + 1j * b2
            print(f"  {gn:14s} phi {a:5.1f} done", flush=True)
    return out


def angle_from_coeffs(z, est):
    """Measured bar angle [deg from the Sun-GC line] from a summed complex m = 2 coefficient.

    THE one place this conversion lives. Three pieces, all of which have been got wrong at least
    once: the <v_R> extrema sit 45 deg from its node, so `off` removes that; ring_coeffs rotates
    its input by +BAR_ANGLE_DEG, so that must come off before wrapping, not after; and the result
    is MINUS the bar angle, because the bar lies at negative azimuth in this frame. Check:
    a perfect measurement of a bar at 27 deg returns +27.
    """
    off = 45.0 if est == "vR" else 0.0
    return -float((0.5 * np.degrees(np.angle(z)) - off - BAR_ANGLE_DEG + 90) % 180 - 90)


def phi_from(z, rule, zref, est="absvR"):
    """Measured bar angle from per-ring coefficients z, under one window rule.

    zref is the <v_R> coefficient set, used by the 'vR:' rules to pick the rings."""
    ok = np.isfinite(z)
    kind, _, arg = rule.partition(":")
    if kind == "fixed":
        lo, hi = (float(v) for v in arg.split("-"))
        m = ok & (Rc > lo) & (Rc < hi)
    elif kind in ("self", "vR"):
        amp = np.abs(z if kind == "self" else zref)
        amp = np.where(np.isfinite(amp), amp, np.nan)
        m = ok & (amp > float(arg) * np.nanmax(amp))
    elif kind == "peakfrac":
        lo, hi = (float(v) for v in arg.split("-"))
        rp = Rc[int(np.nanargmax(np.abs(z)))]
        m = ok & (Rc > lo * rp) & (Rc < hi * rp)
    else:
        raise ValueError(rule)
    if not m.any():
        return np.nan
    return angle_from_coeffs(np.nansum(z[m]), est)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--n-real", type=int, default=8)
    args = ap.parse_args()

    if args.rebuild or not CACHE.exists():
        Z = build(SNAPS, ANGLES, args.n_real)
        np.savez_compressed(CACHE, Z=Z, snaps=SNAPS, angles=ANGLES, est=EST, Rc=Rc)
        print(f"wrote {rel(CACHE)}")
    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], list(d["snaps"]), list(d["angles"])

    RULES = (["fixed:1-4", "fixed:1-3", "fixed:0.5-3"]
             + [f"self:{f:g}" for f in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
             + [f"vR:{f:g}" for f in (0.5, 0.6, 0.7, 0.8, 0.9)]
             + ["peakfrac:0.3-1.2", "peakfrac:0.4-1.5", "peakfrac:0.5-2.0"])

    rows = []
    for ei, e in enumerate(EST):
        for rule in RULES:
            P = np.full((len(snaps), len(angles)), np.nan)
            for gi in range(len(snaps)):
                for ai in range(len(angles)):
                    v = [phi_from(Z[gi, ai, s, ei], rule, Z[gi, ai, s, 0], e)
                         for s in range(Z.shape[2])]
                    P[gi, ai] = np.nanmean(v)
            loo = []
            for ci in range(len(snaps)):
                A = np.polyfit(angles, P[ci], 1)
                err = [(P[gi] - A[1]) / A[0] - np.array(angles)
                       for gi in range(len(snaps)) if gi != ci]
                loo.append(np.sqrt(np.nanmean(np.concatenate(err) ** 2)))
            loo = np.array(loo)
            slope = np.polyfit(angles, P.mean(0), 1)[0]
            # Time-ordered residual: correlation of each galaxy's mean offset with snapshot TIME,
            # over the N-body sequence only. Two things were wrong before: [:-1] assumed the only
            # non-N-body galaxy was au18 and sat last, which stopped being true once the Auriga
            # haloes were added; and np.arange indexed the list, whose first entry is t = 7, not
            # the earliest snapshot. Only the N-body run has a time ordering at all - the Auriga
            # galaxies are 13 separate z = 0 galaxies and cannot enter this statistic.
            off = np.nanmean(P - P.mean(0), axis=1)
            nb = [i for i, g in enumerate(snaps) if not str(g).startswith("au")]
            tt = [7.0 if snaps[i] == "nbody" else float(str(snaps[i]).split("_t")[1]) for i in nb]
            tord = np.corrcoef(tt, off[nb])[0, 1] if len(nb) > 2 else np.nan
            rows.append(dict(est=e, rule=rule, slope=slope, loo_mean=loo.mean(),
                             loo_worst=loo.max(), t7_rms=loo[0], time_corr=tord))
    R = pd.DataFrame(rows).sort_values(["est", "loo_mean"])
    R.to_csv(RESULTS_DIR / "bar_angle_window_scan.csv", index=False)
    for e in EST:
        print(f"\n=== {e} ===   (loo_mean = systematic; |time_corr| near 0 = no residual trend)")
        print(f"  {'rule':18s} {'slope':>7s} {'loo mean':>9s} {'loo worst':>10s} "
              f"{'on t=7':>7s} {'time corr':>10s}")
        for _, r in R[R.est == e].iterrows():
            print(f"  {r['rule']:18s} {r.slope:7.3f} {r.loo_mean:9.2f} {r.loo_worst:10.2f} "
                  f"{r.t7_rms:7.2f} {r.time_corr:+10.2f}")
    b = R.loc[R.loo_mean.idxmin()]
    print(f"\nbest overall: {b['est']} / {b['rule']} -> systematic {b.loo_mean:.2f} deg "
          f"(slope {b.slope:.2f}, time corr {b.time_corr:+.2f})")
    print(f"wrote {rel(RESULTS_DIR / 'bar_angle_window_scan.csv')}")


if __name__ == "__main__":
    main()
