"""Remove the bar-strength dependence from the kinematic bar-angle calibration.

The measured m = 2 phase is compressed towards the Sun-GC line by distance smearing, and the amount
of compression depends on the bar: across thirteen galaxies the residual offset from a single
calibration correlates with bar strength A2 at r = -0.73 (p = 0.007). A calibration in phi_true
alone therefore carries a bar-dependent systematic, which is what the colour ordering by snapshot
time shows.

The fix is to calibrate a SURFACE rather than a line, using a strength proxy measured from the same
data. The proxy used here is the amplitude of the observed kinematic m = 2 signal itself - no
photometry, no prior knowledge of the bar:

    phi_obs = a phi_true + b + c (A_kin - <A_kin>)

fitted on the simulations and inverted for a target galaxy using its own A_kin. Scored by
leave-one-out: fit on twelve galaxies, infer the bar angle of the thirteenth.

Usage: python plot_bar_angle_amplitude_correction.py [--est absvR] [--rule self:0.4]
Figure: figures/test2_bar_angle_amplitude_correction.png
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT
from scan_bar_angle_windows import phi_from, CACHE


def observed_A2(gname, frac=0.125, seed=7):
    """m = 2 amplitude of the OBSERVED density over 1-4 kpc: the bar-strength proxy.

    Measured on the same smeared, selection-applied sample the angle comes from, so it needs no
    photometry and no prior knowledge of the bar. It tracks the true A2 at r = +0.86 and recovers
    almost all of the gain that knowing the true A2 would give."""
    import h5py
    from common import MOCK_DIR
    from galaxies import galaxy
    import observe, selection, steady_state_likelihood as ssl
    g = galaxy(gname)
    tab = selection.load_tables(
        MOCK_DIR / f"{g['prefix']}_mock_selection_sigmaDM{frac * 5 / np.log(10):.3f}_N50k.hdf5")
    with h5py.File(g["stars"], "r") as f:
        W = np.column_stack([f[k][:].astype(np.float64)
                             for k in ("x", "y", "z", "vx", "vy", "vz")])
    o = observe.to_galactic(*W.T)
    d = observe.scatter_distance_modulus(o["d"], frac * 5 / np.log(10),
                                         np.random.default_rng(seed))
    w = np.column_stack(observe.to_galactocentric(o["l"], o["b"], d, o["pml"], o["pmb"], o["vlos"]))
    S = selection.evaluate(tab, o["l"], o["b"], d, np.hypot(w[:, 0], w[:, 1]))
    sel = w[np.random.default_rng(seed + 70).random(len(S)) < np.minimum(S, 1.0)]
    R = np.hypot(sel[:, 0], sel[:, 1])
    phi = np.arctan2(sel[:, 1], sel[:, 0])
    m = (R > 1) & (R < 4) & (np.abs(sel[:, 2]) < 1.5)
    return float(np.abs(np.sum(np.exp(2j * phi[m]))) / m.sum())


def amp_of(z, rule, zref):
    """Mean kinematic m = 2 amplitude over the rings the rule selects (kept for reference)."""
    from scan_bar_angle_windows import Rc
    ok = np.isfinite(z)
    kind, _, arg = rule.partition(":")
    if kind == "fixed":
        lo, hi = (float(v) for v in arg.split("-"))
        m = ok & (Rc > lo) & (Rc < hi)
    elif kind in ("self", "vR"):
        a = np.abs(z if kind == "self" else zref)
        m = ok & (a > float(arg) * np.nanmax(a))
    else:
        rp = Rc[int(np.nanargmax(np.abs(z)))]
        lo, hi = (float(v) for v in arg.split("-"))
        m = ok & (Rc > lo * rp) & (Rc < hi * rp)
    return float(np.nanmean(np.abs(z[m]))) if m.any() else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--est", default="absvR")
    ap.add_argument("--rule", default="self:0.4")
    args = ap.parse_args()
    set_style()
    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], list(d["snaps"]), list(d["angles"])
    ei = ["vR", "absvR"].index(args.est)
    ns, na = len(snaps), len(angles)

    P = np.full((ns, na), np.nan)
    A = np.full((ns, na), np.nan)
    for gi in range(ns):
        for ai in range(na):
            P[gi, ai] = np.nanmean([phi_from(Z[gi, ai, s, ei], args.rule, Z[gi, ai, s, 0], args.est)
                                    for s in range(Z.shape[2])])
    Ag = np.array([observed_A2(g) for g in snaps])   # bar-strength proxy, one per galaxy
    ang = np.array(angles)

    def loo(use_amp):
        """Leave-one-out: fit on all but one galaxy, infer the held-out galaxy's bar angle."""
        errs = []
        for gi in range(ns):
            tr = [j for j in range(ns) if j != gi]
            x1 = np.repeat(ang[None, :], len(tr), 0).ravel()
            y = P[tr].ravel()
            if use_amp:
                x2 = np.repeat(Ag[tr][:, None] - Ag.mean(), na, 1).ravel()
                M = np.column_stack([x1, x2, np.ones_like(x1)])
            else:
                M = np.column_stack([x1, np.ones_like(x1)])
            c = np.linalg.lstsq(M, y, rcond=None)[0]
            if use_amp:
                inferred = (P[gi] - c[1] * (Ag[gi] - Ag.mean()) - c[2]) / c[0]
            else:
                inferred = (P[gi] - c[1]) / c[0]
            errs.append(inferred - ang)
        return np.array(errs)

    e0, e1 = loo(False), loo(True)
    r0 = np.sqrt(np.nanmean(e0 ** 2))
    r1 = np.sqrt(np.nanmean(e1 ** 2))
    off = np.nanmean(P - np.nanmean(P, 0), axis=1)
    rr = stats.pearsonr(Ag, off)
    print(f"{args.est} / {args.rule}")
    print(f"  offset vs kinematic amplitude: r = {rr.statistic:+.2f} (p = {rr.pvalue:.4f})")
    print(f"  leave-one-out systematic: {r0:.2f} deg  ->  {r1:.2f} deg with the amplitude term")
    print(f"  worst galaxy:             {np.sqrt(np.nanmean(e0**2,1)).max():.2f}  ->  "
          f"{np.sqrt(np.nanmean(e1**2,1)).max():.2f} deg")

    # viridis along the N-body time sequence, a warm ramp across the Auriga galaxies; Au18 black.
    # (This used to be `[cmap(i/(ns-2)) for i in range(ns-1)] + [INK]`, which assumed the last
    # galaxy in the list was Au18 - true only while Au18 was the only non-N-body galaxy.)
    isau = [str(g).startswith("au") for g in snaps]
    cmap, amap = plt.get_cmap("viridis"), plt.get_cmap("autumn")
    nnb, nau = isau.count(False), isau.count(True)
    inb = iau = 0
    tcol = []
    for g, a in zip(snaps, isau):
        if a:
            tcol.append(INK if g == "au18" else amap(0.82 * iau / max(nau - 1, 1)))
            iau += 1
        else:
            tcol.append(cmap(inb / max(nnb - 1, 1)))
            inb += 1
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.9))

    a = ax[0]
    a.axhline(0, color=INK_MUTED, lw=1, ls="--")
    for gi, g in enumerate(snaps):
        a.plot(Ag[gi], off[gi], "s" if isau[gi] else "o", ms=10 if g == "au18" else 7,
               color=tcol[gi], mec=INK if g == "au18" else "none", mew=1.5)
        a.annotate(f"Au{str(g)[2:]}" if isau[gi]
                   else g.replace("nbody_t", "").replace("nbody", "7"),
                   (Ag[gi], off[gi]), textcoords="offset points", xytext=(0, 9), ha="center",
                   fontsize=8, color=INK_2)
    p = np.polyfit(Ag, off, 1)
    xs = np.linspace(Ag.min(), Ag.max(), 20)
    a.plot(xs, np.polyval(p, xs), color=ACCENT, lw=1.6)
    a.set_xlabel("observed density $m=2$ amplitude  (bar-strength proxy)")
    a.set_ylabel("offset from the mean relation [deg]")
    a.set_title(f"the residual tracks bar strength\nr = {rr.statistic:+.2f} "
                f"(p = {rr.pvalue:.3f})", color=INK, fontsize=10.5)

    for k, (E, lab, rms) in enumerate(((e0, "calibrated on $\\phi_b$ only", r0),
                                       (e1, "+ bar-strength term", r1))):
        b = ax[1 + k]
        b.axhline(0, color=INK_MUTED, lw=1, ls="--")
        for gi, g in enumerate(snaps):
            b.plot(ang, E[gi], "-s" if isau[gi] else "-o", ms=7 if g == "au18" else 4,
                   lw=2.4 if g == "au18" else 1.1, color=tcol[gi],
                   zorder=6 if g == "au18" else 2,
                   label="Au18" if g == "au18" else None)
        b.set_xlabel(r"true bar angle $\phi_b$ [deg]")
        b.set_ylabel("inferred $-$ true [deg]")
        b.set_ylim(-11, 11)
        b.set_title(f"{lab}\nleave-one-out rms {rms:.2f}$^\\circ$", color=INK, fontsize=10.5)
        if k == 1:
            b.legend(frameon=False, fontsize=9)
    fig.suptitle(f"Removing the bar-strength dependence — {args.est}, window {args.rule}, "
                 f"12.5% distance error, {nnb} N-body snapshots (circles, viridis = time) + "
                 f"{nau} Auriga galaxies (squares, black = Au18)",
                 color=INK, y=0.99, fontsize=11)
    fig.tight_layout()
    savefig(fig, "test2_bar_angle_amplitude_correction")


if __name__ == "__main__":
    main()
