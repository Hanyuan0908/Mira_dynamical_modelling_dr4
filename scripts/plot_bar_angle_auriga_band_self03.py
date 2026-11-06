"""Measured against true bar angle for the 13 Auriga galaxies, with the 16-84% band.

This is the calibration a survey would actually invert: at each true bar angle, the median measured
angle over the 13 Auriga galaxies and the range the middle 68% of them span. The spread between
galaxies - not the mean bias - is the irreducible part, because the compression depends on a bar
property that is still unidentified (PROGRESS 2026-11-04/05).

Both estimators are shown on the same axes scale so the difference in slope is visible: <v_R>
compresses hard (slope ~0.4, the measured angle barely responds to the truth) while
<|v_R|/v_tot> stays near 1:1.

Window: self:0.3, adopted 2026-11-05 - rings whose own m = 2 amplitude exceeds 0.3 of its peak.
It gives the best blind transfer to the held-out N-body set (1.80 deg rms against 4.70 for
self:0.7) and is the only rule whose formal band width and blind error agree, so its error bar
can be quoted as it stands.

Usage: python plot_bar_angle_auriga_band_self03.py [--rule kin0.3] [--show-nbody]
Figure: figures/test2_bar_angle_auriga_band.png
Table:  results/bar_angle_auriga_band.csv
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from scan_bar_angle_windows import phi_from, CACHE

EST = {"vR": r"$\langle v_R\rangle$", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rule", default="self:0.3")
    ap.add_argument("--show-nbody", action="store_true",
                    help="overlay the 12 N-body snapshots, which are not in the calibration")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()

    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], [str(x) for x in d["snaps"]], list(d["angles"])
    ang = np.array(angles)
    iau = [i for i, g in enumerate(snaps) if g.startswith("au")]
    inb = [i for i, g in enumerate(snaps) if not g.startswith("au")]

    fig, AX = plt.subplots(1, 2, figsize=(12.0, 5.6), sharex=True, sharey=True)
    rows = []
    for k, (e, elab) in enumerate(EST.items()):
        ax = AX[k]
        P = np.array([[np.nanmean([phi_from(Z[g, a, s, k], args.rule, Z[g, a, s, 0], e)
                                   for s in range(Z.shape[2])])
                       for a in range(len(ang))] for g in range(len(snaps))])
        A = P[iau]
        med = np.nanmedian(A, 0)
        lo, hi = np.nanpercentile(A, 16, 0), np.nanpercentile(A, 84, 0)
        width = float(np.mean(hi - lo))
        slope = float(np.polyfit(ang, med, 1)[0])
        for a, m, l, h in zip(ang, med, lo, hi):
            rows.append(dict(est=e, rule=args.rule, phi_true=a, p16=l, median=m, p84=h))

        ax.plot(ang, ang, ":", color=INK_MUTED, lw=1.2, zorder=1, label="1:1 (no bias)")
        if args.show_nbody:
            for i in inb:
                ax.plot(ang, P[i], "-", color=ACCENT, lw=0.8, alpha=0.7, zorder=2,
                        label="N-body (not in the calibration)" if i == inb[0] else None)
        for i in iau:
            ax.plot(ang, P[i], "-", color=INK_MUTED, lw=0.7, alpha=0.6, zorder=3,
                    label="individual Auriga" if i == iau[0] else None)
        ax.fill_between(ang, lo, hi, color=BLUE_RAMP[2], alpha=0.6, lw=0, zorder=4,
                        label="16–84% over the 13 galaxies")
        ax.plot(ang, med, "-o", color=BLUE_RAMP[6], lw=2.6, ms=6, zorder=6, label="median")

        ax.set_title(f"{elab}\nband {width:.2f}$^\\circ$, slope {slope:.2f}  "
                     f"$\\Rightarrow\\ \\pm${0.5 * width / slope:.2f}$^\\circ$ on $\\phi_b$",
                     color=INK, fontsize=11)
        ax.set_xlabel(r"true bar angle $\phi_b$ [deg]")
        ax.set_xlim(ang[0] - 1.5, ang[-1] + 1.5)
        if k == 0:
            ax.set_ylabel("measured bar angle [deg]")
            ax.legend(frameon=False, fontsize=8.5, loc="upper left")
        print(f"  {e:6s} band {width:5.2f}  slope {slope:5.2f}  "
              f"+-phi_b {0.5 * width / slope:5.2f}")

    fig.suptitle(f"Bar-angle calibration from 13 Auriga galaxies at a 12.5% distance error "
                 f"— window {args.rule}", color=INK, y=0.98, fontsize=12)
    fig.tight_layout()
    savefig(fig, f"test2_bar_angle_auriga_band{args.tag}")
    T = pd.DataFrame(rows)
    T.to_csv(RESULTS_DIR / f"bar_angle_auriga_band{args.tag}.csv", index=False)
    print(f"\n  {'est':7s} {'phi_true':>8s} {'p16':>7s} {'median':>7s} {'p84':>7s}")
    for _, r in T.iterrows():
        print(f"  {r.est:7s} {r.phi_true:8.1f} {r.p16:7.2f} {r['median']:7.2f} {r.p84:7.2f}")
    print(f"\nwrote {rel(RESULTS_DIR / f'bar_angle_auriga_band{args.tag}.csv')}")


if __name__ == "__main__":
    main()
