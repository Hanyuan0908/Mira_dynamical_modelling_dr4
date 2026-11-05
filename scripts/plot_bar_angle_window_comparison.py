"""The kinematic-window panels of the Auriga calibration, remade for several window rules.

Same construction as plot_bar_angle_auriga_calibration.py - band from the 13 Auriga, N-body
inferred through it as a blind check - but reading the per-ring coefficient cache instead of the
transfer table, so any window rule can be drawn without re-measuring anything.

Each panel reports the two numbers that matter and the one that decides between them:
  +-phi_b     the formal calibration uncertainty, half the band width divided by the slope
  blind rms   the actual error on the 12 held-out N-body snapshots, which INCLUDES the
              suite-to-suite offset that the band does not see
  rho(R_CR)   Spearman of the per-Auriga-galaxy offset against corotation radius: the
              bar-property dependence itself. This is what we are trying to make small.

Usage: python plot_bar_angle_window_comparison.py [--rules self:0.7 self:0.3 peakfrac:0.2-1]
Figure: figures/test2_bar_angle_window_comparison.png
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from scan_bar_angle_windows import phi_from, CACHE

EST = {"vR": r"$\langle v_R\rangle$", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}


def label(rule):
    kind, _, arg = rule.partition(":")
    if kind == "self":
        return f"rings above {arg} of the peak $m=2$ amplitude"
    if kind == "vR":
        return rf"rings above {arg} of the $\langle v_R\rangle$ peak"
    if kind == "peakfrac":
        lo, hi = arg.split("-")
        return rf"rings at {lo}–{hi} $\times R_{{\rm peak}}$"
    return f"fixed {arg} kpc"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rules", nargs="+", default=["self:0.7", "self:0.3", "peakfrac:0.2-1"])
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()

    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], [str(x) for x in d["snaps"]], list(d["angles"])
    ang = np.array(angles)
    iau = [i for i, g in enumerate(snaps) if g.startswith("au")]
    inb = [i for i, g in enumerate(snaps) if not g.startswith("au")]
    Q = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    rcr = np.array([float(Q.loc[int(snaps[i][2:]), "R_CR"]) for i in iau])

    def measured(ei, rule, est):
        return np.array([[np.nanmean([phi_from(Z[g, a, s, ei], rule, Z[g, a, s, 0], est)
                                      for s in range(Z.shape[2])])
                          for a in range(len(ang))] for g in range(len(snaps))])

    def invert(curve, m):
        return (np.interp(m, curve[::-1], ang[::-1]) if curve[0] > curve[-1]
                else np.interp(m, curve, ang))

    nr, nc = len(EST), len(args.rules)
    fig, AX = plt.subplots(nr, nc, figsize=(4.6 * nc, 4.6 * nr), squeeze=False)
    for i, (e, elab) in enumerate(EST.items()):
        for j, rule in enumerate(args.rules):
            ax = AX[i][j]
            P = measured(i, rule, e)
            A = P[iau]
            med = np.nanmedian(A, 0)
            lo, hi = np.nanpercentile(A, 16, 0), np.nanpercentile(A, 84, 0)
            width = float(np.mean(hi - lo))
            slope = float(np.polyfit(ang, med, 1)[0])
            err = np.array([invert(med, P[k]) - ang for k in inb]).ravel()
            rho = float(stats.spearmanr(rcr, np.nanmean(A - A.mean(0), axis=1)).statistic)

            ax.fill_between(ang, lo, hi, color=BLUE_RAMP[2], alpha=0.55, lw=0,
                            label="Auriga 16–84% (13)")
            ax.plot(ang, med, "-", color=BLUE_RAMP[6], lw=2.4, zorder=5, label="Auriga median")
            for k in iau:
                ax.plot(ang, P[k], "-", color=INK_MUTED, lw=0.6, alpha=0.6, zorder=2,
                        label="individual Auriga" if k == iau[0] else None)
            for k in inb:
                ax.plot(ang, P[k], "-", color=ACCENT, lw=0.9, alpha=0.85, zorder=3,
                        label="N-body (blind)" if k == inb[0] else None)
            ax.plot(ang, ang, ":", color=INK_MUTED, lw=1.1, zorder=1, label="1:1")
            ax.set_title(f"{elab}   {rule}\n{label(rule)}\n"
                         f"$\\pm${0.5 * width / slope:.2f}$^\\circ$ formal   "
                         f"blind rms {np.sqrt((err ** 2).mean()):.2f}$^\\circ$   "
                         f"$\\rho(R_{{\\rm CR}})$ {rho:+.2f}", color=INK, fontsize=9.5)
            ax.set_xlabel(r"true $\phi_b$ [deg]")
            if j == 0:
                ax.set_ylabel("measured bar angle [deg]")
            ax.set_xlim(ang[0] - 1, ang[-1] + 1)
            if i == 0 and j == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="upper left")
            print(f"  {e:6s} {rule:16s} band {width:5.2f}  slope {slope:5.2f}  "
                  f"+-phi {0.5 * width / slope:5.2f}  blind rms {np.sqrt((err**2).mean()):5.2f}  "
                  f"rho(R_CR) {rho:+.2f}")

    fig.suptitle("Kinematic window and the bar-property dependence — band from 13 Auriga, "
                 "N-body held out as a blind check (12.5% distance error)",
                 color=INK, y=0.995, fontsize=11.5)
    fig.tight_layout()
    savefig(fig, f"test2_bar_angle_window_comparison{args.tag}")


if __name__ == "__main__":
    main()
