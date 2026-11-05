"""The bar-angle calibration as a band, and how to invert it for a real measurement.

This is the product a survey actually uses. Across a suite of simulated galaxies, each observed at
a 12.5% distance error with a DR3-Mira-like selection, the measured m = 2 phase is a compressed,
biased version of the true bar angle. The spread between galaxies - not the bias - is the
irreducible part, so the calibration is a BAND, not a line:

  - at each true bar angle, the median and the 16-84% range of the measured angle over all galaxies;
  - a horizontal line at the angle measured in the real data;
  - the true angle is where that line crosses the median, and its uncertainty is where it crosses
    the 16th and 84th percentile curves.

The band is narrowed by correcting for bar strength, which the same data provides: the m = 2
amplitude of the observed density. Both versions are drawn so the gain is visible.

Usage: python plot_bar_angle_calibration_band.py [--est absvR] [--rule self:0.4] [--observed 20.9]
Figure: figures/test2_bar_angle_calibration_band.png
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from scan_bar_angle_windows import phi_from, CACHE
from plot_bar_angle_amplitude_correction import observed_A2


def invert(ang, med, lo, hi, obs):
    """True angle and its 16-84 range, by crossing the observed value with the band."""
    f = lambda y: float(np.interp(obs, y, ang)) if y[0] < y[-1] else float(np.interp(obs, y[::-1], ang[::-1]))
    return f(med), f(hi), f(lo)      # a higher curve is crossed at a SMALLER true angle


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--est", default="absvR")
    ap.add_argument("--rule", default="self:0.4")
    ap.add_argument("--observed", type=float, default=None,
                    help="measured bar angle to invert (default: the band median at 27 deg)")
    args = ap.parse_args()
    set_style()
    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], list(d["snaps"]), list(d["angles"])
    ei = ["vR", "absvR"].index(args.est)
    ang = np.array(angles)
    ns, na = len(snaps), len(angles)

    P = np.array([[np.nanmean([phi_from(Z[gi, ai, s, ei], args.rule, Z[gi, ai, s, 0], args.est)
                               for s in range(Z.shape[2])]) for ai in range(na)]
                  for gi in range(ns)])
    A2 = np.array([observed_A2(g) for g in snaps])

    # bar-strength corrected version: remove the fitted dependence on the observed m = 2 amplitude
    x1 = np.repeat(ang[None, :], ns, 0).ravel()
    x2 = np.repeat(A2[:, None] - A2.mean(), na, 1).ravel()
    c = np.linalg.lstsq(np.column_stack([x1, x2, np.ones_like(x1)]), P.ravel(), rcond=None)[0]
    Pc = P - c[1] * (A2[:, None] - A2.mean())

    fig, AX = plt.subplots(1, 2, figsize=(12.6, 5.4), sharey=True)
    for ax, M, lab in ((AX[0], P, "no correction"),
                       (AX[1], Pc, "corrected for the observed $m=2$ amplitude")):
        med = np.nanmedian(M, 0)
        lo = np.nanpercentile(M, 16, 0)
        hi = np.nanpercentile(M, 84, 0)
        ax.fill_between(ang, lo, hi, color=BLUE_RAMP[2], alpha=0.55, lw=0,
                        label="16–84% over the galaxies")
        ax.plot(ang, med, "-", color=BLUE_RAMP[6], lw=2.2, label="median")
        for gi in range(ns):
            ax.plot(ang, M[gi], "-", color=INK_MUTED, lw=0.6, alpha=0.7, zorder=1)
        ax.plot(ang, ang, ":", color=INK_MUTED, lw=1.1, label="1:1")

        obs = args.observed if args.observed is not None else float(np.interp(27.0, ang, med))
        t, tlo, thi = invert(ang, med, lo, hi, obs)
        ax.axhline(obs, color=ACCENT, lw=1.6, ls="--")
        ax.plot([t], [obs], "o", color=ACCENT, ms=9, zorder=6)
        ax.annotate(f"observed {obs:.1f}$^\\circ$", (ang[0] + 0.4, obs), xytext=(0, 5),
                    textcoords="offset points", fontsize=9, color=ACCENT)
        ax.plot([tlo, thi], [obs, obs], color=ACCENT, lw=5, alpha=0.35, solid_capstyle="butt")
        ax.annotate(f"$\\phi_b = {t:.1f}^{{+{thi-t:.1f}}}_{{-{t-tlo:.1f}}}$ deg",
                    (t, obs), xytext=(6, -22), textcoords="offset points", fontsize=11,
                    color=INK, fontweight="bold")
        w = float(np.mean(hi - lo))
        sl = np.polyfit(ang, med, 1)[0]
        ax.set_title(f"{lab}\nband width {w:.2f}$^\\circ$, slope {sl:.2f}  "
                     f"$\\Rightarrow$ $\\pm${0.5*w/sl:.2f}$^\\circ$ on $\\phi_b$",
                     color=INK, fontsize=10.5)
        ax.set_xlabel(r"true bar angle $\phi_b$ [deg]")
        ax.set_xlim(ang[0] - 1, ang[-1] + 1)
    AX[0].set_ylabel("measured bar angle [deg]")
    AX[0].legend(frameon=False, fontsize=9, loc="upper left")
    fig.suptitle(f"Bar-angle calibration band at a 12.5% distance error — {args.est}, "
                 f"window {args.rule}, {ns} galaxies", color=INK, y=0.98)
    fig.tight_layout()
    savefig(fig, "test2_bar_angle_calibration_band")

    for M, lab in ((P, "no correction"), (Pc, "corrected")):
        med = np.nanmedian(M, 0)
        w = float(np.mean(np.nanpercentile(M, 84, 0) - np.nanpercentile(M, 16, 0)))
        sl = np.polyfit(ang, med, 1)[0]
        print(f"  {lab:12s}: band width {w:5.2f} deg, slope {sl:.3f} -> "
              f"phi_b uncertainty +-{0.5*w/sl:.2f} deg")


if __name__ == "__main__":
    main()
