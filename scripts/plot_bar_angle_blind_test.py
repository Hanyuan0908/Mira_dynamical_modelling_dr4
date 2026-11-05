"""Blind test of the bar-angle calibration: build it on the N-body suite, apply it to Auriga 18.

This is the protocol a survey would follow. From a suite of simulated galaxies, each observed at the
survey's distance error with its selection function, build the relation between the true bar angle
and the measured m = 2 phase as a BAND: median and 16-84% over the galaxies. Measure the phase in
the real data, draw it as a horizontal line, and read the true angle off where the line crosses the
median, with its uncertainty from where it crosses the percentile curves.

Au18 is held out of the band entirely. It is a cosmological zoom with gas, a live dark halo and a
star-formation history, against an isolated two-component N-body disc, and its bar is structurally
different (A2 = 0.365 against 0.40-0.50, m = 2 density peak at 1.4 kpc against 2.4-3.4). If the
calibration recovers its bar angle, the relation is not specific to the simulation it came from.

Both estimators are shown, each with and without the bar-strength correction, which uses the m = 2
amplitude of the observed density measured from the same sample.

Usage: python plot_bar_angle_blind_test.py [--rule self:0.4]
Figure: figures/test2_bar_angle_blind_test.png;  table: results/bar_angle_blind_test.csv
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from scan_bar_angle_windows import phi_from, CACHE
from plot_bar_angle_amplitude_correction import observed_A2

EST = {"vR": r"$\langle v_R\rangle$", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}


def invert(ang, curve, obs):
    """True angle where a monotonic calibration curve crosses the observed value."""
    if curve[0] > curve[-1]:
        ang, curve = ang[::-1], curve[::-1]
    return float(np.interp(obs, curve, ang))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rule", default="self:0.4")
    args = ap.parse_args()
    set_style()
    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], list(d["snaps"]), list(d["angles"])
    ang = np.array(angles)
    # calibrate on the N-body sequence, test blind on the Auriga galaxies
    test = [i for i, g in enumerate(snaps) if g.startswith("au")]
    cal = [i for i, g in enumerate(snaps) if not g.startswith("au")]
    ia = test[0]

    A2 = np.array([observed_A2(g) for g in snaps])
    rows = []
    fig, AX = plt.subplots(2, 2, figsize=(12.8, 10.2))
    for col, est in enumerate(EST):
        ei = ["vR", "absvR"].index(est)
        P = np.array([[np.nanmean([phi_from(Z[gi, ai, s, ei], args.rule, Z[gi, ai, s, 0], est)
                                   for s in range(Z.shape[2])]) for ai in range(len(ang))]
                      for gi in range(len(snaps))])
        # bar-strength term fitted on the CALIBRATORS only, then applied to Au18
        x1 = np.repeat(ang[None, :], len(cal), 0).ravel()
        x2 = np.repeat(A2[cal][:, None] - A2[cal].mean(), len(ang), 1).ravel()
        c = np.linalg.lstsq(np.column_stack([x1, x2, np.ones_like(x1)]), P[cal].ravel(),
                            rcond=None)[0]
        Pc = P - c[1] * (A2[:, None] - A2[cal].mean())

        for row, (M, lab) in enumerate(((P, "no correction"),
                                        (Pc, "bar-strength corrected"))):
            ax = AX[row, col]
            B = M[cal]
            med, lo, hi = (np.nanmedian(B, 0), np.nanpercentile(B, 16, 0),
                           np.nanpercentile(B, 84, 0))
            ax.fill_between(ang, lo, hi, color=BLUE_RAMP[2], alpha=0.6, lw=0,
                            label=f"16–84% of the {len(cal)} N-body calibrators")
            ax.plot(ang, med, "-", color=BLUE_RAMP[6], lw=2.2, label="calibration median")
            for gi in cal:
                ax.plot(ang, M[gi], "-", color=INK_MUTED, lw=0.5, alpha=0.6, zorder=1)
            ax.plot(ang, ang, ":", color=INK_MUTED, lw=1.1, label="1:1")
            for j, ti in enumerate(test):
                ax.plot(ang, M[ti], "o-", color=ACCENT, lw=1.6, ms=5, zorder=6, alpha=0.85,
                        label="Auriga (held out)" if j == 0 else None)

            err = []
            for ti in test:
                for k, a in enumerate(ang):
                    t = invert(ang, med, M[ti, k])
                    tl, th = invert(ang, hi, M[ti, k]), invert(ang, lo, M[ti, k])
                    err.append(t - a)
                    rows.append(dict(galaxy=snaps[ti], est=est, corrected=row == 1, phi_true=a,
                                     phi_meas=M[ti, k], inferred=t, lo=min(tl, th),
                                     hi=max(tl, th), error=t - a))
            err = np.array(err)
            w = float(np.mean(hi - lo))
            sl = np.polyfit(ang, med, 1)[0]
            ax.set_title(f"{EST[est]} — {lab}\nband $\\pm${0.5*w/sl:.2f}$^\\circ$   |   "
                         f"**Auriga blind error: mean {err.mean():+.2f}$^\\circ$, "
                         f"rms {np.sqrt((err**2).mean()):.2f}$^\\circ$**",
                         color=INK, fontsize=10.5)
            ax.set_xlabel(r"true bar angle $\phi_b$ [deg]")
            ax.set_ylabel("measured bar angle [deg]")
            if row == 0 and col == 0:
                ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    fig.suptitle(f"Blind test — calibration built on {len(cal)} N-body snapshots, applied to "
                 f"{len(test)} independent Auriga galaxies (12.5% distance error)",
                 color=INK, y=0.995, fontsize=12)
    fig.tight_layout()
    savefig(fig, "test2_bar_angle_blind_test")

    R = pd.DataFrame(rows)
    R.to_csv(RESULTS_DIR / "bar_angle_blind_test.csv", index=False)
    print(f"\nAu18 blind inference (true angle recovered through the N-body calibration)\n")
    for est in EST:
        for corr in (False, True):
            g = R[(R.est == est) & (R.corrected == corr)]
            tag = "corrected" if corr else "raw      "
            print(f"  {est:6s} {tag}  mean error {g.error.mean():+6.2f}  "
                  f"rms {np.sqrt((g.error**2).mean()):5.2f}  max |err| {g.error.abs().max():5.2f} deg")
    print()
    g = R[(R.est == "absvR") & (R.corrected)]
    print(f"\n  per-galaxy blind error, absvR corrected:")
    for gal, sub in g.groupby("galaxy"):
        print(f"    {gal:8s} mean {sub.error.mean():+6.2f}  rms "
              f"{np.sqrt((sub.error**2).mean()):5.2f} deg")
    print(f"\nwrote {rel(RESULTS_DIR / 'bar_angle_blind_test.csv')}")


if __name__ == "__main__":
    main()
