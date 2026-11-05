"""Calibrate the bar-angle relation on the 13 Auriga galaxies; use the N-body run as a blind check.

The measured m = 2 phase is a compressed, biased version of the true bar angle, and the compression
depends on the galaxy - on 2026-11-04 the per-galaxy offset was shown to track corotation radius
(Spearman -0.84) and NOT bar strength (r = -0.09), with no mechanism yet identified (three were
tested and refuted). Since the dependence is on some bar property we cannot yet name, the honest
calibration is a BAND over a sample that spans a wide range of bar properties, not a line fitted to
one galaxy.

That argues for calibrating on Auriga rather than on the N-body run:
  - 13 separate galaxies formed in cosmological environments, against 12 snapshots of one secularly
    evolving isolated bar whose properties are all correlated with each other through time;
  - a much wider spread of pattern speed (16.8-65.1 against 26-44 km/s/kpc), corotation radius
    (4.8-11.3 kpc), bar length and bar strength.

So: the band is the 16-84% range of the measured angle over the 13 Auriga at each true angle, and
the 12 N-body snapshots are then inferred THROUGH that band and compared with their known truth.
Nothing about the N-body enters the calibration, which is what makes it a blind check.

Four panels: {<v_R>, <|v_R|/v_tot>} x {fixed 1-4 kpc, kinematic window}, reading
results/bar_angle_transfer_12p5.csv (written by run_transfer_parallel.py).

Usage: python plot_bar_angle_auriga_calibration.py [--drop au5]
Figure: figures/test2_bar_angle_auriga_calibration.png
Table:  results/bar_angle_auriga_calibration.csv
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP

EST = {"vR": r"$\langle v_R\rangle$", "absvR": r"$\langle |v_R|/v_{\rm tot}\rangle$"}
WINS = {"fixed": "fixed 1–4 kpc", "kin": "kinematic window"}


def invert(ang, curve, measured):
    """True angle implied by a measured angle, read off a monotonic calibration curve."""
    if curve[0] > curve[-1]:
        return float(np.interp(measured, curve[::-1], ang[::-1]))
    return float(np.interp(measured, curve, ang))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drop", nargs="*", default=[],
                    help="Auriga galaxies to exclude from the calibration, e.g. --drop au5")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()

    D = pd.read_csv(RESULTS_DIR / "bar_angle_transfer_12p5.csv")
    ang = np.array(sorted(D.phi_true.unique()))
    au = sorted({g for g in D.galaxy if str(g).startswith("au")} - set(args.drop))
    nb = sorted({g for g in D.galaxy if not str(g).startswith("au")})
    print(f"calibration: {len(au)} Auriga galaxies" + (f" (dropped {', '.join(args.drop)})" if args.drop else ""))
    print(f"blind check: {len(nb)} N-body snapshots\n")

    def series(g, col):
        d = D[D.galaxy == g].sort_values("phi_true")
        return d[col].values

    fig, AX = plt.subplots(2, 2, figsize=(11.6, 9.0))
    rows = []
    for i, e in enumerate(EST):
        for j, w in enumerate(WINS):
            col = f"{e}_{w}"
            ax = AX[i, j]
            M = np.array([series(g, col) for g in au])          # (n_auriga, n_angle)
            med = np.nanmedian(M, 0)
            lo = np.nanpercentile(M, 16, 0)
            hi = np.nanpercentile(M, 84, 0)

            # blind check: infer every N-body snapshot through the Auriga median
            err = []
            for g in nb:
                meas = series(g, col)
                err.extend([invert(ang, med, m) - a for m, a in zip(meas, ang)])
            err = np.array(err)

            width = float(np.mean(hi - lo))
            slope = float(np.polyfit(ang, med, 1)[0])
            dphi = 0.5 * width / slope
            rows.append(dict(est=e, window=w, band_width=width, slope=slope, phi_err=dphi,
                             blind_mean=err.mean(), blind_rms=float(np.sqrt((err ** 2).mean())),
                             blind_max=float(np.abs(err).max())))

            ax.fill_between(ang, lo, hi, color=BLUE_RAMP[2], alpha=0.55, lw=0,
                            label=f"Auriga 16–84% ({len(au)} galaxies)")
            ax.plot(ang, med, "-", color=BLUE_RAMP[6], lw=2.4, zorder=5, label="Auriga median")
            for k, g in enumerate(au):
                ax.plot(ang, M[k], "-", color=INK_MUTED, lw=0.6, alpha=0.65, zorder=2,
                        label="individual Auriga" if k == 0 else None)
            for k, g in enumerate(nb):
                ax.plot(ang, series(g, col), "-", color=ACCENT, lw=0.9, alpha=0.85, zorder=3,
                        label="N-body (blind check)" if k == 0 else None)
            ax.plot(ang, ang, ":", color=INK_MUTED, lw=1.1, zorder=1, label="1:1")

            ax.set_title(f"{EST[e]} / {WINS[w]}\n"
                         f"band {width:.2f}$^\\circ$, slope {slope:.2f} "
                         f"$\\Rightarrow\\pm${dphi:.2f}$^\\circ$ on $\\phi_b$  |  "
                         f"blind N-body rms {rows[-1]['blind_rms']:.2f}$^\\circ$",
                         color=INK, fontsize=10)
            ax.set_xlabel(r"true bar angle $\phi_b$ [deg]")
            ax.set_ylabel("measured bar angle [deg]")
            ax.set_xlim(ang[0] - 1, ang[-1] + 1)
            if i == 0 and j == 0:
                ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.suptitle("Bar-angle calibration from the 13 Auriga galaxies, N-body as a blind check "
                 "(12.5% distance error)", color=INK, y=0.995, fontsize=12)
    fig.tight_layout()
    savefig(fig, f"test2_bar_angle_auriga_calibration{args.tag}")

    T = pd.DataFrame(rows).sort_values("blind_rms")
    T.to_csv(RESULTS_DIR / f"bar_angle_auriga_calibration{args.tag}.csv", index=False)
    print(f"  {'estimator / window':34s} {'band':>6s} {'slope':>6s} {'+-phi_b':>8s} | "
          f"{'blind mean':>10s} {'rms':>6s} {'max':>6s}")
    print("  " + "-" * 82)
    for _, r in T.iterrows():
        print(f"  {EST[r.est] + ' / ' + WINS[r.window]:34s} {r.band_width:6.2f} {r.slope:6.2f} "
              f"{r.phi_err:8.2f} | {r.blind_mean:+10.2f} {r.blind_rms:6.2f} {r.blind_max:6.2f}")
    b = T.iloc[0]
    print(f"\n  best blind-check: {EST[b.est]} / {WINS[b.window]} -> rms {b.blind_rms:.2f} deg, "
          f"calibration uncertainty +-{b.phi_err:.2f} deg")
    print(f"wrote {rel(RESULTS_DIR / f'bar_angle_auriga_calibration{args.tag}.csv')}")


if __name__ == "__main__":
    main()
