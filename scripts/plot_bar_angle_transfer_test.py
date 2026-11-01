"""Does the bar-angle calibration transfer between snapshots?

`plot_bar_angle_inversion.py` fits node = a * phi_true + b at one snapshot. That calibration is only useful
if it can be applied to a galaxy it was not derived from - otherwise it is a fit, not a method. Here the
relation fitted at the reference snapshot is applied to the measured nodes of a *different* snapshot, and the
inferred bar angle is compared with the truth.

The two snapshots are genuinely different bars: at t = 7 the bar is weaker and longer (A2 0.397, R_bar
5.0 kpc, Omega_b 25.2) than at t = 4 (A2 0.472, R_bar 4.25 kpc, Omega_b 33.3), so this is a real test of
transferability rather than a repeat of the same galaxy.

Needs both inversion tables; run plot_bar_angle_inversion.py for each snapshot first.
Usage: python plot_bar_angle_transfer_test.py [--ref nbody] [--test nbody_t4] [--frac 0.10]
Figure: figures/test2_bar_angle_transfer_<ref>_to_<test>.png
"""
import argparse
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, ACCENT, BLUE_RAMP
from galaxies import galaxy


def table(prefix, frac):
    p = RESULTS_DIR / f"{prefix}_bar_angle_inversion_{100*frac:g}pct.csv"
    if not p.exists():
        raise SystemExit(f"missing {rel(p)} - run plot_bar_angle_inversion.py --galaxy ... first")
    return pd.read_csv(p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="nbody", help="snapshot the calibration is fitted on")
    ap.add_argument("--test", default="nbody_t4", help="snapshot it is applied to")
    ap.add_argument("--frac", type=float, default=0.10)
    ap.add_argument("--col", default="node_sf", choices=["node_sf", "node_nosf"])
    args = ap.parse_args()
    set_style()
    gr, gt = galaxy(args.ref), galaxy(args.test)
    ref, tst = table(gr["prefix"], args.frac), table(gt["prefix"], args.frac)
    sd = args.col + "_sd"

    a, b = np.polyfit(ref.phi_true, ref[args.col], 1)          # calibration from the reference
    at, bt = np.polyfit(tst.phi_true, tst[args.col], 1)        # the test snapshot's own relation
    inferred = (tst[args.col] - b) / a                          # apply the reference calibration
    err = inferred - tst.phi_true
    sig = tst[sd] / abs(a)

    print(f"calibration fitted on {gr['label']}: node = {a:.3f} phi + {b:+.2f}")
    print(f"the {gt['label']} relation itself:   node = {at:.3f} phi + {bt:+.2f}")
    print(f"\n{'phi_true':>9s} {'node':>8s} {'inferred phi':>13s} {'error':>8s}")
    for _, r in tst.iterrows():
        print(f"{r.phi_true:9.1f} {r[args.col]:8.2f} {(r[args.col]-b)/a:13.2f} "
              f"{(r[args.col]-b)/a - r.phi_true:+8.2f}")
    print(f"\nmean error {err.mean():+.2f} deg, rms {np.sqrt((err**2).mean()):.2f} deg, "
          f"expected scatter {sig.median():.2f} deg")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.2, 5.2), gridspec_kw=dict(wspace=0.26))
    lim = (tst.phi_true.min() - 3, tst.phi_true.max() + 3)
    a1.plot(lim, lim, color=INK, lw=1.3, ls="--", label="perfect recovery")
    a1.errorbar(tst.phi_true, inferred, yerr=sig, fmt="o", color=ACCENT, ms=7, capsize=3, lw=1.4,
                label=f"{gt['label']} nodes, {gr['label']} calibration")
    a1.set_xlim(*lim); a1.set_ylim(*lim); a1.set_aspect("equal")
    a1.set_xlabel("true bar angle [deg]"); a1.set_ylabel("inferred bar angle [deg]")
    a1.legend(fontsize=8.5, loc="upper left")
    a1.set_title("transferring the calibration to another snapshot", color=INK_2, fontsize=10)

    a2.axhline(0, color=INK, lw=1.3, ls="--")
    a2.errorbar(tst.phi_true, err, yerr=sig, fmt="o", color=ACCENT, ms=7, capsize=3, lw=1.4)
    a2.axhspan(-sig.median(), sig.median(), color=BLUE_RAMP[0], alpha=0.5, lw=0,
               label=f"expected scatter ±{sig.median():.1f}°")
    a2.set_xlim(*lim)
    a2.set_xlabel("true bar angle [deg]"); a2.set_ylabel("inferred − true [deg]")
    a2.legend(fontsize=8.5)
    a2.set_title(f"residuals: mean {err.mean():+.2f}°, rms {np.sqrt((err**2).mean()):.2f}°",
                 color=INK_2, fontsize=10)
    fig.suptitle(f"Does the bar-angle calibration transfer? {gr['label']} → {gt['label']}, "
                 f"{args.frac:.0%} distance error + selection", fontsize=11.5, y=0.99)
    savefig(fig, f"test2_bar_angle_transfer_{gr['prefix']}_to_{gt['prefix']}")


if __name__ == "__main__":
    main()
