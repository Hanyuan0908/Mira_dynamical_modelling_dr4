"""Does pooling K distance-error draws per library sample improve the emulator? (yes)

The orbit library is a Monte Carlo estimate of the model's predicted moments, and the runs make it
with ONE distance-error draw per sample. chi2 is quadratic in the coefficient error, so a one-sample
library biases chi2 HIGH by an amount set by the local coefficient variance - which varies over
parameter space, so it does not cancel in a posterior. Pooling K draws into a single fit cuts that
variance by K and removes the bias; averaging K separate chi2 values would not.

This compares the full 257-point design of one run scored at K = 1 and K = 10.

Usage: python plot_pooled_draws_test.py [--run test2_t7_sf_12p5pct]
Figure: figures/emulator_validation/pooled_draws_<run>.png
"""
import argparse
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import ROOT, RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from nbody import SNAPSHOT_DIR

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
# measured at three design points x 5 seeds, scripts in PROGRESS 2026-11-03
SCALING = {1: (0.65, 0.00), 2: (0.44, 1.09), 5: (0.23, 1.66), 10: (0.19, 1.84)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default="test2_t7_sf_12p5pct")
    ap.add_argument("--time", type=float, default=7.0)
    args = ap.parse_args()
    set_style()
    ev = pd.read_csv(RESULTS_DIR / args.run / "evals.csv")
    a, b = ev[f"logLs_{B}_sfratio"], ev[f"logLs_{B}_sfratioK10"]
    g = np.isfinite(a) & np.isfinite(b)
    A, Bv = a[g].values, b[g].values
    d = Bv - A
    J = {k: json.load(open(RESULTS_DIR / args.run / f"emulator_results_logLs_{B}_{c}.json"))
         for k, c in (("K1", "sfratio"), ("K10", "sfratioK10"))}
    tab = pd.read_csv(SNAPSHOT_DIR / "bar_strength_and_pattern_speed.csv")
    T = float(tab[np.isclose(tab.snap_time, args.time)]["pattern_speed_fd"].iloc[0])

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.6))

    # 1. the measured scaling
    ax = axes[0]
    K = np.array(sorted(SCALING)); sd = np.array([SCALING[k][0] for k in K])
    ax.plot(K, sd, "o-", color=ACCENT, ms=7, lw=1.6, label="measured scatter")
    ax.plot(K, sd[0] / np.sqrt(K), "--", color=INK_MUTED, lw=1.4, label=r"$1/\sqrt{K}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xticks(K); ax.set_xticklabels(K)
    ax.set_xlabel("K draws per library sample")
    ax.set_ylabel("seed-to-seed scatter in logL")
    ax.set_title("the noise falls as $1/\\sqrt{K}$", color=INK)
    ax.legend(frameon=False, fontsize=9)

    # 2. the shift across the whole design
    ax = axes[1]
    ax.plot(A, d, "o", ms=3.5, color=BLUE_RAMP[3])
    ax.axhline(d.mean(), color=ACCENT, lw=1.6, label=f"mean {d.mean():+.2f}")
    ax.axhline(0, color=INK_MUTED, lw=1, ls="--")
    ax.fill_between([A.min(), A.max()], d.mean() - d.std(), d.mean() + d.std(),
                    color=ACCENT, alpha=0.15, lw=0, label=f"sd {d.std():.2f}")
    ax.set_xlim(A.max() - 60, A.max() + 2)
    ax.set_xlabel("logL at K = 1"); ax.set_ylabel("logL(K=10) $-$ logL(K=1)")
    ax.set_title("chi$^2$ falls, by a point-dependent amount", color=INK)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    # 3. what it does to the emulator and the answer
    ax = axes[2]
    lab = ["LOO rms\n(emulator error)", r"$|$bias$|$ in $\Omega_b$" + "\n(% of truth)"]
    k1 = [J["K1"]["loo_rms"], abs(J["K1"]["posterior"]["Omega_abs"]["p50"] - T) / T * 100]
    k10 = [J["K10"]["loo_rms"], abs(J["K10"]["posterior"]["Omega_abs"]["p50"] - T) / T * 100]
    x = np.arange(2); w = 0.36
    ax.bar(x - w / 2, k1, w, color=BLUE_RAMP[4], label="K = 1 (the runs as made)")
    ax.bar(x + w / 2, k10, w, color=ACCENT, label="K = 10 (pooled draws)")
    for xi, v in zip(x - w / 2, k1):
        ax.text(xi, v + 0.04, f"{v:.2f}", ha="center", fontsize=9, color=INK_2)
    for xi, v in zip(x + w / 2, k10):
        ax.text(xi, v + 0.04, f"{v:.2f}", ha="center", fontsize=9, color="#b4480f")
    ax.axhline(0.5, color=INK_MUTED, ls=":", lw=1.2)
    ax.text(1.42, 0.55, r"1$\sigma$ of the posterior", fontsize=8, color=INK_2, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(lab)
    ax.set_ylabel("value")
    ax.set_title("the emulator improves 2.3x", color=INK)
    ax.legend(frameon=False, fontsize=9)

    p1, p10 = J["K1"]["posterior"]["Omega_abs"], J["K10"]["posterior"]["Omega_abs"]
    fig.suptitle(f"Pooling K distance-error draws per library sample — {args.run}, scheme B, "
                 f"S estimated by density ratio\n"
                 f"$\\Omega_b$ = {p1['p50']:.2f} (K=1) → {p10['p50']:.2f} (K=10), truth {T:.2f}",
                 color=INK, y=1.04, fontsize=11)
    fig.tight_layout()
    savefig(fig, f"emulator_validation/pooled_draws_{args.run}")

    print(f"\n{'':8s} {'LOO rms':>8s} {'Omega_b':>22s} {'bias':>7s}")
    for k, lb in (("K1", "K = 1"), ("K10", "K = 10")):
        p = J[k]["posterior"]["Omega_abs"]
        print(f"{lb:8s} {J[k]['loo_rms']:8.2f} {p['p50']:7.2f} [{p['p16']:.2f}-{p['p84']:.2f}] "
              f"{(p['p50']-T)/T*100:+6.1f}%")
    print(f"\nlogL shift: mean {d.mean():+.2f}, sd {d.std():.2f} over {len(d)} design points")


if __name__ == "__main__":
    main()
