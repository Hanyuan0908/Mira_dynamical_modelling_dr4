"""Pattern speed recovery at two distance-error levels, side by side.

Summary figure for the 12.5% campaign (2026-11-02/03): the same five selection-function mocks, the same
257-point focused Sobol design and the same chi2 scheme, run at sigma_DM = 0.217 (10%) and 0.271 (12.5%),
so the only difference between the two sets of points is the distance error of the mock.

  left   recovered against true Omega_b, both levels, with the 1:1 line
  right  fractional bias against true Omega_b, which is what the left panel cannot resolve

Scheme B (the dispersion ring averages dropped) and the sf library weighting throughout -- the combination
the 10% campaign settled on; scheme A trades chi2 against library temperature and reads ~8% high.
Truth is the simulation directory's finite-difference column, the convention adopted 2026-11-01.

Usage: python plot_test2_error_level_comparison.py [--truth tabulated_fd] [--tag NAME]
Figure: figures/test2_error_level_comparison<tag>.png
"""
import argparse
import json

import numpy as np
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from plot_test2_pattern_speed_recovery import SCHEME_B, tabulated

# (label, 10% run, 12.5% run, model time)
PAIRS = [("t = 2",   "test2_t2_sf_10pct",   "test2_t2_sf_12p5pct",   2.0),
         ("t = 3",   "test2_t3_sf_10pct",   "test2_t3_sf_12p5pct",   3.0),
         ("t = 4",   "test2_t4_sf_10pct",   "test2_t4_sf_12p5pct",   4.0),
         ("t = 5.5", "test2_t5.5_sf_10pct", "test2_t5.5_sf_12p5pct", 5.5),
         ("t = 7",   "test2_t7_sf_10pct",   "test2_t7_sf_12p5pct",   7.0)]
LEVELS = [("10%",   BLUE_RAMP[5], "o", 1), ("12.5%", ACCENT, "D", 2)]


def posterior(run, scheme, variant="sf"):
    """Omega_b posterior median and 16-84% range of one run, or None if it has not been emulated."""
    f = RESULTS_DIR / run / f"emulator_results_logLs_{scheme}_{variant}.json"
    if not f.exists():
        return None
    return json.load(open(f))["posterior"]["Omega_abs"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default="tabulated_fd", choices=["measured", "tabulated", "tabulated_fd"])
    ap.add_argument("--scheme", default=SCHEME_B)
    ap.add_argument("--variant", default="sf",
                    help="library weighting: sf (S known exactly), sfratio / sfemp (S estimated from "
                         "the data -- what a real survey can do), asis (no weighting)")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()

    rows = []
    for lbl, r10, r125, t_model in PAIRS:
        if args.truth == "measured":
            om = float(json.load(open(RESULTS_DIR / r10 / "run_meta.json"))["omega_truth"])
        else:
            om = tabulated(t_model, "pattern_speed_fd" if args.truth == "tabulated_fd"
                           else "pattern_speed_dehnen")
        p10, p125 = posterior(r10, args.scheme, args.variant), posterior(r125, args.scheme, args.variant)
        if p10 is None or p125 is None:
            print(f"  skipping {lbl}: missing {'10%' if p10 is None else '12.5%'} posterior")
            continue
        rows.append((lbl, om, p10, p125))

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.6, 5.0))
    lim = (22, 48)
    ax.plot(lim, lim, color=INK_MUTED, lw=1, ls="--", zorder=1)
    bx.axhline(0, color=INK_MUTED, lw=1, ls="--", zorder=1)

    for j, (name, colour, marker, _) in enumerate(LEVELS):
        xs, ys, los, his, bias, blo, bhi = [], [], [], [], [], [], []
        for lbl, om, p10, p125 in rows:
            p = (p10, p125)[j]
            xs.append(om); ys.append(p["p50"])
            los.append(p["p50"] - p["p16"]); his.append(p["p84"] - p["p50"])
            bias.append((p["p50"] - om) / om * 100)
            blo.append((p["p50"] - p["p16"]) / om * 100); bhi.append((p["p84"] - p["p50"]) / om * 100)
        kw = dict(color=colour, marker=marker, ms=7, lw=0, elinewidth=1.4, capsize=3, zorder=3 + j)
        ax.errorbar(xs, ys, yerr=[los, his], label=f"{name} distance error", **kw)
        bx.errorbar(xs, bias, yerr=[blo, bhi], label=f"{name}  (mean |bias| "
                    f"{np.mean(np.abs(bias)):.1f}%)", **kw)

    # join each pair, so the eye reads the shift at fixed snapshot rather than two clouds
    for lbl, om, p10, p125 in rows:
        ax.plot([om, om], [p10["p50"], p125["p50"]], color=INK_MUTED, lw=0.8, zorder=2)
        bx.plot([om, om], [(p10["p50"] - om) / om * 100, (p125["p50"] - om) / om * 100],
                color=INK_MUTED, lw=0.8, zorder=2)
        top = max((p125["p84"] - om) / om * 100, (p10["p84"] - om) / om * 100)
        bx.annotate(lbl, (om, top), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=9, color=INK_2)

    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel(r"true $\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]")
    ax.set_ylabel(r"recovered $\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]")
    ax.set_title("recovered against true", color=INK)
    ax.legend(frameon=False, loc="upper left", fontsize=9)

    bx.set_xlim(*lim)
    bx.set_xlabel(r"true $\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]")
    bx.set_ylabel(r"fractional bias [%]")
    bx.set_title("fractional bias", color=INK)
    bx.legend(frameon=False, loc="lower right", fontsize=9)

    VNAME = {"sf": "S known exactly (upper bound)", "sfratio": r"S estimated: $\rho_{\rm obs}/\rho_{\rm lib}$",
             "sfemp": "S estimated: separable fit", "asis": "no selection weighting"}
    fig.suptitle("Pattern speed recovery at 10% and 12.5% distance error  (scheme B, "
                 + VNAME.get(args.variant, args.variant) + ")", color=INK, y=0.99)
    fig.tight_layout()
    savefig(fig, f"test2_error_level_comparison{args.tag}")

    print(f"\n{'snapshot':9s} {'truth':>7s} {'10%':>8s} {'bias':>7s} {'12.5%':>8s} {'bias':>7s} {'change':>8s}")
    print("-" * 60)
    b10, b125 = [], []
    for lbl, om, p10, p125 in rows:
        f1, f2 = (p10["p50"] - om) / om * 100, (p125["p50"] - om) / om * 100
        b10.append(f1); b125.append(f2)
        print(f"{lbl:9s} {om:7.2f} {p10['p50']:8.2f} {f1:+6.1f}% {p125['p50']:8.2f} {f2:+6.1f}% "
              f"{f2 - f1:+7.1f}%")
    print("-" * 60)
    print(f"{'mean |bias|':9s} {'':7s} {'':8s} {np.mean(np.abs(b10)):6.1f}% {'':8s} "
          f"{np.mean(np.abs(b125)):6.1f}%")
    w10 = np.mean([(p["p84"] - p["p16"]) / 2 for _, _, p, _ in rows])
    w125 = np.mean([(p["p84"] - p["p16"]) / 2 for _, _, _, p in rows])
    print(f"mean 68% half-width: {w10:.2f} -> {w125:.2f} km/s/kpc (x{w125 / w10:.2f})")


if __name__ == "__main__":
    main()
