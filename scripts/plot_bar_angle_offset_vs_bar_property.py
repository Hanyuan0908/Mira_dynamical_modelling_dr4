"""Evidence for the claim that the bar-angle calibration depends on a bar property.

The claim being tested: each galaxy's departure from the Auriga calibration is not random scatter
but tracks a property of its bar. The quantity plotted is

    offset_g = mean over the six true angles of [ phi_measured,g(angle) - median_Auriga(angle) ]

i.e. how far that galaxy sits from the calibration a survey would actually invert. It is plotted
against three candidate properties, with both suites shown, because the interesting result is not
only that the correlation exists but that the two suites follow DIFFERENT relations - which is
exactly why an R_CR correction fitted on Auriga makes the blind transfer to the N-body worse
(1.80 -> 2.40 deg rms, PROGRESS 2026-11-05).

The N-body properties come from results/nbody_pattern_speed_evolution.csv (measured per snapshot)
and the snapshot table; the Auriga ones from results/auriga_pattern_speed.csv and
results/auriga_bar_survey.csv.

Usage: python plot_bar_angle_offset_vs_bar_property.py [--rule self:0.3] [--est absvR]
Figure: figures/test2_bar_angle_offset_vs_bar_property.png
Table:  results/bar_angle_offset_vs_bar_property.csv
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from common import (RESULTS_DIR, set_style, savefig, rel, INK, INK_2, INK_MUTED, ACCENT,
                    BLUE_RAMP)
from nbody import SNAPSHOT_DIR
from scan_bar_angle_windows import phi_from, CACHE

PROPS = [("R_CR", r"corotation radius $R_{\rm CR}$ [kpc]"),
         ("Omega_p", r"pattern speed $\Omega_b$ [km s$^{-1}$ kpc$^{-1}$]"),
         ("A2_max", r"bar strength $A_2^{\rm max}$")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rule", default="self:0.3")
    ap.add_argument("--est", default="absvR", choices=["vR", "absvR"])
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    set_style()
    ei = ["vR", "absvR"].index(args.est)

    d = np.load(CACHE, allow_pickle=True)
    Z, snaps, angles = d["Z"], [str(x) for x in d["snaps"]], list(d["angles"])
    ang = np.array(angles)
    P = np.array([[np.nanmean([phi_from(Z[g, a, s, ei], args.rule, Z[g, a, s, 0], args.est)
                               for s in range(Z.shape[2])])
                   for a in range(len(ang))] for g in range(len(snaps))])
    iau = [i for i, g in enumerate(snaps) if g.startswith("au")]
    inb = [i for i, g in enumerate(snaps) if not g.startswith("au")]
    med = np.nanmedian(P[iau], 0)                 # the calibration a survey would invert
    offset = np.nanmean(P - med[None, :], axis=1)

    AQ = pd.read_csv(RESULTS_DIR / "auriga_pattern_speed.csv").set_index("halo")
    AS = pd.read_csv(RESULTS_DIR / "auriga_bar_survey.csv").set_index("halo")
    EV = pd.read_csv(RESULTS_DIR / "nbody_pattern_speed_evolution.csv")
    rows = []
    for i, g in enumerate(snaps):
        if g.startswith("au"):
            h = int(g[2:])
            rows.append(dict(galaxy=g, suite="Auriga", offset=offset[i],
                             R_CR=float(AQ.loc[h, "R_CR"]), Omega_p=float(AQ.loc[h, "Omega_p"]),
                             A2_max=float(AS.loc[h, "A2_max"])))
        else:
            t = 7.0 if g == "nbody" else float(g.split("_t")[1])
            k = (EV.t_model - t).abs().idxmin()
            rows.append(dict(galaxy=g, suite="N-body", offset=offset[i],
                             R_CR=float(EV.loc[k, "R_CR_kpc"]),
                             Omega_p=float(EV.loc[k, "Omega_p"]),
                             A2_max=float(EV.loc[k, "A2_max"])))
    D = pd.DataFrame(rows)
    D.to_csv(RESULTS_DIR / f"bar_angle_offset_vs_bar_property{args.tag}.csv", index=False)

    STY = {"Auriga": (BLUE_RAMP[6], "o", 9), "N-body": (ACCENT, "s", 7)}
    fig, AX = plt.subplots(1, 3, figsize=(15.2, 5.0))
    print(f"  {args.est} / {args.rule}:  offset from the Auriga median calibration\n")
    for ax, (col, lab) in zip(AX, PROPS):
        ax.axhline(0, color=INK_MUTED, lw=1.1, ls="--", zorder=1)
        txt = []
        for suite in ("Auriga", "N-body"):
            S = D[D.suite == suite]
            c, mk, ms = STY[suite]
            ax.plot(S[col], S.offset, mk, ms=ms, color=c, mec=INK, mew=0.8, zorder=4,
                    label=f"{suite} ({len(S)})")
            r = stats.pearsonr(S[col], S.offset)
            rho = stats.spearmanr(S[col], S.offset)
            b, a0 = np.polyfit(S[col], S.offset, 1)
            xs = np.linspace(S[col].min(), S[col].max(), 20)
            ax.plot(xs, b * xs + a0, "-", color=c, lw=1.7, alpha=0.9, zorder=3)
            txt.append(f"{suite}: $r$ = {r.statistic:+.2f} ($p$ = {r.pvalue:.3f}), "
                       f"$\\rho$ = {rho.statistic:+.2f}")
            print(f"    {lab.split('$')[0].strip():22s} {suite:7s} r = {r.statistic:+.2f} "
                  f"(p = {r.pvalue:.4f}), rho = {rho.statistic:+.2f}, "
                  f"slope = {b:+.3f} deg per unit")
        for i, g in enumerate(snaps):
            if g.startswith("au"):
                ax.annotate(g[2:], (D[col][i], D.offset[i]), textcoords="offset points",
                            xytext=(0, 9), ha="center", fontsize=7, color=INK_2)
        ax.set_xlabel(lab)
        ax.set_title("\n".join(txt), color=INK, fontsize=9)
    AX[0].set_ylabel("offset from the Auriga calibration [deg]")
    AX[0].legend(frameon=False, fontsize=9, loc="lower left")

    fig.suptitle(f"Does the calibration offset track a bar property?  "
                 f"{args.est}, window {args.rule}, 12.5% distance error — "
                 f"strong in $R_{{\\rm CR}}$, absent in bar strength, and the two suites "
                 f"follow different lines", color=INK, y=0.99, fontsize=11.5)
    fig.tight_layout()
    savefig(fig, f"test2_bar_angle_offset_vs_bar_property{args.tag}")

    # the suite-to-suite gap, which is what breaks the transfer
    a, n = D[D.suite == "Auriga"], D[D.suite == "N-body"]
    print(f"\n    mean offset: Auriga {a.offset.mean():+.2f} deg, N-body {n.offset.mean():+.2f} deg"
          f"  ->  gap {n.offset.mean() - a.offset.mean():+.2f} deg")
    ov = (min(a.R_CR.max(), n.R_CR.max()) - max(a.R_CR.min(), n.R_CR.min()))
    print(f"    R_CR ranges: Auriga {a.R_CR.min():.1f}-{a.R_CR.max():.1f}, "
          f"N-body {n.R_CR.min():.1f}-{n.R_CR.max():.1f} kpc (overlap {ov:.1f} kpc)")
    print(f"\nwrote {rel(RESULTS_DIR / f'bar_angle_offset_vs_bar_property{args.tag}.csv')}")


if __name__ == "__main__":
    main()
