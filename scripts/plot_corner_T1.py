"""Corner plots of a run's posterior at T = 1 (logL = −Δχ²/2), from the chains of rerun_emulators_T1.py.

Two figures:
  <run>_corner_T1.png              the T = 1 posterior on its own scale, with 16/50/84% titles; prior box edges that
                                   fall inside a panel are drawn as dashed lines, since a posterior piled against one
                                   is cut by the prior, not by the data
  <run>_corner_T1_vs_tempered.png  T = 1 over the √(2N) tempered posterior on the full prior box
The truth is marked when the run has one.

Usage: python plot_corner_T1.py --run osarg_zhang2024_mw_hunter24
"""
import argparse
import json

import numpy as np
import corner
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from common import RESULTS_DIR, set_style, savefig, INK, INK_MUTED, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    args = ap.parse_args()
    set_style()
    ge.use_run(args.run)
    ch = np.load(RESULTS_DIR / args.run / "posterior_chains_T1.npz")
    res = json.load(open(RESULTS_DIR / args.run / "emulator_results_T1.json"))
    box = np.column_stack([ge.LO, ge.HI])
    box[2:] = np.exp(box[2:])                     # the emulator box holds ln S and ln bar_frac
    truth = res["truth"]
    tv = [truth[k] for k in ge.KEYS] if truth is not None else None

    post = ch["T1"].astype(float)
    lo, hi = np.percentile(post, [0.05, 99.95], axis=0)
    zoom = [(a - 0.25 * (b - a), b + 0.25 * (b - a)) for a, b in zip(lo, hi)]   # padded past any box edge it touches
    fig = corner.corner(post, labels=ge.NAMES, range=zoom, truths=tv, truth_color=ACCENT, color=BLUE_RAMP[5],
                        quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".2f", bins=30,
                        title_kwargs=dict(fontsize=9), label_kwargs=dict(fontsize=10), plot_datapoints=False)
    AX = np.array(fig.axes).reshape(4, 4)
    for i in range(4):
        for j in range(i + 1):
            ax = AX[i, j]
            for e, side in zip(box[j], (-1, 1)):
                x0, x1 = ax.get_xlim()
                if x0 <= e <= x1:
                    ax.axvline(e, ls="--", color=INK, lw=1.2)
                    ax.axvspan(*((x0, e) if side < 0 else (e, x1)), color=INK_MUTED, alpha=0.25, lw=0)
            if i != j:
                for e, side in zip(box[i], (-1, 1)):
                    y0, y1 = ax.get_ylim()
                    if y0 <= e <= y1:
                        ax.axhline(e, ls="--", color=INK, lw=1.2)
                        ax.axhspan(*((y0, e) if side < 0 else (e, y1)), color=INK_MUTED, alpha=0.25, lw=0)
    sc = res["T1"]
    fig.suptitle(f"{args.run}: posterior at T = 1 (logL = −Δχ²/2), emulator refitted on {res['n_train']} design "
                 f"points;\ndashed line and grey shading = prior box edge and the region outside it; emulator leave-one-out rms {sc['loo_rms']:.1f} in logL"
                 + ("; orange = truth" if tv else ""), y=1.03, fontsize=10, color=INK)
    savefig(fig, f"{args.run}_corner_T1")

    rng = [tuple(b) for b in box]
    temp = ch["tempered"].astype(float)
    fig = corner.corner(temp, labels=ge.NAMES, range=rng, color=INK_MUTED, plot_datapoints=False, bins=30,
                        smooth=1.5, label_kwargs=dict(fontsize=10), truths=tv, truth_color=ACCENT)
    corner.corner(post, range=rng, color=BLUE_RAMP[5], plot_datapoints=False, bins=60, smooth=1.0, fig=fig)
    fig.legend(handles=[Line2D([], [], color=INK_MUTED, lw=2, label=f"tempered, √(2N) = {res['tempered']['scale']:.1f}"),
                        Line2D([], [], color=BLUE_RAMP[5], lw=2, label="T = 1")],
               loc="upper right", bbox_to_anchor=(0.97, 0.95), frameon=False, fontsize=11)
    fig.suptitle(f"{args.run}: T = 1 against the √(2N) tempered posterior, on the full prior box", y=1.01,
                 fontsize=10, color=INK)
    savefig(fig, f"{args.run}_corner_T1_vs_tempered")
    for tag in ("T1", "tempered"):
        p = res[tag]["posterior"]
        print(f"[{tag}] " + ", ".join(f"{k} {p[k]['p50']:.3f} ({p[k]['p16']:.3f}-{p[k]['p84']:.3f})" for k in ge.KEYS))
    print("prior box:", {k: tuple(np.round(b, 3)) for k, b in zip(ge.KEYS, box)})


if __name__ == "__main__":
    main()
