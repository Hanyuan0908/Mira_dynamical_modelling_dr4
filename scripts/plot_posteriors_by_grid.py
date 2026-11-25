"""All four marginal posteriors for every Fourier grid, one row per grid, with a drawing of that grid's design.

For a run, reads the posterior samples that compare_grid_strategies.py saved (results/<run>/grid_strategies_chains.npz;
each grid at its own T1) and its summary (grid_strategies.csv). Each row:
  - left: the grid design in the (R, |z|) plane, the cells drawn to scale. Blue cells are fitted (rings inside
    6 kpc); the grey 0.5 kpc rings at 6-8 kpc exist in the grid but are not in the fitted scheme. Stars above the top
    layer are not used. Under it: N_data, T1, best chi2/N.
  - right: |Omega_b|, phi_b, S (bar length scale) and bar_frac (bar amplitude). Light histogram, the 16-84% interval
    darker, the median as a line, and its value written above. The columns share their x and y axes, so a narrower
    posterior has a taller peak. For a mock, the truth is a dashed line and the |Omega_b| bias is given.
The prior box is the plotted x-range of every column (flat in |Omega_b| and phi_b; flat in ln S and ln bar_frac).

Usage: python plot_posteriors_by_grid.py --run osarg_zhang2024_sormani22
       python plot_posteriors_by_grid.py --run mw_hunter24_mock_from_osarg_ics --truth 32.5 25 1 1
Figure: figures/<run>_all_posteriors_by_grid_design.png
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import test1_nbody_qmc_focused_fourier_stars as fs

# (compare_grid_strategies label, set_grid arguments, row title, design note)
GRIDS = [
    ("default", (0.3, 3.0, 1.0, None), "1. Current default: 0.3 kpc rings, 3 layers",
     "20 rings, each 0.3 kpc wide (R 0-6 kpc)\n3 |z| layers: 0-1, 1-2, 2-3 kpc"),
    ("_r0.5_z2", (0.5, 2.0, 1.0, None), "2. 0.5 kpc rings, 2 layers",
     "12 rings, each 0.5 kpc wide\n2 |z| layers: 0-1, 1-2 kpc"),
    ("_r0.5_z2_dz2", (0.5, 2.0, 2.0, None), "3. 0.5 kpc rings, one layer |z| < 2",
     "12 rings, each 0.5 kpc wide\n1 |z| layer: 0-2 kpc"),
    ("_r1_z2_dz2", (1.0, 2.0, 2.0, None), "4. 1 kpc rings, one layer |z| < 2",
     "6 rings, each 1 kpc wide\n1 |z| layer: 0-2 kpc"),
    ("_r0.5_z3_dz3", (0.5, 3.0, 3.0, None), "5. 0.5 kpc rings, one layer |z| < 3",
     "12 rings, each 0.5 kpc wide\n1 |z| layer: 0-3 kpc"),
    ("_r0.3_z2_dz2", (0.3, 2.0, 2.0, None), "6. 0.3 kpc rings, one layer |z| < 2",
     "20 rings, each 0.3 kpc wide\n1 |z| layer: 0-2 kpc"),
    ("_geo20_z2_dz2", (0.3, 2.0, 2.0, "geo20"), "7. Geometric, 20 rings, one layer |z| < 2",
     "ring 0-0.2 kpc, then 19 rings to 6 kpc,\neach 1.196x wider than the one inside.\n"
     "Width 0.05 kpc at R 0.25, 0.28 at 1.6,\n0.58 at 3, 0.82 at 5 (last: 5.02-6.0)\n1 |z| layer: 0-2 kpc"),
    ("_geo32_z2_dz2", (0.3, 2.0, 2.0, "geo32"), "8. Geometric, 32 rings, one layer |z| < 2",
     "ring 0-0.2 kpc, then 31 rings to 6 kpc,\neach 1.116x wider than the one inside.\n"
     "Width 0.03 kpc at R 0.25, 0.17 at 1.6,\n0.32 at 3, 0.56 at 5 (last: 5.38-6.0)\n1 |z| layer: 0-2 kpc"),
    ("_ilr01_z2_dz2", (0.3, 2.0, 2.0, "ilr01"), "9. Refined at the Milky Way ILR, one layer |z| < 2",
     "28 rings: 0.3 kpc wide over 0-1.2 kpc (4),\n0.1 kpc over 1.2-2.4 kpc (12),\n"
     "0.3 kpc over 2.4-6 kpc (12)\n1 |z| layer: 0-2 kpc"),
]
PARAMS = [(r"$|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]", (15, 45), 1.0, "{:.1f}"),
          (r"$\phi_b$ [deg]", (20, 30), 0.4, "{:.1f}"),
          (r"$S$ (bar length scale)", (0.5, 1.5), 0.025, "{:.2f}"),
          (r"bar_frac (bar amplitude)", (0.5, 1.5), 0.025, "{:.2f}")]


def draw_design(ax, edges, layers):
    for z0, z1 in zip(layers[:-1], layers[1:]):
        for r0, r1 in zip(edges[:-1], edges[1:]):
            fitted = r1 <= 6.0 + 1e-9
            ax.add_patch(Rectangle((r0, z0), r1 - r0, z1 - z0, facecolor=BLUE_RAMP[1] if fitted else "#e6e5e1",
                                   edgecolor=BLUE_RAMP[5] if fitted else INK_MUTED, lw=0.5))
    ax.set(xlim=(0, 8), ylim=(0, 3), xticks=[0, 2, 4, 6, 8], yticks=[0, 1, 2, 3])
    ax.set_xlabel("R [kpc]", fontsize=8, labelpad=1)
    ax.set_ylabel("|z| [kpc]", fontsize=8, labelpad=1)
    ax.tick_params(labelsize=7.5, length=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--truth", type=float, nargs=4, default=None, metavar=("OMEGA", "PHI", "S", "F"))
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    set_style()
    chains = np.load(RESULTS_DIR / args.run / "grid_strategies_chains.npz")
    summ = pd.read_csv(RESULTS_DIR / args.run / "grid_strategies.csv").set_index("tag")
    rows = [g for g in GRIDS if g[0] in chains.files]

    H = 2.15 * len(rows) + 1.9
    fig = plt.figure(figsize=(18.5, H))
    gs = fig.add_gridspec(len(rows), 5, width_ratios=[1.8, 1, 1, 1, 1], hspace=0.95, wspace=0.12,
                          left=0.035, right=0.99, top=1 - 1.75 / H, bottom=0.035)
    # common y-scale per column
    ymax = np.zeros(4)
    hists = {}
    for lab, *_ in rows:
        for k, (_, rng, dx, _) in enumerate(PARAMS):
            bins = np.arange(rng[0], rng[1] + dx / 2, dx)
            h, _ = np.histogram(chains[lab][:, k], bins=bins, density=True)
            hists[lab, k] = (h, bins)
            ymax[k] = max(ymax[k], h.max())

    for i, (lab, gargs, title, note) in enumerate(rows):
        fs.set_grid(*gargs)
        s = summ.loc[lab]
        ax = fig.add_subplot(gs[i, 0])
        draw_design(ax, fs.RING_EDGES, fs.LAYER_EDGES)
        ax.set_title(title, loc="left", fontsize=10, color=INK, fontweight="bold", pad=4)
        ax.text(1.04, 1.0, note, transform=ax.transAxes, fontsize=8, color=INK_2, va="top", ha="left",
                linespacing=1.2)
        ax.text(1.04, 0.0, f"N_data {int(s.N_data)}   T1 {s.T1:.1f}   best χ²/N {s.best_chi2_per_N:.3f}",
                transform=ax.transAxes, fontsize=8.5, color=INK, va="bottom", ha="left")
        # make room for the text on the right of the drawing
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * 0.42, pos.height])
        ax.tick_params(top=False, right=False)
        for k, (plab, rng, dx, fmt) in enumerate(PARAMS):
            a = fig.add_subplot(gs[i, k + 1])
            h, bins = hists[lab, k]
            x = chains[lab][:, k]
            q16, q50, q84 = np.percentile(x, [16, 50, 84])
            c = 0.5 * (bins[:-1] + bins[1:])
            inside = (c >= q16) & (c <= q84)
            a.bar(c, h, width=dx, color=BLUE_RAMP[1], lw=0)
            a.bar(c[inside], h[inside], width=dx, color=BLUE_RAMP[4], lw=0)
            a.axvline(q50, color=INK, lw=1.4)
            if args.truth:
                a.axvline(args.truth[k], color=ACCENT, lw=1.4, ls="--")
            txt = f"{fmt.format(q50)}  ({fmt.format(q16)}–{fmt.format(q84)})"
            if k == 0:
                txt += f"\nwidth {100 * (q84 - q16) / q50:.0f}%"
                if args.truth:
                    txt += f", bias {100 * (q50 / args.truth[0] - 1):+.1f}%"
            a.set_title(txt, fontsize=9, color=INK, pad=3, loc="left")
            a.set(xlim=rng, ylim=(0, 1.08 * ymax[k]), yticks=[])
            a.tick_params(labelsize=7.5, length=2, top=False, right=False, left=False)
            if i == len(rows) - 1:
                a.set_xlabel(plab, fontsize=10)
            for side in ("top", "right", "left"):
                a.spines[side].set_visible(False)
    fs.set_grid()

    # column headers
    for k, (plab, *_) in enumerate(PARAMS):
        ref = fig.axes[1 + k]
        bb = ref.get_position()
        fig.text(bb.x0, 1 - 0.9 / H, plab, fontsize=12, color=INK, ha="left", va="top", fontweight="bold")
    fig.text(0.035, 1 - 0.9 / H, "Grid design in (R, |z|), to scale", fontsize=12, color=INK, ha="left", va="top",
             fontweight="bold")
    title = args.title or args.run
    sub = ("Each grid at its own temperature T1. Dark band: 16–84%; black line: median" +
           ("; orange dashed: truth" if args.truth else "") +
           ". Grid drawings: blue cells are fitted (R < 6 kpc); grey rings at 6–8 kpc are not fitted; "
           "stars above the top layer are not used.")
    fig.text(0.035, 1 - 0.12 / H, title, fontsize=14, color=INK, ha="left", va="top", fontweight="bold")
    fig.text(0.035, 1 - 0.45 / H, sub, fontsize=9, color=INK_2, ha="left", va="top")
    savefig(fig, f"{args.run}_all_posteriors_by_grid_design")


if __name__ == "__main__":
    main()
