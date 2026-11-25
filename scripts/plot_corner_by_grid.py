"""Corner plot of the full 4-parameter posterior for every Fourier grid, one figure per grid, in one folder per run.

Posterior samples: results/<run>/grid_strategies_chains.npz from compare_grid_strategies.py (each grid at its own
temperature T1; emulator and sampler as the pipeline). Every corner is drawn on the full prior box (|Omega_b| 15-45,
phi_b 20-30 deg, S and bar_frac 0.5-1.5), so all grids share the same axes and a posterior cut by a box edge shows
as one piled against the frame. Titles give the median and 16-84%. The empty upper-right corner carries the grid's
name, its design in words, N_data, T1 and best chi2/N, and a drawing of the grid in (R, |z|) to scale (blue: fitted
cells, R < 6 kpc; grey: rings at 6-8 kpc that exist but are not fitted; stars above the top layer are not used).
For a mock, the truth is marked in orange.

Usage: python plot_corner_by_grid.py --run osarg_zhang2024_sormani22 --label "OSARG data, Sormani+22 family"
       python plot_corner_by_grid.py --run mw_hunter24_mock_from_osarg_ics --truth 32.5 25 1 1 \
           --label "steady-state mock, Hunter+24 family"
Figures: figures/corner_plots_by_grid/<run>/<NN>_<grid>.png
"""
import argparse
import re

import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP
import test1_nbody_qmc_focused_fourier_stars as fs
from plot_posteriors_by_grid import GRIDS, draw_design

NAMES = [r"$|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]", r"$\phi_b$ [deg]", r"$S$ (bar length)", r"bar_frac (bar amplitude)"]
BOX = [(15, 45), (20, 30), (0.5, 1.5), (0.5, 1.5)]
SLUG = {"default": "01_current_default_0.3kpc_rings_3_layers", "_r0.5_z2": "02_0.5kpc_rings_2_layers",
        "_r0.5_z2_dz2": "03_0.5kpc_rings_one_layer_z2", "_r1_z2_dz2": "04_1kpc_rings_one_layer_z2",
        "_r0.5_z3_dz3": "05_0.5kpc_rings_one_layer_z3", "_r0.3_z2_dz2": "06_0.3kpc_rings_one_layer_z2",
        "_geo20_z2_dz2": "07_geometric_20_rings_one_layer_z2", "_geo32_z2_dz2": "08_geometric_32_rings_one_layer_z2",
        "_ilr01_z2_dz2": "09_0.1kpc_rings_at_ILR_one_layer_z2",
        "_geo20_z2": "10_geometric_20_rings_two_layers_z0to1_z1to2"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--truth", type=float, nargs=4, default=None, metavar=("OMEGA", "PHI", "S", "F"))
    args = ap.parse_args()
    set_style()
    chains = np.load(RESULTS_DIR / args.run / "grid_strategies_chains.npz")
    summ = pd.read_csv(RESULTS_DIR / args.run / "grid_strategies.csv").set_index("tag")
    for lab, gargs, title, note in GRIDS:
        if lab not in chains.files:
            continue
        post = chains[lab].astype(float)
        s = summ.loc[lab]
        fig = plt.figure(figsize=(10, 10))
        corner.corner(post, fig=fig, labels=NAMES, range=BOX, truths=args.truth, truth_color=ACCENT,
                      color=BLUE_RAMP[5], quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".2f", bins=30,
                      titles=[r"$|\Omega_b|$", r"$\phi_b$", r"$S$", "bar_frac"],
                      smooth=1.0, title_kwargs=dict(fontsize=10), label_kwargs=dict(fontsize=11),
                      plot_datapoints=False, levels=(0.393, 0.865), fill_contours=True)
        for ax in fig.axes:
            ax.tick_params(top=False, right=False, labelsize=8)
        q = np.percentile(post[:, 0], [16, 50, 84])
        width = f"|Ω_b| width (16–84%) {100 * (q[2] - q[0]) / q[1]:.0f}%"
        if args.truth:
            width += f", bias {100 * (q[1] / args.truth[0] - 1):+.1f}%"
        fig.text(0.56, 0.975, args.label, fontsize=11, color=INK_2, ha="left", va="top")
        fig.text(0.56, 0.95, title, fontsize=13, color=INK, ha="left", va="top", fontweight="bold")
        design = re.sub(r"\n(?=\d \|z\|)", "; ", note).replace("\n", " ")
        fig.text(0.56, 0.915, "Grid design: " + design, fontsize=9.5, color=INK,
                 ha="left", va="top", wrap=True)
        fig.text(0.56, 0.80, f"N_data {int(s.N_data)}    T1 {s.T1:.1f}    best χ²/N {s.best_chi2_per_N:.3f}\n{width}"
                 + ("\norange: truth" if args.truth else ""), fontsize=10, color=INK, ha="left", va="top",
                 linespacing=1.4)
        fs.set_grid(*gargs)
        ax = fig.add_axes([0.60, 0.585, 0.33, 0.12])
        draw_design(ax, fs.RING_EDGES, fs.LAYER_EDGES)
        ax.tick_params(top=False, right=False)
        ax.set_title("the grid in (R, |z|), to scale: blue fitted, grey not fitted", fontsize=8.5, color=INK_2,
                     loc="left", pad=3)
        fs.set_grid()
        savefig(fig, f"corner_plots_by_grid/{args.run}/{SLUG[lab]}")
        plt.close(fig)


if __name__ == "__main__":
    main()
