"""Where in the sample does the second (low-|Omega_b|, weak-bar) mode come from? Maps of the per-block Delta chi2.

Reads results/bimodality_origin_<case>.npz/.json (bimodality_spatial_origin.py). One row per Fourier grid:
  left:   Delta chi2(low mode - main mode) per (ring, |z| layer) block, drawn to scale in (R, |z|), summed over the four
          fields of scheme B. Blue: the block's data prefer the LOW mode; red: they prefer the main mode. Grey: block
          not valid in every run. The number in each R band and layer is the summed Delta, +- its standard error
          over the 12 replicate runs.
  middle: the same for (main - truth), mocks only.
  right:  Delta per field, for each pair.
Usage: python plot_bimodality_spatial_origin.py --case h24mock
Figure: figures/bimodality_spatial_origin_<case>.png
"""
import argparse
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED

CMAP = LinearSegmentedColormap.from_list("div", ["#1f5fbf", "#9ec5f4", "#f2f1ee", "#f4a582", "#b2182b"])
TITLES = {"h24mock": "Steady-state mock in Hunter+24 (truth 32.5, 25, 1, 1)", "osarg_h24": "OSARG data in Hunter+24",
          "osarg_s22": "OSARG data in Sormani+22"}
PAIR_LABEL = {"low-main": "low mode − main mode", "second-main": "second mode − main mode",
              "second-peak": "second peak (Ω_b 37.6) − main peak (Ω_b 30.0)",
              "main-peak": "strong-bar ridge (Ω_b 34.5) − main peak (Ω_b 30.0)",
              "main-truth": "main mode − truth", "low-truth": "low mode − truth"}


def draw_map(ax, D, edges, layers, vmax, band_text):
    norm = TwoSlopeNorm(0, -vmax, vmax)
    nl = len(layers) - 1
    for i in range(len(D)):
        r, z = i // nl, i % nl
        if edges[r + 1] > 6.0 + 1e-9:
            continue
        x0, w, y0, h = edges[r], edges[r + 1] - edges[r], layers[z], layers[z + 1] - layers[z]
        col = "#d9d8d4" if not np.isfinite(D[i]) else CMAP(norm(D[i]))
        ax.add_patch(Rectangle((x0, y0), w, h, facecolor=col, edgecolor="white", lw=0.6))
    for (lo, hi, z), (v, s) in band_text.items():
        ax.text(0.5 * (lo + hi), 0.5 * (layers[z] + layers[z + 1]), f"{v:+.0f}\n±{s:.0f}", ha="center", va="center",
                fontsize=7.5, color=INK)
    ax.set(xlim=(0, 6), ylim=(0, layers[-1]), xlabel="R [kpc]", ylabel="|z| [kpc]")
    ax.tick_params(top=False, right=False, labelsize=8)
    return norm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True)
    args = ap.parse_args()
    set_style()
    A = np.load(RESULTS_DIR / f"bimodality_origin_{args.case}.npz")
    J = json.load(open(RESULTS_DIR / f"bimodality_origin_{args.case}.json"))
    grids = list(J["grids"])
    sec = "low" if "low" in J["points"] else "second"
    pairs = ([f"{sec}-peak", "main-peak"] if "peak" in J["points"] else [f"{sec}-main"]) \
        + (["main-truth"] if "truth" in J["points"] else [])
    fields = ["<v_R> m2,4", "<v_phi> m0,2,4", "sigma_R^2 m2,4", "sigma_phi^2 m2,4"]
    vmax = 0
    for g in grids:
        for p in pairs:
            D = np.nansum(A[f"{g}__{p}__delta_block_field"], axis=1)
            vmax = max(vmax, np.nanmax(np.abs(np.where(np.isnan(A[f"{g}__{p}__delta_block_field"][:, 0]), np.nan, D))))
    ncol = len(pairs) + 1
    fig, axes = plt.subplots(len(grids), ncol, figsize=(5.6 * ncol, 3.3 * len(grids) + 1.2), squeeze=False,
                             gridspec_kw=dict(width_ratios=[1.3] * len(pairs) + [0.9]), layout="constrained")
    P = J["points"]
    pt = lambda k: f"Ω_b {P[k]['Omega']:.1f}, φ_b {P[k]['phi']:.1f}, S {P[k]['S']:.2f}, bar_frac {P[k]['f']:.2f}"
    for row, g in enumerate(grids):
        G = J["grids"][g]
        edges, layers = A[f"{g}__ring_edges"], A[f"{g}__layer_edges"]
        for col, p in enumerate(pairs):
            ax = axes[row, col]
            D = A[f"{g}__{p}__delta_block_field"]
            Dsum = np.where(np.isnan(D[:, 0]), np.nan, np.nansum(D, axis=1))
            bt = {}
            for k, (v, s) in G["pairs"][p]["by_band_layer"].items():
                rr, zz = k.split("_z")
                lo, hi = map(float, rr[1:].split("-"))
                bt[(lo, hi, int(zz))] = (v, s)
            norm = draw_map(ax, Dsum, edges, layers, vmax, bt)
            tot = G["pairs"][p]
            T1 = G.get("T1")
            ax.set_title((f"{G['desc']}\n" if col == 0 else "\n") + f"{PAIR_LABEL[p]}: Δχ² = {tot['delta']:+.1f} ± {tot['se']:.1f}"
                         + (f" ({tot['delta'] / T1:+.2f} T1)" if T1 else ""), fontsize=9, color=INK, loc="left")
        ax = axes[row, -1]
        y = np.arange(len(fields))
        for j, p in enumerate(pairs):
            v = [G["pairs"][p]["by_field"][f] for f in fields]
            ax.barh(y + (j - (len(pairs) - 1) / 2) * 0.36, v, height=0.34, color=["#4a7fd4", INK_MUTED, "#e07b28"][j],
                    label=PAIR_LABEL[p])
        ax.axvline(0, color=INK_2, lw=0.8)
        ax.set_yticks(y, ["⟨v_R⟩ m2,4", "⟨v_φ⟩ m0,2,4", "σ_R² m2,4", "σ_φ² m2,4"], fontsize=8)
        ax.tick_params(top=False, right=False, labelsize=8)
        ax.set_xlabel("Δχ² by field", fontsize=9)
        if row == 0:
            ax.legend(frameon=False, fontsize=8, loc="lower right")
    sm = plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax), cmap=CMAP)
    cb = fig.colorbar(sm, ax=axes[:, len(pairs) - 1].ravel().tolist(), shrink=0.5, pad=0.02)
    cb.set_label("Δχ² per block (blue: the data prefer the first point)", fontsize=8)
    fig.suptitle(f"{TITLES.get(args.case, args.case)}: which blocks prefer which mode (12 replicate runs per point)\n"
                 f"{sec} mode: {pt(sec)};  main mode: {pt('main')}"
                 + (f";  posterior main peak: {pt('peak')}" if "peak" in P else "")
                 + (f";  truth: {pt('truth')}" if "truth" in P else ""), fontsize=10, color=INK, x=0.01, ha="left")
    savefig(fig, f"bimodality_spatial_origin_{args.case}")


if __name__ == "__main__":
    main()
