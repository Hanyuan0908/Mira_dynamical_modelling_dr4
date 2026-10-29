"""Spatial distribution of a mock catalogue against the DR3 Miras: face-on (x-y), edge-on (x-z) and sky (l, b).

Miras: data/gaia_dr3_lpv_mira.csv with DM_err < 0.4 and Galactocentric R < 8 kpc (the target of make_mock_selection.py);
catalogue x, y, z (astropy Galactocentric frame, Sun at x = -R0). Mock: observed positions and (l, b) of the given file.
Each panel shows log10 of the fraction of the sample per square pixel, with the same pixels and colour range for both
samples in a column.
Usage: python plot_mock_vs_dr3_mira_spatial.py [mock file name in mock/]
Figure: figures/test2_spatial_mock_vs_dr3_miras_<mock file stem>.png
"""
import sys

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, MIRA_DR3_FILE, R0_KPC, Z_SUN_KPC, set_style, savefig, ACCENT, CMAP_SEQ_W

SIGMA_DM_CUT = 0.4
R_MAX = 8.0
PIX = 0.2          # kpc
SKY_PIX = 1.0      # deg


def panels():
    xy = np.arange(-8.6, 8.6 + 1e-9, PIX)
    z = np.arange(-4.0, 4.0 + 1e-9, PIX)
    l = np.arange(-90.0, 90.0 + 1e-9, SKY_PIX)        # all stars with R < 8 kpc lie within |l| < 90 deg
    b = np.arange(-20.0, 20.0 + 1e-9, SKY_PIX)
    return [("x [kpc]", "y [kpc]", xy, xy, "face-on (x–y)"), ("x [kpc]", "z [kpc]", xy, z, "edge-on (x–z)"),
            ("l [deg]", "b [deg]", l, b, "sky (l, b)")]


def main():
    mock_name = sys.argv[1] if len(sys.argv) > 1 else "nbody_mock_selection_sigmaDM0.163_N50k.hdf5"
    set_style()
    df = pd.read_csv(MIRA_DR3_FILE, usecols=["l", "b", "DM_err", "x", "y", "z"]).dropna()
    df = df[(df.DM_err < SIGMA_DM_CUT) & (np.hypot(df.x, df.y) < R_MAX)]
    with h5py.File(MOCK_DIR / mock_name, "r") as f:
        m = {k: f[k][:] for k in ("x", "y", "z", "l", "b")}
        sigma = f.attrs.get("sigma_DM", np.nan)
    wrap = lambda l: ((np.asarray(l) + 180.0) % 360.0) - 180.0
    samples = [(f"DR3 Miras (σ_DM < {SIGMA_DM_CUT}, R < {R_MAX:g} kpc), {len(df):,}",
                (df.x.values, df.y.values, df.z.values, wrap(df.l.values), df.b.values)),
               (f"mock (observed), {len(m['x']):,}", (m["x"], m["y"], m["z"], wrap(m["l"]), m["b"]))]
    fig, axes = plt.subplots(3, 2, figsize=(15, 15.5), gridspec_kw=dict(hspace=0.28, wspace=0.12,
                                                                        height_ratios=[1.6, 0.85, 0.85]))
    axes = axes.T                                       # axes[sample, projection]
    for c, (xl, yl, ex, ey, title) in enumerate(panels()):
        hs = []
        for _, (x, y, z, l, b) in samples:
            u, v = ((x, y), (x, z), (l, b))[c]
            hs.append(np.histogram2d(u, v, bins=[ex, ey])[0].T / len(u))
        pos = np.concatenate([h[h > 0] for h in hs])
        lo, hi = np.log10(np.percentile(pos, 1)), np.log10(pos.max())
        for r, ((label, _), h) in enumerate(zip(samples, hs)):
            ax = axes[r, c]
            im = ax.pcolormesh(ex, ey, np.log10(np.where(h > 0, h, np.nan)), cmap=CMAP_SEQ_W, vmin=lo, vmax=hi,
                               rasterized=True)
            if c < 2:
                ax.plot(-R0_KPC, Z_SUN_KPC if c == 1 else 0.0, marker="*", color=ACCENT, ms=11, ls="none")
                ax.set_aspect("equal")
            else:
                ax.set_xlim(ex[-1], ex[0])        # l increasing to the left
                ax.set_aspect("equal")
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            ax.set_title(f"{label}\n{title}", fontsize=9.5)
            cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
            cb.set_label("log10 fraction per pixel", fontsize=8)
            cb.ax.tick_params(labelsize=7)
    fig.suptitle(f"Spatial distribution: DR3 Miras vs {mock_name} (σ_DM = {sigma:.3f} mag); {PIX} kpc and "
                 f"{SKY_PIX:g}° pixels; star = Sun", fontsize=11, y=0.995)
    savefig(fig, f"test2_spatial_mock_vs_dr3_miras_{mock_name.replace('.hdf5', '')}")


if __name__ == "__main__":
    main()
