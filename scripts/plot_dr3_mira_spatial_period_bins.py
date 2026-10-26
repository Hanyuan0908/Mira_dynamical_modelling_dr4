"""Face-on (x-y) and edge-on (x-z) surface density of DR3 Miras in five equal-number
period bins, with the Sun and the adopted bar orientation marked."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm

from common import (load_dr3_mira, set_style, savefig, CMAP_SEQ, ACCENT, INK, INK_2,
                    INK_MUTED, R0_KPC, Z_SUN_KPC, BAR_ANGLE_DEG, BAR_HALF_LENGTH_KPC)

N_PERIOD_BINS = 5
XY_LIM = 15.0     # kpc
Z_LIM = 4.0       # kpc
DXY = 0.25        # kpc pixel in x and y
DZ = 0.1          # kpc pixel in z
Z_STRETCH = 2.0   # vertical exaggeration of the x-z panels

OUTLINE = [pe.withStroke(linewidth=3.5, foreground="white")]


def period_edges(period, n):
    """Equal-number bin edges, rounded to 5 d (outer edges widened to include all stars)."""
    q = np.percentile(period, np.linspace(0, 100, n + 1))
    e = 5.0 * np.round(q / 5.0)
    e[0], e[-1] = np.floor(q[0] / 5.0) * 5.0, np.ceil(q[-1] / 5.0) * 5.0
    return e


def surface_density(a, b, a_edges, b_edges):
    H, _, _ = np.histogram2d(a, b, bins=[a_edges, b_edges])
    return H / (np.diff(a_edges)[0] * np.diff(b_edges)[0])   # stars kpc^-2


def mark_sun(ax, x, y, label=True):
    ax.plot(x, y, "o", ms=9, mfc="white", mec=INK, mew=1.2, zorder=5)
    ax.plot(x, y, "o", ms=2.2, color=INK, zorder=6)
    if label:
        ax.annotate("Sun", (x, y), xytext=(0, -11), textcoords="offset points",
                    ha="center", va="top", fontsize=8, color=INK, path_effects=OUTLINE)


def draw_bar_xy(ax):
    phi = np.radians(BAR_ANGLE_DEG)
    a = BAR_HALF_LENGTH_KPC
    xn, yn = -a * np.cos(phi), a * np.sin(phi)          # near end, l > 0
    ax.plot([-R0_KPC, 0], [0, 0], color=INK_MUTED, lw=0.8, ls=":", zorder=4)
    ax.plot([xn, -xn], [yn, -yn], color=ACCENT, lw=1.4, ls=(0, (4, 2)), zorder=5)
    ax.plot([xn, -xn], [yn, -yn], "o", color=ACCENT, ms=3.5, zorder=5)   # bar ends
    ax.plot(0, 0, "+", color=INK, ms=6, mew=1.0, zorder=6)
    ax.annotate(f"bar, {BAR_ANGLE_DEG:.0f}°", (xn, yn), xytext=(-4, 4),
                textcoords="offset points", ha="right", va="bottom", fontsize=8,
                color=INK, path_effects=OUTLINE)


def draw_bar_xz(ax):
    xe = BAR_HALF_LENGTH_KPC * np.cos(np.radians(BAR_ANGLE_DEG))
    ax.plot([-xe, xe], [0, 0], color=ACCENT, lw=1.4, ls=(0, (4, 2)), zorder=5)
    ax.plot([-xe, xe], [0, 0], "o", color=ACCENT, ms=3.5, zorder=5)
    ax.plot(0, 0, "+", color=INK, ms=6, mew=1.0, zorder=6)


def main():
    set_style()
    df = load_dr3_mira()
    p_edges = period_edges(df["period"].values, N_PERIOD_BINS)

    xy_edges = np.arange(-XY_LIM, XY_LIM + 1e-9, DXY)
    z_edges = np.arange(-Z_LIM, Z_LIM + 1e-9, DZ)

    maps_xy, maps_xz, info = [], [], []
    for lo, hi in zip(p_edges[:-1], p_edges[1:]):
        s = df[(df["period"] >= lo) & (df["period"] < hi)]
        maps_xy.append(surface_density(s["x"], s["y"], xy_edges, xy_edges))
        maps_xz.append(surface_density(s["x"], s["z"], xy_edges, z_edges))
        info.append((lo, hi, len(s), np.nanmedian(s["age_both"])))

    # one colour scale per projection, shared across period bins; floor = 1 star per pixel
    norm_xy = LogNorm(vmin=1 / DXY**2, vmax=max(m.max() for m in maps_xy))
    norm_xz = LogNorm(vmin=1 / (DXY * DZ), vmax=max(m.max() for m in maps_xz))

    w = 3.0
    h_xz = w * (2 * Z_LIM * Z_STRETCH) / (2 * XY_LIM)
    fig, axes = plt.subplots(2, N_PERIOD_BINS, figsize=(N_PERIOD_BINS * w + 1.2, w + h_xz + 1.0),
                             sharex=True, gridspec_kw=dict(height_ratios=[w, h_xz],
                                                           hspace=0.06, wspace=0.06))

    for j, (Hxy, Hxz, (lo, hi, n, age)) in enumerate(zip(maps_xy, maps_xz, info)):
        ax, axz = axes[0, j], axes[1, j]
        mxy = ax.pcolormesh(xy_edges, xy_edges, np.ma.masked_equal(Hxy, 0).T, cmap=CMAP_SEQ,
                            norm=norm_xy, rasterized=True)
        mxz = axz.pcolormesh(xy_edges, z_edges, np.ma.masked_equal(Hxz, 0).T, cmap=CMAP_SEQ,
                             norm=norm_xz, rasterized=True)

        draw_bar_xy(ax)
        mark_sun(ax, -R0_KPC, 0.0)
        draw_bar_xz(axz)
        mark_sun(axz, -R0_KPC, Z_SUN_KPC, label=False)

        ax.set_title(f"{lo:.0f}–{hi:.0f} d\nN = {n:,}, median age {age:.1f} Gyr")
        ax.set_xticks([-10, -5, 0, 5, 10])
        ax.set_yticks([-10, -5, 0, 5, 10])
        ax.set_xlim(-XY_LIM, XY_LIM)
        ax.set_ylim(-XY_LIM, XY_LIM)
        ax.set_aspect("equal")
        axz.set_ylim(-Z_LIM, Z_LIM)
        axz.set_aspect(Z_STRETCH)
        axz.set_xlabel(r"$x$ [kpc]")
        if j == 0:
            ax.set_ylabel(r"$y$ [kpc]")
            axz.set_ylabel(r"$z$ [kpc]" + "\n" + rf"($\times${Z_STRETCH:.0f} stretch)")
        else:
            ax.tick_params(labelleft=False)
            axz.tick_params(labelleft=False)

    for row, m, label in ((0, mxy, r"$\Sigma$ [stars kpc$^{-2}$]"),
                          (1, mxz, r"$\Sigma$ [stars kpc$^{-2}$]")):
        cax = axes[row, -1].inset_axes([1.05, 0.0, 0.05, 1.0])
        cb = fig.colorbar(m, cax=cax)
        cb.set_label(label, fontsize=9)
        cb.outline.set_linewidth(0.5)

    fig.suptitle("Gaia DR3 Miras in equal-number period bins (Galactocentric, Sun at "
                 rf"$x=-{R0_KPC}$ kpc; dashed: bar major axis, half-length {BAR_HALF_LENGTH_KPC:.0f} kpc)",
                 y=1.0, fontsize=11, color=INK)
    savefig(fig, "dr3_mira_spatial_period_bins")


if __name__ == "__main__":
    main()
