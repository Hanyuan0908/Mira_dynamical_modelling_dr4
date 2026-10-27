"""Face-on (x-y) and edge-on (x-z) stellar surface density of Au18 in the MW-like frame,
wide view (same extent and stretch as the DR3 Mira maps) and a zoom on the bar."""
import argparse

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter

from common import (set_style, savefig, CMAP_SEQ_W, INK_2, ACCENT, R0_KPC, Z_SUN_KPC,
                    BAR_ANGLE_DEG)
from galaxies import galaxy, add_argument as add_galaxy_argument
from plot_dr3_mira_spatial_period_bins import mark_sun

VIEWS = {  # name: (x half-width, z half-width, pixel [kpc], z stretch)
    "wide": (15.0, 4.0, 0.1, 2.0),
    "zoom": (6.0, 3.0, 0.05, 1.0),
}
DYNAMIC_RANGE_DEX = 3.0
CONTOUR_STEP_DEX = 0.25


def surface_density(a, b, mass, a_lim, b_lim, pix):
    ea = np.arange(-a_lim, a_lim + 1e-9, pix)
    eb = np.arange(-b_lim, b_lim + 1e-9, pix)
    H, _, _ = np.histogram2d(a, b, bins=[ea, eb], weights=mass)
    return ea, eb, H.T / pix**2      # Msun kpc^-2


def draw_density(ax, ea, eb, H, norm, pix):
    im = ax.pcolormesh(ea, eb, np.ma.masked_equal(H, 0), cmap=CMAP_SEQ_W, norm=norm, rasterized=True)
    smooth = gaussian_filter(H, sigma=0.15 / pix)                  # ~150 pc smoothing
    levels = 10 ** np.arange(np.log10(norm.vmin) + 1.0, np.log10(norm.vmax), CONTOUR_STEP_DEX)
    ac, bc = 0.5 * (ea[1:] + ea[:-1]), 0.5 * (eb[1:] + eb[:-1])
    ax.contour(ac, bc, smooth, levels=levels, colors=INK_2, linewidths=0.45, alpha=0.7)
    return im


def draw_bar(ax, half_length, edge_on):
    phi = np.radians(BAR_ANGLE_DEG)
    xn, yn = -half_length * np.cos(phi), half_length * np.sin(phi)
    xs, ys = [xn, -xn], ([0, 0] if edge_on else [yn, -yn])
    ax.plot(xs, ys, color=ACCENT, lw=1.4, ls=(0, (4, 2)), zorder=5)
    ax.plot(xs, ys, "o", color=ACCENT, ms=3.5, zorder=5)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    add_galaxy_argument(ap)
    g = galaxy(ap.parse_args().galaxy)
    set_style()
    with h5py.File(g["stars"], "r") as f:
        x, y, z, m = (f[k][:] for k in ("x", "y", "z", "mass"))
        R_bar = float(f.attrs["bar_half_length_kpc"])

    w = 4.2
    hz = {k: w * (2 * zl * st) / (2 * xl) for k, (xl, zl, _, st) in VIEWS.items()}
    fig = plt.figure(figsize=(2 * w + 2.2, w + max(hz.values()) + 0.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[w, max(hz.values())], hspace=0.08, wspace=0.42)

    for j, (name, (xl, zl, pix, stretch)) in enumerate(VIEWS.items()):
        ax, axz = fig.add_subplot(gs[0, j]), fig.add_subplot(gs[1, j])
        ex, ey, Hxy = surface_density(x, y, m, xl, xl, pix)
        ex, ez, Hxz = surface_density(x, z, m, xl, zl, pix)
        vmax = max(Hxy.max(), Hxz.max())
        norm = LogNorm(vmin=vmax / 10**DYNAMIC_RANGE_DEX, vmax=vmax)
        im = draw_density(ax, ex, ey, Hxy, norm, pix)
        draw_density(axz, ex, ez, Hxz, norm, pix)

        draw_bar(ax, R_bar, edge_on=False)
        draw_bar(axz, R_bar, edge_on=True)
        if R0_KPC < xl:
            mark_sun(ax, -R0_KPC, 0.0)
            mark_sun(axz, -R0_KPC, Z_SUN_KPC, label=False)
        for a in (ax, axz):
            a.plot(0, 0, "+", color="white", ms=7, mew=1.2, zorder=6)

        ax.set_xlim(-xl, xl)
        ax.set_ylim(-xl, xl)
        ax.set_aspect("equal")
        axz.set_xlim(-xl, xl)
        axz.set_ylim(-zl, zl)
        axz.set_aspect(stretch)
        ax.tick_params(labelbottom=False)
        ax.set_ylabel(r"$y$ [kpc]")
        axz.set_xlabel(r"$x$ [kpc]")
        axz.set_ylabel(r"$z$ [kpc]" + (rf" ($\times${stretch:.0f})" if stretch != 1 else ""))
        ax.set_title(("Wide view" if name == "wide" else "Bar region") +
                     f" ({pix * 1e3:.0f} pc pixels)", color=INK_2)

        cax = ax.inset_axes([1.04, 0.0, 0.045, 1.0])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(r"$\Sigma_\star$ [M$_\odot$ kpc$^{-2}$]", fontsize=9)
        cb.outline.set_linewidth(0.5)

    fig.suptitle(f"{g['label']} stars in a MW-like frame\n"
                 f"bar at {BAR_ANGLE_DEG:.0f}° to the Sun–GC line, half-length {R_bar:.2f} kpc "
                 rf"(dashed); Sun at $x=-{R0_KPC}$ kpc; contours every {CONTOUR_STEP_DEX} dex",
                 y=1.02, fontsize=10.5)
    savefig(fig, f"{g['prefix']}_stellar_density_mwframe")


if __name__ == "__main__":
    main()
