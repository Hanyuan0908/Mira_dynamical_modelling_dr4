"""Overview of a mock catalogue: face-on and edge-on star counts, and face-on mean v_R and
v_phi (all stars, no |z| cut). Usage: python plot_mock_overview.py mock/<file>.hdf5"""
import argparse
from pathlib import Path

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from common import set_style, savefig, CMAP_SEQ, CMAP_SEQ_W, INK_2, R0_KPC, Z_SUN_KPC
from plot_au18_density import draw_bar
from plot_au18_velocity_maps import binned_mean, CMAP_DIV
from plot_dr3_mira_spatial_period_bins import mark_sun
from plot_au18_vR_distance_error import quadrupole_node

LIM = 9.0          # kpc, axis half-width (shows the Sun)
Z_LIM = 4.5        # kpc
Z_STRETCH = 2.0
PIX_N = 0.2        # kpc, pixels of the count maps
PIX_V = 0.5        # kpc, pixels of the velocity maps
MIN_STARS = 10     # minimum stars per velocity pixel
V_R_LIM = 40.0     # km/s
N_BOOT = 200


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mock", type=Path)
    args = p.parse_args()

    set_style()
    with h5py.File(args.mock, "r") as f:
        x, y, z, vx, vy = (f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy"))
        R_bar = float(f.attrs["bar_half_length_kpc"])
        R_max = float(f.attrs["R_max"])
    n = len(x)
    R = np.hypot(x, y)
    vR = (x * vx + y * vy) / R
    vphi = -(x * vy - y * vx) / R           # positive along (clockwise) rotation

    # v_R quadrupole (|z| < 1 kpc, 1 < R < 4 kpc) with bootstrap errors
    node, amp = quadrupole_node(x, y, z, vx, vy, np.ones(n))
    rng = np.random.default_rng(0)
    boot = np.array([quadrupole_node(x[i], y[i], z[i], vx[i], vy[i], np.ones(n))
                     for i in (rng.integers(0, n, n) for _ in range(N_BOOT))])
    print(f"v_R m=2 (|z|<1, 1<R<4 kpc): amplitude {amp:.1f} +/- {boot[:, 1].std():.1f} km/s, "
          f"node {node:.1f} +/- {boot[:, 0].std():.1f} deg")

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 10.4), gridspec_kw=dict(hspace=0.22, wspace=0.38))
    (ax_xy, ax_xz), (ax_vR, ax_vphi) = axes

    e = np.arange(-LIM, LIM + 1e-9, PIX_N)
    ez = np.arange(-Z_LIM, Z_LIM + 1e-9, PIX_N / Z_STRETCH)
    Nxy, _, _ = np.histogram2d(x, y, bins=[e, e])
    Nxz, _, _ = np.histogram2d(x, z, bins=[e, ez])
    norm = LogNorm(vmin=1, vmax=max(Nxy.max(), Nxz.max()))
    im_n = ax_xy.pcolormesh(e, e, np.ma.masked_equal(Nxy, 0).T, cmap=CMAP_SEQ, norm=norm, rasterized=True)
    ax_xz.pcolormesh(e, ez, np.ma.masked_equal(Nxz, 0).T, cmap=CMAP_SEQ, norm=norm, rasterized=True)

    ev = np.arange(-LIM, LIM + 1e-9, PIX_V)
    w = np.ones_like(x)
    m_vR, _ = binned_mean(x, y, w, vR, ev, min_count=MIN_STARS)
    m_vphi, _ = binned_mean(x, y, w, vphi, ev, min_count=MIN_STARS)
    im_vR = ax_vR.pcolormesh(ev, ev, m_vR, cmap=CMAP_DIV, vmin=-V_R_LIM, vmax=V_R_LIM, rasterized=True)
    vphi_max = np.ceil(np.nanpercentile(m_vphi, 99) / 10) * 10
    im_vphi = ax_vphi.pcolormesh(ev, ev, m_vphi, cmap=CMAP_SEQ_W, vmin=0, vmax=vphi_max, rasterized=True)

    t = np.linspace(0, 2 * np.pi, 400)
    for ax in (ax_xy, ax_vR, ax_vphi):
        ax.plot(R_max * np.cos(t), R_max * np.sin(t), color=INK_2, lw=0.6, ls=":")
        draw_bar(ax, R_bar, edge_on=False)
        mark_sun(ax, -R0_KPC, 0.0)
        ax.plot(0, 0, "+", color="k", ms=7, mew=1.0, zorder=6)
        ax.set_ylim(-LIM, LIM)
        ax.set_aspect("equal")
        ax.set_ylabel(r"$y$ [kpc]")
    draw_bar(ax_xz, R_bar, edge_on=True)
    mark_sun(ax_xz, -R0_KPC, Z_SUN_KPC, label=False)
    ax_xz.plot(0, 0, "+", color="k", ms=7, mew=1.0, zorder=6)
    ax_xz.set_ylim(-Z_LIM, Z_LIM)
    ax_xz.set_aspect(Z_STRETCH)
    ax_xz.set_ylabel(rf"$z$ [kpc] ($\times${Z_STRETCH:.0f})")
    for ax in axes.flat:
        ax.set_xlim(-LIM, LIM)
        ax.set_xlabel(r"$x$ [kpc]")

    for ax, im, label, title, ext in (
            (ax_xy, im_n, f"stars per ({PIX_N * 1e3:.0f} pc)$^2$", "Face-on star counts", "neither"),
            (ax_xz, im_n, f"stars per {PIX_N * 1e3:.0f}×{PIX_N / Z_STRETCH * 1e3:.0f} pc pixel",
             "Edge-on star counts", "neither"),
            (ax_vR, im_vR, r"$\langle v_R \rangle$ [km s$^{-1}$]",
             f"Radial velocity ({PIX_V * 1e3:.0f} pc pixels, N ≥ {MIN_STARS})", "both"),
            (ax_vphi, im_vphi, r"$\langle v_\phi \rangle$ [km s$^{-1}$]",
             f"Azimuthal velocity ({PIX_V * 1e3:.0f} pc pixels, N ≥ {MIN_STARS})", "max")):
        ax.set_title(title, color=INK_2)
        cax = ax.inset_axes([1.03, 0.0, 0.045, 1.0])
        cb = fig.colorbar(im, cax=cax, extend=ext)
        cb.set_label(label, fontsize=9)
        cb.outline.set_linewidth(0.5)

    fig.suptitle(f"{args.mock.stem}: {n:,} stars, R < {R_max:g} kpc (dotted); all z; "
                 r"bar dashed, $v_\phi > 0$ along rotation", y=0.95, fontsize=11)
    savefig(fig, f"{args.mock.stem}_overview")


if __name__ == "__main__":
    main()
