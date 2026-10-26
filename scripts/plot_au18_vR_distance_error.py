"""Effect of heliocentric distance errors on the face-on <v_R> map of Au18.

Stars are observed from the Sun (l, b, d, pm_l*, pm_b, v_los), the distance is scattered with
a Gaussian fractional error while the angles, proper motions and v_los are kept, and the stars
are converted back to Galactocentric coordinates. Maps use the |z| slab in each frame."""
import resource

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK_2, INK_MUTED, R0_KPC
from prepare_au18_mwframe import OUT_FILE
from plot_au18_density import draw_bar
from plot_au18_velocity_maps import binned_mean, CMAP_DIV, XY_LIM, PIX, Z_SLAB, MIN_PART
from plot_dr3_mira_spatial_period_bins import mark_sun
import observe

FRAC_DIST_ERR = 0.2
SEED = 42
LOS_LONGITUDES = (-30, -15, 0, 15, 30)   # deg, guide lines of sight


def vR_map(x, y, z, vx, vy, m, edges):
    s = np.abs(z) < Z_SLAB
    R = np.hypot(x[s], y[s])
    vR = (x[s] * vx[s] + y[s] * vy[s]) / R
    return binned_mean(x[s], y[s], m[s], vR, edges)


def quadrupole_node(x, y, z, vx, vy, m, R_range=(1.0, 4.0)):
    """Angle [deg] of the v_R quadrupole node that lies 45 deg behind the maximum of the m=2 term
    (for a bar this node is the bar major axis); also returns the m=2 amplitude [km/s]."""
    s = (np.abs(z) < Z_SLAB) & (np.hypot(x, y) > R_range[0]) & (np.hypot(x, y) < R_range[1])
    phi = np.arctan2(y[s], x[s])
    vR = (x[s] * vx[s] + y[s] * vy[s]) / np.hypot(x[s], y[s])
    w = m[s] / m[s].sum()
    a, b = 2 * np.sum(w * vR * np.cos(2 * phi)), 2 * np.sum(w * vR * np.sin(2 * phi))
    phi_max = 0.5 * np.degrees(np.arctan2(b, a))
    return (phi_max - 45 + 90) % 180 - 90, np.hypot(a, b)


def draw_los(ax):
    for l in LOS_LONGITUDES:
        t = np.radians(l)
        ax.plot([-R0_KPC, -R0_KPC + 30 * np.cos(t)], [0, 30 * np.sin(t)], color=INK_MUTED,
                lw=0.6, ls=":", zorder=4)


def main():
    set_style()
    with h5py.File(OUT_FILE, "r") as f:
        x, y, z, vx, vy, vz, m = (f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz", "mass"))
        R_bar = float(f.attrs["bar_half_length_kpc"])

    obs = observe.to_galactic(x, y, z, vx, vy, vz)
    rng = np.random.default_rng(SEED)
    d_obs = observe.scatter_distance(obs["d"], FRAC_DIST_ERR, rng)
    xo, yo, zo, vxo, vyo, vzo = observe.to_galactocentric(obs["l"], obs["b"], d_obs, obs["pml"],
                                                         obs["pmb"], obs["vlos"])

    # sanity: the noise-free round trip must return the input
    chk = observe.to_galactocentric(obs["l"], obs["b"], obs["d"], obs["pml"], obs["pmb"], obs["vlos"])
    err = max(np.abs(chk[0] - x).max(), np.abs(chk[4] - vy).max())
    print(f"round-trip max error: {err:.2e}")
    assert err < 1e-6

    for name, args in (("true", (x, y, z, vx, vy, m)), ("observed", (xo, yo, zo, vxo, vyo, m))):
        node, amp = quadrupole_node(*args)
        print(f"{name:9s} v_R m=2 (1 < R < 4 kpc): amplitude {amp:5.1f} km/s, node at {node:6.1f} deg")

    edges = np.arange(-XY_LIM, XY_LIM + 1e-9, PIX)
    vR_true, S_true = vR_map(x, y, z, vx, vy, m, edges)
    vR_obs, S_obs = vR_map(xo, yo, zo, vxo, vyo, m, edges)
    diff = vR_obs - vR_true

    lim = 40.0
    panels = [
        (vR_true, S_true, "True positions and velocities", r"$\langle v_R \rangle$ [km s$^{-1}$]"),
        (vR_obs, S_obs, rf"With $\sigma_d/d$ = {FRAC_DIST_ERR:.1f} distance errors",
         r"$\langle v_R \rangle$ [km s$^{-1}$]"),
        (diff, None, "Observed minus true", r"$\Delta\langle v_R \rangle$ [km s$^{-1}$]"),
    ]
    centres = 0.5 * (edges[1:] + edges[:-1])
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), gridspec_kw=dict(wspace=0.38))
    for ax, (mp, S, title, label) in zip(axes, panels):
        im = ax.pcolormesh(edges, edges, mp, cmap=CMAP_DIV, vmin=-lim, vmax=lim, rasterized=True)
        if S is not None:
            with np.errstate(divide="ignore"):
                lS = np.log10(S / PIX**2)
            levels = np.arange(np.floor(lS.max()) - 1.5, lS.max(), 0.25)
            ax.contour(centres, centres, lS, levels=levels, colors=INK_2, linewidths=0.45, alpha=0.6)
        draw_los(ax)
        draw_bar(ax, R_bar, edge_on=False)
        mark_sun(ax, -R0_KPC, 0.0)
        ax.plot(0, 0, "+", color="k", ms=7, mew=1.0, zorder=6)
        ax.set_xlim(-XY_LIM, XY_LIM)
        ax.set_ylim(-XY_LIM, XY_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$ [kpc]")
        ax.set_title(title, color=INK_2)
        cax = ax.inset_axes([1.03, 0.0, 0.045, 1.0])
        cb = fig.colorbar(im, cax=cax, extend="both")
        cb.set_label(label, fontsize=9)
        cb.outline.set_linewidth(0.5)
    axes[0].set_ylabel(r"$y$ [kpc]")
    fig.suptitle(f"Au18, |z| < {Z_SLAB:.0f} kpc, {PIX * 1e3:.0f} pc pixels (N ≥ {MIN_PART}): Gaussian "
                 f"heliocentric distance errors, $\\sigma_d/d$ = {FRAC_DIST_ERR}, with $l, b, \\mu_{{l*}}, "
                 f"\\mu_b, v_{{\\rm los}}$ kept; dotted lines of sight at $l$ = "
                 f"{', '.join(str(l) for l in LOS_LONGITUDES)}°", y=1.0, fontsize=10)
    savefig(fig, "au18_vR_distance_error")
    print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
