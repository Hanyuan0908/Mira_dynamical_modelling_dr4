"""Face-on maps of the mean Galactocentric radial and azimuthal velocity of Au18 stars in
the MW-like frame (bar at the MW angle). v_phi is defined positive in the sense of rotation."""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

from common import set_style, savefig, CMAP_SEQ_W, INK_2, R0_KPC, BAR_ANGLE_DEG
from prepare_au18_mwframe import OUT_FILE
from plot_au18_density import draw_bar
from plot_dr3_mira_spatial_period_bins import mark_sun

XY_LIM = 10.0      # kpc
PIX = 0.2          # kpc
Z_SLAB = 1.0       # kpc, |z| cut
MIN_PART = 20      # minimum particles per pixel
DR_AZ = 0.25       # kpc, annulus width for the azimuthal-mean v_phi

# diverging blue <-> red with a neutral grey midpoint (negative = inflow, blue)
CMAP_DIV = LinearSegmentedColormap.from_list(
    "blue_red", ["#0d366b", "#256abf", "#6da7ec", "#b7d3f6", "#f0efec",
                 "#f6bfba", "#ec7f78", "#d23a39", "#7e1a1a"])
CMAP_DIV.set_bad("white")


def binned_mean(x, y, w, q, edges, min_count=MIN_PART):
    W, _, _ = np.histogram2d(x, y, bins=[edges, edges], weights=w)
    Q, _, _ = np.histogram2d(x, y, bins=[edges, edges], weights=w * q)
    N, _, _ = np.histogram2d(x, y, bins=[edges, edges])
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = Q / W
    mean[N < min_count] = np.nan
    return mean.T, W.T


def main():
    set_style()
    with h5py.File(OUT_FILE, "r") as f:
        x, y, z, vx, vy, m = (f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "mass"))
        R_bar = float(f.attrs["bar_half_length_kpc"])

    s = np.abs(z) < Z_SLAB
    x, y, vx, vy, m = x[s], y[s], vx[s], vy[s], m[s]
    R = np.hypot(x, y)
    vR = (x * vx + y * vy) / R
    vphi = -(x * vy - y * vx) / R            # disc spins with L_z < 0; make rotation positive

    # azimuthal-mean v_phi(R), to isolate the non-axisymmetric part
    R_edges = np.arange(0, np.sqrt(2) * XY_LIM + DR_AZ, DR_AZ)
    idx = np.digitize(R, R_edges) - 1
    vphi_az = np.bincount(idx, weights=m * vphi) / np.bincount(idx, weights=m)
    dvphi = vphi - vphi_az[idx]

    edges = np.arange(-XY_LIM, XY_LIM + 1e-9, PIX)
    maps = {}
    for name, q in (("vR", vR), ("vphi", vphi), ("dvphi", dvphi)):
        maps[name], Sigma = binned_mean(x, y, m, q, edges)

    inner = np.hypot(*np.meshgrid(0.5 * (edges[1:] + edges[:-1]),
                                  0.5 * (edges[1:] + edges[:-1]))) < 8.0
    lim_R = np.ceil(np.nanpercentile(np.abs(maps["vR"][inner]), 98) / 5) * 5
    lim_dphi = np.ceil(np.nanpercentile(np.abs(maps["dvphi"][inner]), 98) / 5) * 5
    lim_phi = np.ceil(np.nanpercentile(maps["vphi"], 99.5) / 10) * 10

    panels = [
        ("vR", r"$\langle v_R \rangle$ [km s$^{-1}$]", CMAP_DIV, -lim_R, lim_R),
        ("vphi", r"$\langle v_\phi \rangle$ [km s$^{-1}$]", CMAP_SEQ_W, 0, lim_phi),
        ("dvphi", r"$\langle v_\phi \rangle - \overline{v_\phi}(R)$ [km s$^{-1}$]", CMAP_DIV,
         -lim_dphi, lim_dphi),
    ]
    titles = ["Radial velocity", "Azimuthal velocity", "Azimuthal velocity minus azimuthal mean"]

    sm = gaussian_filter(Sigma / PIX**2, sigma=1.0)
    with np.errstate(divide="ignore"):
        lsm = np.log10(sm)
    levels = np.arange(np.floor(np.nanmax(lsm)) - 1.5, np.nanmax(lsm), 0.25)
    centres = 0.5 * (edges[1:] + edges[:-1])

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), gridspec_kw=dict(wspace=0.38))
    for ax, (key, label, cmap, vmin, vmax), title in zip(axes, panels, titles):
        im = ax.pcolormesh(edges, edges, maps[key], cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
        ax.contour(centres, centres, lsm, levels=levels, colors=INK_2, linewidths=0.45, alpha=0.6)
        draw_bar(ax, R_bar, edge_on=False)
        mark_sun(ax, -R0_KPC, 0.0)
        ax.plot(0, 0, "+", color="k", ms=7, mew=1.0, zorder=6)
        ax.set_xlim(-XY_LIM, XY_LIM)
        ax.set_ylim(-XY_LIM, XY_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$ [kpc]")
        ax.set_title(title, color=INK_2)
        cax = ax.inset_axes([1.03, 0.0, 0.045, 1.0])
        cb = fig.colorbar(im, cax=cax, extend="both" if vmin < 0 else "max")
        cb.set_label(label, fontsize=9)
        cb.outline.set_linewidth(0.5)
    axes[0].set_ylabel(r"$y$ [kpc]")

    fig.suptitle(f"Au18 face-on stellar kinematics, |z| < {Z_SLAB:.0f} kpc, mass-weighted, "
                 f"{PIX * 1e3:.0f} pc pixels (N ≥ {MIN_PART}); bar at {BAR_ANGLE_DEG:.0f}° (dashed), "
                 r"disc rotates clockwise; $v_\phi > 0$ along rotation",
                 y=1.0, fontsize=10.5)
    savefig(fig, "au18_velocity_maps_faceon")


if __name__ == "__main__":
    main()
