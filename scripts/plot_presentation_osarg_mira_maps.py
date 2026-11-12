"""Presentation maps (black background): face-on number density and mean radial velocity of the OSARGs, and the
number density of the DR3 Miras.

Samples:
  OSARGs  the Zhang et al. (2024) final sample, 33,704 stars (M_JKs_sigma < 0.4, Gaia DR3 RV), positions and
          velocities exactly as make_osarg_data.final_sample (PL distances, MW-like frame, Sun at x = -8.122 kpc)
  Miras   data/gaia_dr3_lpv_mira.csv with DM_err < 0.4 (54,573 stars); its x, y, z are in the same frame
Both are cut to |z| < 3 kpc, the vertical range the likelihood uses. Maps: 0.2 kpc pixels over |x|, |y| < 10 kpc,
counts (and v_R sums) smoothed with a Gaussian of SMOOTH kpc before dividing, as in the animations. Density is
stars per kpc^2 on a log scale; <v_R> (Galactocentric cylindrical, positive outward) is shown where the smoothed
count corresponds to at least N_MIN stars inside the kernel. The dashed line is the bar axis at the conventional
27 deg from the Sun-centre line, drawn as a guide (it is not fitted here).

Colours: density on the reference blue ramp, darkest step #184f95 (2.6:1 on black, validated --ordinal --mode
dark); <v_R> on a blue-grey-red diverging map with a dark grey midpoint, so zero velocity is distinct from no data.

Usage: python plot_presentation_osarg_mira_maps.py
Figures: figures/presentation/{osarg_zhang2024_face_on_number_density_and_mean_vR_black,
         dr3_mira_face_on_number_density_black}.{png,pdf}
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy.ndimage import gaussian_filter

from common import ROOT, R0_KPC, BAR_ANGLE_DEG
from make_osarg_data import final_sample
from plot_presentation_figures import dark_style, save, BG, TEXT, TEXT_2, AXIS

EXTENT, PIX, SMOOTH, N_MIN, Z_MAX, VLIM = 10.0, 0.2, 0.3, 15, 3.0, 40.0
EDGES = np.arange(-EXTENT, EXTENT + 1e-9, PIX)
DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "blue_ramp_dark", ["#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb", "#f0f6fe"])
DEX = 2.5                  # each density map spans this many decades below its own 99.9th percentile
DENSITY_CMAP.set_bad(BG)
VR_CMAP = LinearSegmentedColormap.from_list(
    "blue_grey_red_dark", ["#9ec5f4", "#3987e5", "#3a3a37", "#e66767", "#f4b4b4"])
VR_CMAP.set_bad(BG)


def smoothed(x, y, w=None):
    h = np.histogram2d(x, y, bins=[EDGES, EDGES], weights=w)[0]
    return gaussian_filter(h, SMOOTH / PIX, mode="constant").T


def density(x, y):
    n = smoothed(x, y) / PIX**2                       # stars per kpc^2
    return np.where(n > 0.5, n, np.nan)


def mean_vR(x, y, vR):
    n, s = smoothed(x, y), smoothed(x, y, vR)
    n_eff = n * 2 * np.pi * (SMOOTH / PIX) ** 2         # stars inside the kernel's effective area
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n_eff >= N_MIN, s / n, np.nan)


def decorate(ax, title, subtitle):
    ax.plot(-R0_KPC, 0, "o", ms=11, mfc="none", mec=TEXT, mew=1.6, zorder=5)
    ax.plot(-R0_KPC, 0, "o", ms=2.5, color=TEXT, zorder=5)
    ax.text(-R0_KPC, -1.0, "Sun", color=TEXT, fontsize=12, ha="center", va="top")
    ax.plot(0, 0, "+", ms=10, mew=1.6, color=TEXT, zorder=5)
    a = np.radians(BAR_ANGLE_DEG)
    ax.plot([-5 * np.cos(a), 5 * np.cos(a)], [5 * np.sin(a), -5 * np.sin(a)], ls=(0, (4, 3)), lw=1.2,
            color=TEXT, alpha=0.7, zorder=4)
    ax.text(5.3 * np.cos(a), -5.3 * np.sin(a) - 0.3, f"bar, {BAR_ANGLE_DEG:g}°", color=TEXT_2, fontsize=11,
            ha="left", va="top")
    ax.set(xlim=(-EXTENT, EXTENT), ylim=(-EXTENT, EXTENT), xlabel=r"$x$  [kpc]", ylabel=r"$y$  [kpc]")
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_title(f"{title}\n", color=TEXT, fontsize=16)
    ax.text(0.5, 1.015, subtitle, transform=ax.transAxes, color=TEXT_2, fontsize=12, ha="center", va="bottom")
    ax.tick_params(length=4, width=1.0)


def colorbar(fig, ax, im, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label, color=TEXT)
    cb.outline.set_edgecolor(AXIS)
    cb.ax.tick_params(colors=TEXT_2, which="both")


def main():
    dark_style()
    o = final_sample()
    keep = np.abs(o["z"]) < Z_MAX
    x, y, vx, vy = (o[k][keep] for k in ("x", "y", "vx", "vy"))
    vR = (x * vx + y * vy) / np.hypot(x, y)
    print(f"OSARGs: {len(o['x']):,} in the final sample, {keep.sum():,} at |z| < {Z_MAX:g} kpc")

    m = pd.read_csv(ROOT / "data" / "gaia_dr3_lpv_mira.csv", usecols=["x", "y", "z", "DM_err"])
    m = m[(m.DM_err < 0.4) & (np.abs(m.z) < Z_MAX)]
    print(f"Miras: {len(m):,} with DM_err < 0.4 and |z| < {Z_MAX:g} kpc")

    d_osarg, d_mira = density(x, y), density(m.x.values, m.y.values)
    # The two samples differ ~2.5x in size, so each map gets its own scale (a shared one washed both out).
    norms = {k: LogNorm(vmin=10 ** (np.log10(t) - DEX), vmax=t) for k, t in
             (("osarg", float(np.nanpercentile(d_osarg, 99.9))), ("mira", float(np.nanpercentile(d_mira, 99.9))))}
    ext = [-EXTENT, EXTENT, -EXTENT, EXTENT]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16.5, 7.0), gridspec_kw=dict(wspace=0.5))
    im = a1.imshow(d_osarg, origin="lower", extent=ext, cmap=DENSITY_CMAP, norm=norms["osarg"],
                   interpolation="nearest")
    colorbar(fig, a1, im, r"stars per kpc$^2$")
    decorate(a1, "OSARGs: number density", f"{keep.sum():,} stars with 6D, |z| < {Z_MAX:g} kpc")
    v = mean_vR(x, y, vR)
    im = a2.imshow(v, origin="lower", extent=ext, cmap=VR_CMAP, vmin=-VLIM, vmax=VLIM, interpolation="nearest")
    colorbar(fig, a2, im, r"$\langle v_R \rangle$  [km s$^{-1}$]")
    decorate(a2, r"OSARGs: mean radial velocity $\langle v_R \rangle$",
             f"Gaussian-smoothed, σ = {SMOOTH:g} kpc; shown where ≥ {N_MIN} stars")
    save(fig, "osarg_zhang2024_face_on_number_density_and_mean_vR_black")

    fig, ax = plt.subplots(figsize=(7.8, 7.0))
    im = ax.imshow(d_mira, origin="lower", extent=ext, cmap=DENSITY_CMAP, norm=norms["mira"],
                   interpolation="nearest")
    colorbar(fig, ax, im, r"stars per kpc$^2$")
    decorate(ax, "Gaia DR3 Miras: number density", f"{len(m):,} stars, DM error < 0.4 mag, |z| < {Z_MAX:g} kpc")
    save(fig, "dr3_mira_face_on_number_density_black")

    r = np.hypot(x, y)
    print(f"peak density: OSARGs {np.nanmax(d_osarg):.0f}, Miras {np.nanmax(d_mira):.0f} stars/kpc^2; colour "
          f"scales {norms['osarg'].vmin:.1f}-{norms['osarg'].vmax:.0f} and {norms['mira'].vmin:.1f}-"
          f"{norms['mira'].vmax:.0f}")
    print(f"<v_R> map: {np.isfinite(v).sum()} pixels shown, range {np.nanmin(v):.0f} to {np.nanmax(v):.0f} km/s; "
          f"OSARGs at R < 5 kpc: {(r < 5).sum():,}")


if __name__ == "__main__":
    main()
