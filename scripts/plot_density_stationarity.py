"""Does the sample's density stay put? Face-on relative density change of the orbit libraries (time-averaged,
1/tau-weighted mass in the saved 0.5 kpc boxes, |z| < 1 kpc) in trial potentials vs the ground-truth library,
and the data vs the ground-truth library. Bar frame, folded by point symmetry, boxes wholly inside R < 8 kpc,
each map normalised to its total over those boxes. No new integrations.
Figure: figures/test1_nbody_density_stationarity.png
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK_2, ACCENT
from plot_au18_velocity_maps import CMAP_DIV
import steady_state_likelihood as ssl
from test1_nbody_qmc import ARR_DIR, MOCK, TRUTH
from plot_SF_square_maps import grey_bad

Z_CUT = 1.0
E = ssl.BOX_XY
XC = 0.5 * (E[1:] + E[:-1])
XX, YY = np.meshgrid(XC, XC, indexing="ij")
INSIDE = np.hypot(XX, YY) + 0.25 * np.sqrt(2) < ssl.R_SEL
KZ = np.abs(0.5 * (ssl.BOX_Z[1:] + ssl.BOX_Z[:-1])) < Z_CUT
MODELS = [("S 1.25 (longer)", "maps_S1.25_f1.npz"), ("S 0.8 (shorter)", "maps_S0.8_f1.npz"),
          ("f 1.4 (stronger)", "maps_S1_f1.4.npz"), ("f 0.7 (weaker)", "maps_S1_f0.7.npz"),
          ("Sobol #15: S 0.82, f 1.50 (ridge)", "015.npz")]


def lib_density(fn):
    with np.load(ARR_DIR / fn) as a:
        d = a["box_S0"][:, :, KZ].sum(2)
    d = np.where(INSIDE, d, np.nan)
    return d / np.nansum(d)


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        w = ssl.rotate_z(np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")]),
                         np.radians(TRUTH["phi_b_deg"]))
    w = w[np.abs(w[:, 2]) < Z_CUT]
    N = np.histogram2d(np.r_[w[:, 0], -w[:, 0]], np.r_[w[:, 1], -w[:, 1]], bins=[E, E])[0]   # folded counts
    N = np.where(INSIDE, N, np.nan)
    truth = lib_density("000.npz")
    rel = {lab: lib_density(fn) / truth - 1 for lab, fn in MODELS}
    rel["data (Poisson error ±1/√N, median %.0f%%)" % (100 * np.nanmedian(1 / np.sqrt(N)))] = (N / np.nansum(N)) / truth - 1
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 9.6), gridspec_kw=dict(hspace=0.3, wspace=0.35))
    for ax, (lab, r) in zip(axes.flat, rel.items()):
        inner = INSIDE & (np.hypot(XX, YY) < 6)
        rms = np.sqrt(np.nanmean(r[inner]**2))
        im = ax.pcolormesh(E, E, (100 * r).T, cmap=grey_bad(CMAP_DIV), vmin=-20, vmax=20, rasterized=True)
        ax.plot([-5, 5], [0, 0], color=ACCENT, lw=0.9, ls=(0, (4, 2)))
        ax.set_aspect("equal")
        ax.set_xlim(-8, 8)
        ax.set_ylim(-8, 8)
        ax.tick_params(labelsize=7)
        ax.set_title(f"{lab}\nvs ground-truth library; rms (R < 6) {100 * rms:.1f}%", color=INK_2, fontsize=9)
        cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
        cb.set_label("Δ density [%]", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        ax.set_xlabel("x [kpc] (bar frame)", fontsize=8)
        ax.set_ylabel("y [kpc]", fontsize=8)
        print(f"{lab:45s} rms relative density change (R < 6 kpc, |z| < 1): {100 * rms:.1f}%")
    fig.suptitle("Density stationarity: time-averaged library density in trial potentials relative to the "
                 "ground-truth library\n0.5 kpc squares, |z| < 1 kpc, bar frame (bar along x, dashed); "
                 f"Ω_b = {TRUTH['Omega_abs']}, φ_b = {TRUTH['phi_b_deg']:g}° (Sobol #15: Ω_b 26.77, φ_b 27.6°)",
                 y=0.985, fontsize=11)
    savefig(fig, "test1_nbody_density_stationarity")


if __name__ == "__main__":
    main()
