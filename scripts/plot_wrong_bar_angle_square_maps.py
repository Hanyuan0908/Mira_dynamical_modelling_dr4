"""What do the kinematic maps look like with a wrong bar angle?

The potential's bar is always along x. The mock stars are rotated by the trial angle (true 27 deg; wrong 20 deg) and
integrated in the true potential (Omega_b, S, bar_frac at the truth). Square 0.5 kpc pixels, |z| < 1 kpc, folded by point
symmetry (as in plot_bar_length_strength_square_maps.py). Rows:
  data rotated by 27 deg | library, trial 27 deg | data rotated by 20 deg | library, trial 20 deg |
  library 20 - library 27 | data 20 - library 20 | data 27 - library 27
Dashed line: the potential's bar axis (x). Solid line on data rows: where the real bar lies in that frame
(0 deg for 27, -7 deg for 20). Also prints the m = 2 orientation of data and library (polar cells, as in
check_m2_orientation_vs_trial_angle.py). The 20-deg library is integrated once and cached in
~/mira_dr4_runs/test1_nbody/maps_phi20.npz.
Figure: figures/test1_nbody_square_maps_wrong_bar_angle_20deg_vs_27deg.png
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK, INK_2, ACCENT
from plot_au18_velocity_maps import CMAP_DIV
import steady_state_likelihood as ssl
from test1_nbody_qmc import ARR_DIR, MOCK, TRUTH
from test1_nbody_moments import boot_errors, SEED
from plot_bar_length_strength_square_maps import (data_maps, library_maps, grey_bad, QUANTS, VALUE_LIMS, Z_CUT,
                                                  N_MIN)
from check_m2_orientation_vs_trial_angle import fields, rotation_deg, FIELDS

PHI_TRUE, PHI_WRONG = 27.0, 20.0
CACHE = ARR_DIR / "maps_phi20.npz"


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    if not CACHE.exists():
        ssl.ACCURACY = 1e-5
        _, arr = ssl.evaluate(stars_mw, TRUTH["Omega_abs"], PHI_WRONG, TRUTH["S"], TRUTH["bar_frac"], galaxy="nbody",
                              extended=True)
        np.savez(CACHE, **{k: v for k, v in arr.items() if k.startswith(("box_", "lib_"))})
    data, mask = {}, {}
    for phi in (PHI_TRUE, PHI_WRONG):
        w = ssl.rotate_z(stars_mw, np.radians(phi))
        data[phi] = data_maps(w[np.abs(w[:, 2]) < Z_CUT])
        mask[phi] = data[phi]["N"] < N_MIN
    lib = {PHI_TRUE: library_maps(ARR_DIR / "000.npz"), PHI_WRONG: library_maps(CACHE)}

    # m = 2 orientation (polar cells), data and library, relative to the true-angle patterns
    obs = {p: boot_errors(ssl.rotate_z(stars_mw, np.radians(p)), np.random.default_rng(SEED)) for p in lib}
    polar = {PHI_TRUE: np.load(ARR_DIR / "000.npz"), PHI_WRONG: np.load(CACHE)}
    orient = {}
    for f, name in FIELDS.items():
        Zd0, Zl0 = fields(obs[PHI_TRUE], polar[PHI_TRUE]["lib_m1"], polar[PHI_TRUE]["lib_s"], f)
        Zd, Zl = fields(obs[PHI_WRONG], polar[PHI_WRONG]["lib_m1"], polar[PHI_WRONG]["lib_s"], f)
        orient[name] = (rotation_deg(Zd, Zd0), rotation_deg(Zl, Zl0))
        print(f"{name:10s} trial {PHI_WRONG:.0f} deg: data pattern rotated {orient[name][0]:+5.1f} deg, "
              f"library pattern rotated {orient[name][1]:+5.1f} deg (expected data {PHI_WRONG - PHI_TRUE:+.0f})")

    rows = [("data, stars rotated by 27° (true angle)", "value", ("d", PHI_TRUE), None),
            ("library, trial 27° (ground truth)", "value", ("l", PHI_TRUE), None),
            ("data, stars rotated by 20° (wrong angle)", "value", ("d", PHI_WRONG), None),
            ("library, trial 20°", "value", ("l", PHI_WRONG), None),
            ("library 20° − library 27°", "diff", ("l", PHI_WRONG), ("l", PHI_TRUE)),
            ("data 20° − library 20°  (what χ² sees, wrong angle)", "diff", ("d", PHI_WRONG), ("l", PHI_WRONG)),
            ("data 27° − library 27°  (what χ² sees, true angle)", "diff", ("d", PHI_TRUE), ("l", PHI_TRUE))]
    get = lambda key: data[key[1]] if key[0] == "d" else lib[key[1]]
    e = ssl.BOX_XY
    fig, axes = plt.subplots(len(rows), 4, figsize=(15.5, 3.5 * len(rows)),
                             gridspec_kw=dict(hspace=0.45, wspace=0.42, top=1 - 0.9 / (3.5 * len(rows)), bottom=0.02,
                                              left=0.05, right=0.95))
    for i, (label, kind, a, b) in enumerate(rows):
        m = mask[a[1]] if a[0] == "d" else mask[PHI_TRUE if b is None else b[1]]
        for j, (qname, q, stat) in enumerate(QUANTS):
            ax = axes[i, j]
            img = get(a)[(stat, q)] if kind == "value" else get(a)[(stat, q)] - get(b)[(stat, q)]
            img = np.where(m, np.nan, img)
            if kind == "value":
                cmap, lo, hi = VALUE_LIMS[qname]
                t = f"{label}\n{qname}"
            else:
                cmap, lo, hi = CMAP_DIV, -20, 20
                t = f"{label}\n{qname}: rms {np.sqrt(np.nanmean(img**2)):.1f} km/s"
            im = ax.pcolormesh(e, e, img, cmap=grey_bad(cmap), vmin=lo, vmax=hi, rasterized=True)
            ax.plot([-5, 5], [0, 0], color=ACCENT, lw=1.1, ls=(0, (4, 2)))
            if a[0] == "d" and kind == "value":
                ang = np.radians(a[1] - PHI_TRUE)            # real bar direction in this frame
                ax.plot([-5 * np.cos(ang), 5 * np.cos(ang)], [-5 * np.sin(ang), 5 * np.sin(ang)], color=INK, lw=1.1)
            ax.set_aspect("equal")
            ax.set_xlim(-8, 8)
            ax.set_ylim(-8, 8)
            ax.tick_params(labelsize=7)
            ax.set_title(t, color=INK_2, fontsize=8.3)
            cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
            cb.set_label("km/s" if kind == "value" else "Δ km/s", fontsize=7.5)
            cb.ax.tick_params(labelsize=7)
            if j == 0:
                ax.set_ylabel("y [kpc]", fontsize=8)
            if i == len(rows) - 1:
                ax.set_xlabel("x [kpc] (bar frame)", fontsize=8)
    o = ", ".join(f"{k} {v[0]:+.0f}° / {v[1]:+.0f}°" for k, v in orient.items())
    fig.suptitle("Wrong bar angle: stars rotated by 20° instead of 27°, integrated with the potential's bar on x "
                 "(dashed); solid = real bar in that frame\n0.5 kpc pixels, |z| < 1 kpc, grey = < 20 stars; "
                 f"m = 2 rotation of the 20° patterns vs the 27° ones, data / library: {o}",
                 y=1 - 0.12 / (3.5 * len(rows)), fontsize=10.5, va="top")
    savefig(fig, "test1_nbody_square_maps_wrong_bar_angle_20deg_vs_27deg")


if __name__ == "__main__":
    main()
