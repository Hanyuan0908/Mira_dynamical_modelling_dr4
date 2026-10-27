"""Diagnostic: v_R, v_phi, v_z distributions of the Test 1 mock (observed, t = 0) vs the orbit library
(t = 0.25-1 Gyr, 4-fold folded, 1/tau weights) for three pattern speeds, in cells along the bar major
axis (x = 1..6 kpc, y = 0) and along the 45-degree line (where the v_R quadrupole peaks).
Cells: cylinders of radius R_CELL around (x, y) and its point-symmetric image, |z| < Z_CELL; the data
are folded the same way as the library. Each panel is annotated with the summed per-star
Delta log p(v|x) of the stars in the cell (from results/test1_qmc/per_star_logp_omega.npz).
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import ROOT, MOCK_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl

OMEGAS = (36.27, 45.0, 30.0)
COLS = {36.27: INK, 45.0: ACCENT, 30.0: BLUE_RAMP[3]}
R_CELL, Z_CELL = 0.45, 0.5
RADII = (1, 2, 3, 4, 5, 6)
LINES = {"bar major axis (y = 0)": 0.0, "45° from the bar": 45.0}
BINS = {"vR": np.arange(-400, 401, 20), "vphi": np.arange(-250, 451, 20), "vz": np.arange(-300, 301, 20)}
LABELS = {"vR": r"$v_R$", "vphi": r"$v_\phi$", "vz": r"$v_z$"}


def cyl(w):
    R = np.hypot(w[:, 0], w[:, 1])
    return dict(vR=(w[:, 0] * w[:, 3] + w[:, 1] * w[:, 4]) / R,
                vphi=-(w[:, 0] * w[:, 4] - w[:, 1] * w[:, 3]) / R, vz=w[:, 5] * np.sign(w[:, 2]))


def in_cell(w, x0, y0):
    """Cylinder around (x0, y0) or its point-symmetric image (-x0, -y0), |z| < Z_CELL."""
    d1 = np.hypot(w[:, 0] - x0, w[:, 1] - y0)
    d2 = np.hypot(w[:, 0] + x0, w[:, 1] + y0)
    return (np.minimum(d1, d2) < R_CELL) & (np.abs(w[:, 2]) < Z_CELL)


def moments(v, w=None):
    w = np.ones_like(v) if w is None else w
    m = np.average(v, weights=w)
    s = np.sqrt(np.average((v - m)**2, weights=w))
    return m, s, np.average(((v - m) / s)**3, weights=w), np.average(((v - m) / s)**4, weights=w) - 3


def main():
    set_style()
    with h5py.File(MOCK_DIR / "mock_perfect_R8_N50k.hdf5", "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    stars = ssl.rotate_z(stars_mw, np.radians(27.0))
    logp = np.load(ROOT / "results" / "test1_qmc" / "per_star_logp_omega.npz")
    lp = {Om: logp[f"Om{Om:g}"] for Om in OMEGAS}
    cells = [(name, r, r * np.cos(np.radians(a)), r * np.sin(np.radians(a))) for name, a in LINES.items()
             for r in RADII]

    # library samples inside each cell, for each pattern speed
    lib_cells = {}
    pot = ssl.bar_aligned_potential(1.0, 1.0)
    for Om in OMEGAS:
        times, traj = ssl.integrate(pot, stars, -Om)
        lib, _, lib_w, _ = ssl.build_library(times, traj)
        for c in cells:
            s = in_cell(lib, c[2], c[3])
            lib_cells[(Om, c)] = (cyl(lib[s]), lib_w[s])
        del lib, lib_w, traj

    for name in LINES:
        cs = [c for c in cells if c[0] == name]
        fig, axes = plt.subplots(len(cs), 3, figsize=(14, 2.55 * len(cs)), gridspec_kw=dict(hspace=0.55, wspace=0.22))
        for r, c in enumerate(cs):
            s = in_cell(stars, c[2], c[3])
            obs = cyl(stars[s])
            dlp = {Om: float(np.sum(lp[Om][s] - lp[36.27][s])) for Om in OMEGAS[1:]}
            for j, k in enumerate(("vR", "vphi", "vz")):
                ax = axes[r, j]
                b = BINS[k]
                h, _ = np.histogram(obs[k], b, density=True)
                err = np.sqrt(np.histogram(obs[k], b)[0]) / (s.sum() * np.diff(b))
                bc = 0.5 * (b[1:] + b[:-1])
                ax.fill_between(bc, h - err, h + err, step="mid", color=INK_MUTED, alpha=0.25, lw=0)
                ax.step(bc, h, where="mid", color=INK_MUTED, lw=1.2)
                txt = [f"obs: μ={moments(obs[k])[0]:.0f} σ={moments(obs[k])[1]:.0f}"]
                for Om in OMEGAS:
                    v, wt = lib_cells[(Om, c)]
                    hl, _ = np.histogram(v[k], b, weights=wt, density=True)
                    ax.step(bc, hl, where="mid", color=COLS[Om], lw=1.6 if Om == 36.27 else 1.2)
                    m = moments(v[k], wt)
                    txt.append(f"Ω={Om:g}: μ={m[0]:.0f} σ={m[1]:.0f} s={m[2]:+.2f} k={m[3]:+.2f}")
                ax.text(0.02, 0.97, "\n".join(txt), transform=ax.transAxes, fontsize=6.3, va="top",
                        color=INK_2, family="monospace")
                ax.set_xlabel(LABELS[k] + r" [km s$^{-1}$]", fontsize=8.5)
                ax.tick_params(labelsize=7.5)
                if j == 0:
                    ax.set_ylabel(f"R={c[1]} kpc\n(x,y)=({c[2]:.1f},{c[3]:.1f})", fontsize=8.5)
                if j == 2:
                    ax.text(0.98, 0.97, f"N_obs = {s.sum()}\nΣΔlog p (45 − true) = {dlp[45.0]:+.1f}\n"
                            f"ΣΔlog p (30 − true) = {dlp[30.0]:+.1f}", transform=ax.transAxes, fontsize=7,
                            ha="right", va="top", color=INK)
        handles = [plt.Line2D([], [], color=INK_MUTED, lw=6, alpha=0.4, label="observed (t=0, ±1σ Poisson)")]
        handles += [plt.Line2D([], [], color=COLS[Om], lw=1.6, label=f"library, Ω_b = {Om:g}" + (" (truth)" if Om == 36.27 else ""))
                    for Om in OMEGAS]
        fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=9, bbox_to_anchor=(0.5, 0.995))
        fig.suptitle(f"Velocity distributions along the {name}: cylinders r = {R_CELL} kpc, |z| < {Z_CELL} kpc "
                     "(4-fold folded); s = skewness, k = excess kurtosis", y=1.015, fontsize=11)
        savefig(fig, "test1_vdist_" + ("major_axis" if "major" in name else "45deg"))


if __name__ == "__main__":
    main()
