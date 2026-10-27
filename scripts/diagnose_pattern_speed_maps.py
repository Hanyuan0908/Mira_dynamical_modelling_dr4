"""Diagnostic: bar-frame <v_R>, <v_phi> maps of the Test 1 mock (observed, t = 0) vs the time-sampled
orbit library that the steady-state likelihood uses (t = 0.25-1 Gyr, 4-fold folded, weights 1/tau),
for the true and two wrong pattern speeds. Bar angle, strength and length at their true values.
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import MOCK_DIR, set_style, savefig, CMAP_SEQ_W, INK_2
from plot_au18_velocity_maps import CMAP_DIV
from plot_au18_vR_distance_error import quadrupole_node
import steady_state_likelihood as ssl

OMEGAS = (36.27, 45.0, 30.0)
LIM, PIX, ZCUT = 8.0, 0.5, 1.0


def maps(w, weight=None):
    s = np.abs(w[:, 2]) < ZCUT
    x, y, vx, vy = w[s, 0], w[s, 1], w[s, 3], w[s, 4]
    wt = np.ones(s.sum()) if weight is None else weight[s]
    R = np.hypot(x, y)
    vR, vphi = (x * vx + y * vy) / R, -(x * vy - y * vx) / R
    e = np.arange(-LIM, LIM + 1e-9, PIX)
    W = np.histogram2d(x, y, bins=[e, e], weights=wt)[0]
    N = np.histogram2d(x, y, bins=[e, e])[0]
    out = {}
    for k, q in (("vR", vR), ("vphi", vphi)):
        with np.errstate(invalid="ignore"):
            out[k] = np.where(N >= 20, np.histogram2d(x, y, bins=[e, e], weights=wt * q)[0] / W, np.nan).T
    return e, out


def main():
    set_style()
    with h5py.File(MOCK_DIR / "mock_perfect_R8_N50k.hdf5", "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    stars = ssl.rotate_z(stars_mw, np.radians(27.0))         # true bar -> x-axis
    pot = ssl.bar_aligned_potential(1.0, 1.0)
    e, obs = maps(stars)
    q_obs = quadrupole_node(stars[:, 0], stars[:, 1], stars[:, 2], stars[:, 3], stars[:, 4], np.ones(len(stars)))
    print(f"observed: v_R quadrupole amp {q_obs[1]:.1f} km/s, node {q_obs[0]:.1f} deg")

    rows = []
    for Om in OMEGAS:
        times, traj = ssl.integrate(pot, stars, -Om)
        lib, lib_oid, lib_w, tau = ssl.build_library(times, traj)
        _, m = maps(lib, lib_w)
        _, m_last = maps(traj[:, -1].astype(np.float64))         # single snapshot at t = 1 Gyr, unfolded
        ql = quadrupole_node(lib[:, 0], lib[:, 1], lib[:, 2], lib[:, 3], lib[:, 4], lib_w)
        qe = quadrupole_node(*(traj[:, -1, i].astype(np.float64) for i in (0, 1, 2, 3, 4)), np.ones(len(traj)))
        inner = np.hypot(*np.meshgrid(0.5 * (e[1:] + e[:-1]), 0.5 * (e[1:] + e[:-1]))) < 6
        rms = {k: float(np.sqrt(np.nanmean((m[k] - obs[k])[inner]**2))) for k in m}
        print(f"Omega_b = {Om:5.2f}: library v_R quadrupole amp {ql[1]:.1f} km/s node {ql[0]:.1f} deg; "
              f"t=1 Gyr snapshot amp {qe[1]:.1f} node {qe[0]:.1f}; rms(lib - obs, R<6) "
              f"v_R {rms['vR']:.1f}, v_phi {rms['vphi']:.1f} km/s; tau<1 for {np.mean(tau < 1):.1%} of orbits")
        rows.append((Om, m, m_last, ql, rms))
        del lib, lib_oid, lib_w, traj

    fig, axes = plt.subplots(1 + len(OMEGAS), 5, figsize=(21, 4.1 * (1 + len(OMEGAS))),
                             gridspec_kw=dict(wspace=0.38, hspace=0.3))
    kw_vR = dict(cmap=CMAP_DIV, vmin=-40, vmax=40)
    kw_vphi = dict(cmap=CMAP_SEQ_W, vmin=0, vmax=260)
    kw_d = dict(cmap=CMAP_DIV, vmin=-20, vmax=20)

    def show(ax, data, kw, title, label):
        im = ax.pcolormesh(e, e, data, rasterized=True, **kw)
        ax.plot([-5.75, 5.75], [0, 0], color="#eb6834", lw=1, ls=(0, (4, 2)))
        ax.set_aspect("equal")
        ax.set_title(title, color=INK_2, fontsize=9.5)
        cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
        cb.set_label(label, fontsize=8.5)

    show(axes[0, 0], obs["vR"], kw_vR, f"observed (t = 0) ⟨v_R⟩, quad. amp {q_obs[1]:.1f}", "km/s")
    show(axes[0, 1], obs["vphi"], kw_vphi, "observed (t = 0) ⟨v_φ⟩", "km/s")
    for ax in axes[0, 2:]:
        ax.axis("off")
    for r, (Om, m, m_last, ql, rms) in enumerate(rows, start=1):
        tag = f"Ω_b = {Om:g}" + (" (truth)" if Om == OMEGAS[0] else "")
        show(axes[r, 0], m["vR"], kw_vR, f"{tag}: library ⟨v_R⟩, amp {ql[1]:.1f}", "km/s")
        show(axes[r, 1], m["vR"] - obs["vR"], kw_d, f"library − observed ⟨v_R⟩ (rms {rms['vR']:.1f})", "km/s")
        show(axes[r, 2], m["vphi"], kw_vphi, f"{tag}: library ⟨v_φ⟩", "km/s")
        show(axes[r, 3], m["vphi"] - obs["vphi"], kw_d, f"library − observed ⟨v_φ⟩ (rms {rms['vphi']:.1f})", "km/s")
        show(axes[r, 4], m_last["vR"], kw_vR, f"{tag}: single snapshot t = 1 Gyr ⟨v_R⟩", "km/s")
    for ax in axes.flat:
        if ax.axison:
            ax.set_xlabel("x [kpc] (bar frame)")
    fig.suptitle("Test 1 mock in the bar frame (|z| < 1 kpc, 0.5 kpc pixels): observed vs orbit library "
                 "(t = 0.25–1 Gyr, 4-fold folded, 1/τ weights); dashed = bar major axis", y=0.995, fontsize=12)
    savefig(fig, "test1_omega_maps")


if __name__ == "__main__":
    main()
