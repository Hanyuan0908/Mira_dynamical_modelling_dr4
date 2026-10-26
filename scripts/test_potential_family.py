"""Checks of potential_family.make_potential (run after potential_family.py has built its files).

1. bar_frac = S = 1 reproduces the ground truth (Au18_full.ini) at t = 0 and, rotated, at t > 0.
2. The azimuthally averaged potential (rotation curve) is unchanged for any bar_frac and S.
3. bar_frac = 0 gives an axisymmetric potential.
4. The non-axisymmetric part scales as bar_frac * Phi_nonaxi(x / S).
5. The bar rotates by Omega_b t (clockwise for Omega_b < 0).
6. An orbit integrated in make_potential (inertial frame) and mapped with to_bar_frame equals the
   same orbit integrated in the static potential with agama.orbit(..., Omega=Omega_b).
"""
import json

import numpy as np
import agama
import matplotlib.pyplot as plt

from common import ROOT, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from build_au18_potential import vcirc
from plot_au18_velocity_maps import CMAP_DIV
from potential_family import (FILES, OMEGA_B_TRUE, TIME_UNIT_GYR, make_potential, load_axi,
                              to_bar_frame)

RESULTS_FILE = ROOT / "results" / "potential_family_checks.json"
PHI = np.linspace(0, 2 * np.pi, 72, endpoint=False)


def ring(R, z=0.0):
    return np.column_stack([R * np.cos(PHI), R * np.sin(PHI), np.full_like(PHI, z)])


def m2(pot, R, t=0.0):
    """Amplitude and phase (deg, potential maximum) of the m = 2 term of Phi on a ring."""
    c = np.mean(pot.potential(ring(R), t=t) * np.exp(2j * PHI))
    return 2 * np.abs(c), np.degrees(0.5 * np.angle(c))


def rot_z(x, a):
    c, s = np.cos(a), np.sin(a)
    return np.column_stack([c * x[:, 0] - s * x[:, 1], s * x[:, 0] + c * x[:, 1], x[:, 2]])


def main():
    set_style()
    res = {}
    truth = agama.Potential(file=str(FILES["full"]))
    axi = load_axi()
    fam = make_potential(pot_axi=axi)
    pts = np.random.default_rng(3).uniform(-15, 15, (5000, 3))
    pts[:, 2] /= 5

    # 1. truth at t = 0 and t > 0
    res["truth_t0_max_rel_dphi"] = float(np.max(np.abs(fam.potential(pts) / truth.potential(pts) - 1)))
    res["truth_t0_max_rel_dforce"] = float(np.max(np.abs(fam.force(pts) - truth.force(pts))) /
                                           np.max(np.abs(truth.force(pts))))
    t1 = 0.05
    res["truth_t0.05_max_rel_dphi"] = float(np.max(np.abs(
        fam.potential(pts, t=t1) / truth.potential(rot_z(pts, -OMEGA_B_TRUE * t1)) - 1)))

    # 2-4. rotation curve, axisymmetric limit, scaling
    R = np.linspace(0.5, 20, 40)
    vc_truth = vcirc(truth, R)
    members = [(0.0, 1.0), (0.5, 1.0), (2.0, 1.0), (1.0, 0.8), (1.0, 1.25)]
    res["vc_max_rel_change"] = {}
    for f, S in members:
        p = make_potential(f, S, pot_axi=axi)
        res["vc_max_rel_change"][f"bar_frac={f},S={S}"] = float(np.max(np.abs(vcirc(p, R) / vc_truth - 1)))
    p0 = make_potential(0.0, 1.0, pot_axi=axi)
    res["bar_frac0_max_phi_variation"] = float(max(np.ptp(p0.potential(ring(r))) for r in (1, 3, 6, 10)))
    nonaxi = lambda p, x: p.potential(x) - axi.potential(x)
    pf = make_potential(1.7, 1.3, pot_axi=axi)
    xs = pts[np.hypot(pts[:, 0], pts[:, 1]) < 10]
    ref = 1.7 * nonaxi(fam, xs / 1.3)
    res["scaling_max_rel_err"] = float(np.max(np.abs(nonaxi(pf, xs) - ref)) / np.max(np.abs(ref)))

    # 5. rotation sense
    ts = np.linspace(0, 0.1, 6)
    ph = np.array([m2(fam, 3.0, t)[1] for t in ts])
    expect = ph[0] + np.degrees(OMEGA_B_TRUE * ts)
    res["rotation_max_phase_err_deg"] = float(np.max(np.abs((ph - expect + 90) % 180 - 90)))

    # 6. inertial integration with make_potential vs rotating-frame integration
    static = make_potential(Omega_b=0.0, pot_axi=axi)
    ic = np.array([[3.0, 0.5, 0.2, 20.0, 180.0, 5.0], [6.0, -2.0, 0.5, -60.0, 220.0, -10.0],
                   [1.0, 1.0, -0.3, 120.0, 40.0, 30.0]])
    T = 0.5
    orb_in = agama.orbit(potential=fam, ic=ic, time=T, trajsize=2)
    orb_rot = agama.orbit(potential=static, ic=ic, time=T, trajsize=2, Omega=OMEGA_B_TRUE)
    end_in = to_bar_frame(np.array([o[1][-1] for o in orb_in]), T)
    end_rot = np.array([o[1][-1] for o in orb_rot])
    res["orbit_frames_max_dpos_kpc"] = float(np.max(np.abs(end_in[:, :3] - end_rot[:, :3])))
    res["orbit_frames_max_dvel_kms"] = float(np.max(np.abs(end_in[:, 3:] - end_rot[:, 3:])))

    with open(RESULTS_FILE, "w") as fh:
        json.dump(res, fh, indent=2)
    print(json.dumps(res, indent=1))

    # figure
    fig = plt.figure(figsize=(15, 9.4))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.5)
    e = np.linspace(-8, 8, 161)
    X, Y = np.meshgrid(0.5 * (e[1:] + e[:-1]), 0.5 * (e[1:] + e[:-1]))
    grid = np.column_stack([X.ravel(), Y.ravel(), 0 * X.ravel()])
    lim = np.max(np.abs(nonaxi(fam, grid)))
    for j, (f, S, t, title) in enumerate([(1.0, 1.0, 0.0, "truth: bar_frac = 1, S = 1, t = 0"),
                                          (1.5, 1.0, 0.0, "bar_frac = 1.5, S = 1"),
                                          (1.0, 1.25, 0.05, f"bar_frac = 1, S = 1.25, t = {0.05 * TIME_UNIT_GYR * 1e3:.0f} Myr")]):
        p = make_potential(f, S, pot_axi=axi)
        ax = fig.add_subplot(gs[0, j])
        val = (p.potential(grid, t=t) - axi.potential(grid)).reshape(X.shape)
        im = ax.pcolormesh(e, e, val, cmap=CMAP_DIV, vmin=-lim, vmax=lim, rasterized=True)
        ax.contour(X, Y, val, levels=np.linspace(-lim, lim, 13), colors=INK_2, linewidths=0.4)
        ang = np.radians(-27) + OMEGA_B_TRUE * t
        ax.plot([-6 * np.cos(ang), 6 * np.cos(ang)], [-6 * np.sin(ang), 6 * np.sin(ang)], color=ACCENT,
                lw=1.2, ls=(0, (4, 2)))
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$ [kpc]")
        if j == 0:
            ax.set_ylabel(r"$y$ [kpc]")
        ax.set_title(title, color=INK_2, fontsize=9.5)
        cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
        cb.set_label(r"$\Phi - \Phi_{\rm axi}$ [(km s$^{-1}$)$^2$]", fontsize=8.5)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(R, vc_truth, color=INK, lw=3, label="truth")
    cols = [INK_MUTED, BLUE_RAMP[2], BLUE_RAMP[5], ACCENT, "#1baf7a"]
    for (f, S), col in zip(members, cols):
        ax.plot(R, vcirc(make_potential(f, S, pot_axi=axi), R), color=col, lw=1.2, ls="--",
                label=f"bar_frac={f}, S={S}")
    ax.set_xlabel(r"$R$ [kpc]")
    ax.set_ylabel(r"$v_c$ [km s$^{-1}$]")
    ax.set_title("Rotation curve is unchanged across the family", color=INK_2, fontsize=9.5)
    ax.legend(fontsize=7.5)

    ax = fig.add_subplot(gs[1, 1])
    Rr = np.linspace(0.3, 15, 60)
    for (f, S), col in zip([(1.0, 1.0)] + members[1:], [INK] + cols[1:]):
        p = make_potential(f, S, pot_axi=axi)
        ax.plot(Rr, [m2(p, r)[0] for r in Rr], color=col, lw=2 if (f, S) == (1.0, 1.0) else 1.3,
                label=f"bar_frac={f}, S={S}")
    ax.set_xlabel(r"$R$ [kpc]")
    ax.set_ylabel(r"$|\Phi_2|(R)$ [(km s$^{-1}$)$^2$]")
    ax.set_title("m = 2 amplitude: bar_frac scales it, S stretches it", color=INK_2, fontsize=9.5)
    ax.legend(fontsize=7.5)

    ax = fig.add_subplot(gs[1, 2])
    tt = np.linspace(0, 0.2, 41)
    phs = np.array([m2(fam, 3.0, t)[1] for t in tt])
    ax.plot(tt * TIME_UNIT_GYR * 1e3, np.degrees(np.unwrap(np.radians(2 * phs))) / 2 - 90, "o", color=BLUE_RAMP[4],
            ms=4, label="bar angle from m = 2 phase of Φ")
    ax.plot(tt * TIME_UNIT_GYR * 1e3, -27 + np.degrees(OMEGA_B_TRUE * tt), color=ACCENT, lw=1.5,
            label=rf"$-27° + \Omega_b t$, $\Omega_b$ = {OMEGA_B_TRUE:.2f} km s$^{{-1}}$ kpc$^{{-1}}$")
    ax.set_xlabel("t [Myr]")
    ax.set_ylabel("bar angle [deg]")
    ax.set_title("Bar rotates clockwise at the true pattern speed", color=INK_2, fontsize=9.5)
    ax.legend(fontsize=7.5)
    fig.suptitle("Potential family: Φ_axi + bar_frac·[Φ_bary(R(−ψ) x/S) − Φ_bary,axi(x/S)],  "
                 "ψ = Ω_b t  (dashed: bar major axis)", y=0.97, fontsize=11)
    savefig(fig, "potential_family_check")


if __name__ == "__main__":
    main()
