"""Diagnostic: what bar length (S) and strength (bar_frac) change in the moment maps (N-body, Test 1).

Omega_b and phi_b at their true values. Models: truth (S=1, f=1), two points on the S-f ridge seen in the
global search, and single-parameter offsets (from the pilot). For each, the library means and dispersions of
v_R and v_phi in the folded bar-frame cells are compared with the observed ones:
  figure 1: midplane (|z| < 0.5 kpc) maps, observed values (top row) and pulls (model - observed) / err_boot
  figure 2: radial profiles of the m = 2 amplitude of <v_R> and <v_phi>, and of the azimuthal means of
            <v_phi>, sigma_R, sigma_phi (|z| < 1 kpc), plus the chi^2 change vs the truth per radius.
"""
import json

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import ROOT, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from plot_au18_velocity_maps import CMAP_DIV
from common import CMAP_SEQ_W
import steady_state_likelihood as ssl
from test1_nbody_moments import boot_errors, N_MIN
from test1_nbody_pilot import MOCK, OUT_DIR

TRUE_OM, TRUE_PHI = 26.02, 27.0
MODELS = [  # label, S, f, source (pilot index or None = integrate now)
    ("truth (S 1, f 1)", 1.0, 1.0, 0),
    ("ridge (S 0.82, f 1.50)", 0.82, 1.50, None),
    ("ridge (S 1.25, f 0.70)", 1.25, 0.70, None),
    ("f 0.7", 1.0, 0.7, 6),
    ("f 1.4", 1.0, 1.4, 7),
    ("S 0.8", 0.8, 1.0, 8),
    ("S 1.25", 1.25, 1.0, 9),
]
NR, NP, NZ = len(ssl.CELL_R) - 1, len(ssl.CELL_PHI) - 1, len(ssl.CELL_Z) - 1


def grid(a):
    return a.reshape(NR, NP, NZ, -1)


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        stars_mw = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    stars = ssl.rotate_z(stars_mw, np.radians(TRUE_PHI))
    N, m1, s, em, es = boot_errors(stars, np.random.default_rng(0))
    obs = dict(m1=grid(m1), s=grid(s), em=grid(em), es=grid(es), N=N.reshape(NR, NP, NZ))

    ssl.ACCURACY = 1e-5
    lib = {}
    for label, S, f, src in MODELS:
        if src is not None:
            a = np.load(OUT_DIR / f"pilot_{src:02d}.npz")
            lib[label] = dict(m1=grid(a["lib_m1"]), s=grid(a["lib_s"]))
        else:
            cache = OUT_DIR / f"diag_SF_S{S:g}_f{f:g}.npz"
            if not cache.exists():
                _, arr = ssl.evaluate(stars_mw, TRUE_OM, TRUE_PHI, S, f, galaxy="nbody", extended=True)
                np.savez(cache, lib_m1=arr["lib_m1"], lib_s=arr["lib_s"])
            a = np.load(cache)
            lib[label] = dict(m1=grid(a["lib_m1"]), s=grid(a["lib_s"]))

    Rc = 0.5 * (ssl.CELL_R[1:] + ssl.CELL_R[:-1])
    good = (obs["N"] >= N_MIN) & (Rc[:, None, None] >= 1.0)
    # chi2 (means + dispersions of v_R, v_phi, v_z; R >= 1) per model, and by radius and moment
    chi = {}
    for label in lib:
        pm = ((lib[label]["m1"] - obs["m1"]) / obs["em"])**2
        ps = ((lib[label]["s"] - obs["s"]) / obs["es"])**2
        pm[~good], ps[~good] = 0, 0
        chi[label] = dict(total=float(np.nansum(pm) + np.nansum(ps)),
                          by_R=np.nansum(pm, axis=(1, 2, 3)) + np.nansum(ps, axis=(1, 2, 3)),
                          by_moment={n: float(np.nansum(x[..., k])) for n, x, k in
                                     (("<vR>", pm, 0), ("<vphi>", pm, 1), ("<vz>", pm, 2),
                                      ("sR", ps, 0), ("sphi", ps, 1), ("sz", ps, 2))})
    t = chi[MODELS[0][0]]
    temp = np.sqrt(2 * 6 * int(good.sum()))      # tempering denominator sqrt(2 N_data), N_data = 6 x n_cells
    print(f"Delta logL = -(chi2 - chi2_truth) / (2 sqrt(2 N_data)), N_data = 6 x {int(good.sum())}, "
          f"sqrt(2 N_data) = {temp:.1f}   [means + dispersions of v_R, v_phi, v_z; R >= 1 kpc]")
    for label in lib:
        d = {k: -0.5 * (v - t["by_moment"][k]) / temp for k, v in chi[label]["by_moment"].items()}
        print(f"  {label:24s} total {-0.5 * (chi[label]['total'] - t['total']) / temp:+7.2f}   " +
              "  ".join(f"{k} {v:+6.2f}" for k, v in d.items()))

    # ---------- figure 1: midplane maps ----------
    phi_e = ssl.CELL_PHI
    R_e = ssl.CELL_R
    PP, RR = np.meshgrid(np.concatenate([phi_e, phi_e[1:] + np.pi]), R_e)
    X, Y = RR * np.cos(PP), RR * np.sin(PP)

    def full(a):          # folded (R, phi in [0, pi)) -> full disc by point symmetry
        return np.concatenate([a, a], axis=1)

    quants = [("<v_R>", "m1", 0, "em"), ("<v_phi>", "m1", 1, "em"), ("sigma_R", "s", 0, "es"),
              ("sigma_phi", "s", 1, "es")]
    fig, axes = plt.subplots(len(MODELS) + 1, 4, figsize=(15, 3.3 * (len(MODELS) + 1)),
                             gridspec_kw=dict(hspace=0.35, wspace=0.35))
    mid = 0
    for j, (name, key, k, ekey) in enumerate(quants):
        ax = axes[0, j]
        val = np.where(obs["N"][:, :, mid] >= N_MIN, obs[key][:, :, mid, k], np.nan)
        lim = (-40, 40) if name == "<v_R>" else (np.nanpercentile(val, 2), np.nanpercentile(val, 98))
        im = ax.pcolormesh(X, Y, full(val), cmap=CMAP_DIV if name == "<v_R>" else CMAP_SEQ_W, vmin=lim[0],
                           vmax=lim[1], shading="flat", rasterized=True)
        ax.set_title(f"observed {name}", color=INK_2, fontsize=9.5)
        fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1])).set_label("km/s", fontsize=8)
        for i, (label, S, f, _) in enumerate(MODELS, start=1):
            ax = axes[i, j]
            pull = (lib[label][key][:, :, mid, k] - obs[key][:, :, mid, k]) / obs[ekey][:, :, mid, k]
            pull = np.where(obs["N"][:, :, mid] >= N_MIN, pull, np.nan)
            im = ax.pcolormesh(X, Y, full(pull), cmap=CMAP_DIV, vmin=-3, vmax=3, shading="flat", rasterized=True)
            ax.set_title(f"{label}: pull, rms {np.sqrt(np.nanmean(pull**2)):.2f}", color=INK_2, fontsize=9)
            fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1])).set_label("(model − obs)/err", fontsize=7.5)
    for ax in axes.flat:
        ax.set_aspect("equal")
        ax.set_xlim(-7, 7)
        ax.set_ylim(-7, 7)
        ax.plot([-5, 5], [0, 0], color=ACCENT, lw=0.8, ls=(0, (4, 2)))
        ax.tick_params(labelsize=7)
    fig.suptitle("N-body Test 1, midplane (|z| < 0.5 kpc), bar frame (bar along x, dashed): observed moments "
                 "and pulls of the orbit library for different bar length S and strength f "
                 f"(Ω_b = {TRUE_OM}, φ_b = {TRUE_PHI}°)", y=0.995, fontsize=11)
    savefig(fig, "test1_nbody_bar_length_strength_polar_pull_maps")

    # ---------- figure 2: radial profiles ----------
    def fourier(a, w):     # m = 2 amplitude and phase-free azimuthal mean over phi bins, |z| < 1 kpc
        lay = a[:, :, :2]
        ww = w[:, :, :2]
        mean_phi = np.nansum(lay * ww, axis=2) / np.nansum(ww, axis=2)
        ph = 0.5 * (phi_e[1:] + phi_e[:-1])
        c = np.nanmean(mean_phi * np.exp(2j * ph)[None, :], axis=1)
        return 2 * np.abs(c), np.nanmean(mean_phi, axis=1)

    wobs = np.where(obs["N"] >= N_MIN, obs["N"], 0).astype(float)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), gridspec_kw=dict(hspace=0.3, wspace=0.28))
    panels = [("m=2 amplitude of <v_R>", "m1", 0, "a2"), ("m=2 amplitude of <v_phi>", "m1", 1, "a2"),
              ("azimuthal mean <v_phi>", "m1", 1, "m0"), ("azimuthal mean sigma_R", "s", 0, "m0"),
              ("azimuthal mean sigma_phi", "s", 1, "m0")]
    cols = [INK, ACCENT, "#e87ba4", BLUE_RAMP[2], BLUE_RAMP[5], "#1baf7a", "#8a8984"]
    styles = ["-", "-", "-", "--", "--", ":", ":"]
    for ax, (title, key, k, kind) in zip(axes.flat, panels):
        a2, m0 = fourier(obs[key][..., k], wobs)
        ax.plot(Rc, a2 if kind == "a2" else m0, "o", color=INK_MUTED, ms=5, label="observed")
        for (label, S, f, _), c, ls in zip(MODELS, cols, styles):
            a2, m0 = fourier(lib[label][key][..., k], wobs)
            ax.plot(Rc, a2 if kind == "a2" else m0, ls, color=c, lw=1.8 if label.startswith("truth") else 1.2,
                    label=label)
        ax.set_title(title + " (|z| < 1 kpc)", color=INK_2, fontsize=9.5)
        ax.set_xlabel("R [kpc]")
        ax.set_ylabel("km/s")
        ax.set_xlim(0, 8)
    axes[0, 0].legend(fontsize=7.5)
    ax = axes[1, 2]
    for (label, S, f, _), c, ls in zip(MODELS[1:], cols[1:], styles[1:]):
        ax.plot(Rc, -0.5 * (chi[label]["by_R"] - t["by_R"]) / temp, ls, color=c, lw=1.4, label=label)
    ax.axhline(0, color=INK_MUTED, lw=0.8)
    ax.set_title(f"Δlog L vs truth per radius (tempered by {temp:.1f}; all |z|, R ≥ 1)", color=INK_2, fontsize=9.5)
    ax.set_xlabel("R [kpc]")
    ax.set_xlim(0, 8)
    ax.legend(fontsize=7.5)
    fig.suptitle("N-body Test 1: what bar length S and strength f change (Ω_b, φ_b at truth)", y=0.97, fontsize=11)
    savefig(fig, "test1_nbody_bar_length_strength_radial_profiles")
    with open(OUT_DIR / "diag_SF_chi2.json", "w") as fh:
        json.dump(dict(tempering_sqrt_2N=float(temp),
                       models={l: dict(dlogL_tempered=-0.5 * (v["total"] - t["total"]) / temp,
                                       dlogL_untempered=-0.5 * (v["total"] - t["total"]),
                                       by_moment_tempered={k: -0.5 * (x - t["by_moment"][k]) / temp
                                                           for k, x in v["by_moment"].items()})
                               for l, v in chi.items()}), fh, indent=2)


if __name__ == "__main__":
    main()
