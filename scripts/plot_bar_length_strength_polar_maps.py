"""Face-on maps of <v_R>, <v_phi>, sigma_R, sigma_phi for the observed N-body mock and the orbit libraries of
models with different bar length S and strength f (Omega_b, phi_b at the truth). No new integrations: library
cell moments come from the pilot (results/test1_nbody/pilot_NN.npz) and diagnose_bar_length_strength.py
(diag_SF_*.npz). Cells: folded bar-frame polar cells (R 0.5 kpc x phi 15 deg), midplane layer |z| < 0.5 kpc.

  test1_nbody_bar_length_strength_polar_maps_values.png    actual values
  test1_nbody_bar_length_strength_polar_maps_nonaxisymmetric.png    value minus its azimuthal (ring) mean: the non-axisymmetric part
  test1_nbody_bar_length_strength_polar_maps_model_minus_truth.png      model minus the truth model (library noise only, ~1-2 km/s)
  test1_nbody_bar_length_strength_m2_profiles.png        m = 2 amplitude and phase vs R of all four quantities
Also prints rms differences and m = 2 amplitudes (numbers in results/test1_nbody/SF_maps_numbers.json).
"""
import json

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP, CMAP_SEQ_W
from plot_au18_velocity_maps import CMAP_DIV
import steady_state_likelihood as ssl
from test1_nbody_pilot import MOCK, OUT_DIR
from diagnose_bar_length_strength import MODELS, TRUE_PHI

NR, NP, NZ = len(ssl.CELL_R) - 1, len(ssl.CELL_PHI) - 1, len(ssl.CELL_Z) - 1
QUANTS = [("<v_R>", "m1", 0), ("<v_phi>", "m1", 1), ("sigma_R", "s", 0), ("sigma_phi", "s", 1)]
LAYER = 0
N_MIN = 20


def load_models():
    out = {}
    for label, S, f, src in MODELS:
        fn = OUT_DIR / (f"pilot_{src:02d}.npz" if src is not None else f"diag_SF_S{S:g}_f{f:g}.npz")
        a = np.load(fn)
        out[label] = {"m1": a["lib_m1"].reshape(NR, NP, NZ, 3), "s": a["lib_s"].reshape(NR, NP, NZ, 3)}
    return out


def ring_mean(a):
    return np.nanmean(a, axis=1, keepdims=True)


def m2(a):
    """m = 2 amplitude and phase (deg, of the maximum, in [0, 90)) per ring of a (NR, NP) map."""
    ph = 0.5 * (ssl.CELL_PHI[1:] + ssl.CELL_PHI[:-1])
    c = np.nanmean(a * np.exp(2j * ph)[None, :], axis=1)
    return 2 * np.abs(c), np.mod(np.degrees(0.5 * np.angle(c)), 90.0) if False else np.degrees(0.5 * np.angle(c))


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        stars = ssl.rotate_z(np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")]),
                             np.radians(TRUE_PHI))
    N, m1, s = ssl.cell_moments(stars)
    obs = {"m1": m1.reshape(NR, NP, NZ, 3), "s": s.reshape(NR, NP, NZ, 3), "N": N.reshape(NR, NP, NZ)}
    mods = load_models()
    labels = [m[0] for m in MODELS]
    truth = labels[0]

    PP, RR = np.meshgrid(np.concatenate([ssl.CELL_PHI, ssl.CELL_PHI[1:] + np.pi]), ssl.CELL_R)
    X, Y = RR * np.cos(PP), RR * np.sin(PP)

    def full(a):
        return np.concatenate([a, a], axis=1)

    def get(src, key, k, layer=LAYER):
        return src[key][:, :, layer, k]

    obs_mask = obs["N"][:, :, LAYER] < N_MIN

    def panel(fig, ax, data, cmap, vmin, vmax, title, label="km/s"):
        im = ax.pcolormesh(X, Y, full(data), cmap=cmap, vmin=vmin, vmax=vmax, shading="flat", rasterized=True)
        ax.plot([-5, 5], [0, 0], color=ACCENT, lw=0.8, ls=(0, (4, 2)))
        ax.set_aspect("equal")
        ax.set_xlim(-8, 8)
        ax.set_ylim(-8, 8)
        ax.tick_params(labelsize=7)
        ax.set_title(title, color=INK_2, fontsize=8.5)
        fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1])).set_label(label, fontsize=7.5)

    rows = [("observed", obs)] + [(l, mods[l]) for l in labels]

    # ---------- 1. values ----------
    lims = {"<v_R>": (-40, 40), "<v_phi>": (0, 240), "sigma_R": (40, 200), "sigma_phi": (30, 190)}
    fig, axes = plt.subplots(len(rows), 4, figsize=(15, 3.2 * len(rows)), gridspec_kw=dict(hspace=0.35, wspace=0.35))
    for i, (lab, src) in enumerate(rows):
        for j, (name, key, k) in enumerate(QUANTS):
            d = get(src, key, k)
            if lab == "observed":
                d = np.where(obs_mask, np.nan, d)
            panel(fig, axes[i, j], d, CMAP_DIV if name == "<v_R>" else CMAP_SEQ_W, *lims[name], f"{lab}: {name}")
    fig.suptitle("Values: midplane (|z| < 0.5 kpc), bar frame, bar along x (dashed); Ω_b = 26.02, φ_b = 27°",
                 y=0.995, fontsize=11)
    savefig(fig, "test1_nbody_bar_length_strength_polar_maps_values")

    # ---------- 2. non-axisymmetric part ----------
    fig, axes = plt.subplots(len(rows), 4, figsize=(15, 3.2 * len(rows)), gridspec_kw=dict(hspace=0.35, wspace=0.35))
    for i, (lab, src) in enumerate(rows):
        for j, (name, key, k) in enumerate(QUANTS):
            d = get(src, key, k)
            if lab == "observed":
                d = np.where(obs_mask, np.nan, d)
            panel(fig, axes[i, j], d - ring_mean(d), CMAP_DIV, -30, 30, f"{lab}: {name} − ring mean")
    fig.suptitle("Non-axisymmetric part (value minus its azimuthal mean in each ring), midplane", y=0.995, fontsize=11)
    savefig(fig, "test1_nbody_bar_length_strength_polar_maps_nonaxisymmetric")

    # ---------- 3. model minus truth model ----------
    fig, axes = plt.subplots(len(labels) - 1, 4, figsize=(15, 3.2 * (len(labels) - 1)),
                             gridspec_kw=dict(hspace=0.35, wspace=0.35))
    numbers = {"rms_model_minus_truth_kms": {}, "typical_data_error_kms": {}, "m2_amplitude_kms": {}}
    Rc = 0.5 * (ssl.CELL_R[1:] + ssl.CELL_R[:-1])
    sel = (obs["N"][:, :, :2] >= N_MIN) & (Rc[:, None, None] >= 1)          # |z| < 1, R >= 1, populated
    for j, (name, key, k) in enumerate(QUANTS):
        err = (obs["s"][..., k] / np.sqrt(np.maximum(obs["N"], 1)))[:, :, :2]
        if key == "s":
            err = err / np.sqrt(2)
        numbers["typical_data_error_kms"][name] = float(np.median(err[sel]))
    for i, lab in enumerate(labels[1:]):
        numbers["rms_model_minus_truth_kms"][lab] = {}
        for j, (name, key, k) in enumerate(QUANTS):
            d = get(mods[lab], key, k) - get(mods[truth], key, k)
            panel(fig, axes[i, j], d, CMAP_DIV, -10, 10, f"{lab} − truth: {name}  (rms {np.sqrt(np.nanmean(d**2)):.1f})")
            dd = (mods[lab][key][..., k] - mods[truth][key][..., k])[:, :, :2]
            numbers["rms_model_minus_truth_kms"][lab][name] = float(np.sqrt(np.mean(dd[sel]**2)))
    fig.suptitle("Model minus truth model (both libraries, so no data noise), midplane; "
                 "typical data error per cell: " + ", ".join(f"{k} {v:.1f}" for k, v in numbers["typical_data_error_kms"].items())
                 + " km/s", y=0.995, fontsize=10.5)
    savefig(fig, "test1_nbody_bar_length_strength_polar_maps_model_minus_truth")

    # ---------- 4. m = 2 profiles of all four quantities (|z| < 0.5 kpc) ----------
    cols = [INK, ACCENT, "#e87ba4", BLUE_RAMP[2], BLUE_RAMP[5], "#1baf7a", "#8a8984"]
    styles = ["-", "-", "-", "--", "--", ":", ":"]
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), gridspec_kw=dict(hspace=0.3, wspace=0.3))
    for j, (name, key, k) in enumerate(QUANTS):
        d = np.where(obs_mask, np.nan, get(obs, key, k))
        a, p = m2(d)
        axes[0, j].plot(Rc, a, "o", color=INK_MUTED, ms=5, label="observed")
        axes[1, j].plot(Rc, p, "o", color=INK_MUTED, ms=5)
        numbers["m2_amplitude_kms"].setdefault("observed", {})[name] = dict(zip(map(str, Rc), map(float, a)))
        for lab, c, ls in zip(labels, cols, styles):
            a, p = m2(get(mods[lab], key, k))
            axes[0, j].plot(Rc, a, ls, color=c, lw=1.8 if lab == truth else 1.2, label=lab)
            axes[1, j].plot(Rc, p, ls, color=c, lw=1.8 if lab == truth else 1.2)
            numbers["m2_amplitude_kms"].setdefault(lab, {})[name] = dict(zip(map(str, Rc), map(float, a)))
        axes[0, j].set_title(f"m = 2 amplitude of {name} (|z| < 0.5)", color=INK_2, fontsize=9.5)
        axes[1, j].set_title(f"m = 2 phase of {name} [deg]", color=INK_2, fontsize=9.5)
        for ax in axes[:, j]:
            ax.set_xlabel("R [kpc]")
            ax.set_xlim(0, 8)
        axes[1, j].set_ylim(-90, 90)
    axes[0, 0].legend(fontsize=7.5)
    fig.suptitle("m = 2 Fourier profiles of the mean and dispersion maps (phase 0 = bar major axis)", y=0.97, fontsize=11)
    savefig(fig, "test1_nbody_bar_length_strength_m2_profiles")

    with open(OUT_DIR / "SF_maps_numbers.json", "w") as fh:
        json.dump(numbers, fh, indent=2)
    print("typical data error per cell (|z|<1, R>=1, N>=20) [km/s]:", {k: round(v, 1) for k, v in numbers["typical_data_error_kms"].items()})
    print("rms(model - truth model) over the same cells [km/s]:")
    for lab, v in numbers["rms_model_minus_truth_kms"].items():
        print(f"  {lab:24s} " + "  ".join(f"{k} {x:5.2f}" for k, x in v.items()))
    print("m = 2 amplitude at R = 2.25 / 3.25 / 4.25 / 5.25 kpc (|z|<0.5) [km/s]:")
    for lab in ["observed"] + labels:
        print(f"  {lab:24s} " + "  ".join(f"{name} " + "/".join(f"{numbers['m2_amplitude_kms'][lab][name][str(r)]:4.1f}" for r in (2.25, 3.25, 4.25, 5.25))
                                          for name, _, _ in QUANTS))


if __name__ == "__main__":
    main()
