"""Square-pixel face-on maps (bar frame, 0.5 kpc pixels, |z| < 1 kpc) of <v_R>, <v_phi>, sigma_R, sigma_phi for the
N-body mock data and for orbit libraries at the true Omega_b and phi_b with different bar length / strength.

  test1_nbody_maps_length.png    rows: data; model at ground truth, S = 1.25, S = 0.8; bootstrap data error per pixel;
                                 each model - data; S = 1.25 - ground-truth model; S = 0.8 - ground-truth model
  test1_nbody_maps_strength.png  same for bar_frac = 1.4 and 0.7
Data are folded with the library's point symmetry ((x, y) -> (-x, -y); v_R, v_phi unchanged), so a pixel holds the stars
of itself and of its mirror pixel; pixels with fewer than N_MIN such stars are grey. Data errors: bootstrap over stars.
Library maps come from the saved 0.5 kpc box sums (run_SF_box_models.py; truth from the global search point 000).

Per panel (all over the data pixels):
  ground truth - data:  rms and chi^2 / pixel with the bootstrap errors
  model - data rows:    rms, Delta chi^2 = chi^2(model) - chi^2(ground truth) and the tempered
                        Delta logL = -Delta chi^2 / (2 sqrt(2 N_data)), N_data = 4 maps x n_pixels
  model - truth rows:   rms and the pixel-scale library noise in that difference (neighbour-difference estimate);
                        sum((model - truth)^2 / err^2) is saved too but includes that library noise
This 2D |z| < 1 chi^2 is a diagnostic of these maps only; the corner plots use the 3D polar-cell likelihood.
Numbers: results/test1_nbody/SF_square_maps_numbers.json.
"""
import json

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK_2, ACCENT, CMAP_SEQ_W
from plot_au18_velocity_maps import CMAP_DIV
import steady_state_likelihood as ssl
from test1_nbody_qmc import ARR_DIR, MOCK, TRUTH, OUT_CSV

Z_CUT = 1.0
N_MIN = 20
N_BOOT = 200
MASK_GREY = "#cfcfcc"
QUANTS = [("<v_R>", 0, "mean"), ("<v_phi>", 1, "mean"), ("sigma_R", 0, "disp"), ("sigma_phi", 1, "disp")]
VALUE_LIMS = {"<v_R>": (CMAP_DIV, -40, 40), "<v_phi>": (CMAP_SEQ_W, 0, 240), "sigma_R": (CMAP_SEQ_W, 40, 200),
              "sigma_phi": (CMAP_SEQ_W, 30, 190)}
PAIR_INDEX = {0: 0, 1: 1}                       # (RR, pp) positions in ssl.PAIRS


def grey_bad(cmap):
    c = cmap.copy()
    c.set_bad(MASK_GREY)
    return c


def data_maps(stars):
    """Folded data maps (N, mean, dispersion of v_R and v_phi) on the box grid; stars already cut to |z| < Z_CUT."""
    R = np.hypot(stars[:, 0], stars[:, 1])
    v = np.column_stack([(stars[:, 0] * stars[:, 3] + stars[:, 1] * stars[:, 4]) / R,
                         -(stars[:, 0] * stars[:, 4] - stars[:, 1] * stars[:, 3]) / R])
    x = np.concatenate([stars[:, 0], -stars[:, 0]])
    y = np.concatenate([stars[:, 1], -stars[:, 1]])
    v = np.concatenate([v, v])
    e = ssl.BOX_XY
    N = np.histogram2d(x, y, bins=[e, e])[0]
    out = {"N": N.T}                              # stars in the pixel plus its mirror pixel
    for a in range(2):
        s1 = np.histogram2d(x, y, bins=[e, e], weights=v[:, a])[0]
        s2 = np.histogram2d(x, y, bins=[e, e], weights=v[:, a]**2)[0]
        with np.errstate(invalid="ignore", divide="ignore"):
            m = s1 / N
            out[("mean", a)] = m.T
            out[("disp", a)] = np.sqrt(np.maximum(s2 / N - m**2, 0)).T
    return out


def boot_errors(stars, rng):
    reps = [data_maps(stars[rng.integers(0, len(stars), len(stars))]) for _ in range(N_BOOT)]
    return {(stat, q): np.nanstd([r[(stat, q)] for r in reps], axis=0) for _, q, stat in QUANTS}


def library_maps(npz):
    a = np.load(npz)
    zc = 0.5 * (ssl.BOX_Z[1:] + ssl.BOX_Z[:-1])
    kz = np.abs(zc) < Z_CUT
    S0 = a["box_S0"][:, :, kz].sum(2)
    S1 = a["box_S1"][:, :, kz].sum(2)
    S2 = a["box_S2"][:, :, kz].sum(2)
    out = {}
    with np.errstate(invalid="ignore", divide="ignore"):
        for q in range(2):
            m = S1[..., q] / S0
            out[("mean", q)] = m.T
            out[("disp", q)] = np.sqrt(np.maximum(S2[..., PAIR_INDEX[q]] / S0 - m**2, 0)).T
    return out


def figure(name, title, rows, data, err, lib, mask, temper, numbers):
    e = ssl.BOX_XY
    fig, axes = plt.subplots(len(rows), 4, figsize=(15.5, 3.45 * len(rows)),
                             gridspec_kw=dict(hspace=0.5, wspace=0.42, top=1 - 0.95 / (3.45 * len(rows)), bottom=0.025,
                                              left=0.05, right=0.95))
    for i, (label, kind, a, b) in enumerate(rows):
        for j, (qname, q, stat) in enumerate(QUANTS):
            ax = axes[i, j]
            key = (stat, q)
            t2 = None
            if kind == "value":
                img = (data if a == "data" else lib[a])[key]
                cmap, lo, hi = VALUE_LIMS[qname]
                unit = "km/s"
            elif kind == "error":
                img = err[key]
                cmap, lo, hi, unit = CMAP_SEQ_W, 0, 20, "km/s"
                t2 = f"median {np.nanmedian(np.where(mask, np.nan, img)):.1f} km/s"
            else:
                A = data if a == "data" else lib[a]
                B = data if b == "data" else lib[b]
                img = A[key] - B[key]
                lim = 20 if b == "data" else 10
                cmap, lo, hi, unit = CMAP_DIV, -lim, lim, "Δ km/s"
                rms = float(np.sqrt(np.nanmean(np.where(mask, np.nan, img)**2)))
                with np.errstate(invalid="ignore", divide="ignore"):
                    chi = np.nansum(np.where(mask, np.nan, (img / err[key])**2))
                if b == "data":       # Delta chi^2 of this model vs the ground-truth model, on the data
                    with np.errstate(invalid="ignore", divide="ignore"):
                        chi0 = np.nansum(np.where(mask, np.nan, ((lib["ground truth"][key] - data[key]) / err[key])**2))
                    dchi = float(chi - chi0)
                    t2 = (f"rms {rms:.1f},  χ²/pixel {chi0 / (~mask).sum():.2f}" if a == "ground truth" else
                          f"rms {rms:.1f},  Δχ² {dchi:+.0f},  ΔlogL {-dchi / (2 * temper):+.2f}")
                    out = dict(rms_kms=rms, delta_chi2=dchi, dlogL_tempered=-dchi / (2 * temper),
                               chi2_per_pixel=float(chi / (~mask).sum()))
                else:                 # model - ground-truth model: no data noise, but library noise
                    D = np.where(mask, np.nan, img)
                    nb = np.concatenate([(D[:, 1:] - D[:, :-1]).ravel(), (D[1:, :] - D[:-1, :]).ravel()])
                    noise = float(np.sqrt(np.nanmedian(nb**2) / 2))
                    t2 = f"rms {rms:.1f},  pixel-scale library noise {noise:.1f}"
                    out = dict(rms_kms=rms, library_noise_kms=noise, sum_sq_over_err_incl_library_noise=float(chi))
                numbers.setdefault(name, {}).setdefault(label, {})[qname] = out
            img = np.where(mask, np.nan, img)
            im = ax.pcolormesh(e, e, img, cmap=grey_bad(cmap), vmin=lo, vmax=hi, rasterized=True)
            ax.plot([-5, 5], [0, 0], color=ACCENT, lw=0.9, ls=(0, (4, 2)))
            ax.set_aspect("equal")
            ax.set_xlim(-8, 8)
            ax.set_ylim(-8, 8)
            ax.tick_params(labelsize=7)
            ax.set_title(f"{label}: {qname}" + (f"\n{t2}" if t2 else ""), color=INK_2, fontsize=8.5)
            cb = fig.colorbar(im, cax=ax.inset_axes([1.03, 0, 0.045, 1]))
            cb.set_label(unit, fontsize=7.5)
            cb.ax.tick_params(labelsize=7)
            if i == len(rows) - 1:
                ax.set_xlabel("x [kpc] (bar frame)", fontsize=8)
            if j == 0:
                ax.set_ylabel("y [kpc]", fontsize=8)
    fig.suptitle(title, y=1 - 0.12 / (3.45 * len(rows)), fontsize=11, va="top")
    savefig(fig, name)


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        stars = ssl.rotate_z(np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")]),
                             np.radians(TRUTH["phi_b_deg"]))
    stars = stars[np.abs(stars[:, 2]) < Z_CUT]
    data = data_maps(stars)
    err = boot_errors(stars, np.random.default_rng(0))
    mask = data["N"] < N_MIN
    n_pix = int((~mask).sum())
    temper = np.sqrt(2 * 4 * n_pix)
    lib = {"ground truth": library_maps(ARR_DIR / "000.npz")}
    for S, f, lab in ((1.25, 1.0, "S 1.25 (longer)"), (0.8, 1.0, "S 0.8 (shorter)"),
                      (1.0, 1.4, "f 1.4 (stronger)"), (1.0, 0.7, "f 0.7 (weaker)")):
        lib[lab] = library_maps(ARR_DIR / f"maps_S{S:g}_f{f:g}.npz")
    numbers = {"n_pixels": n_pix, "n_stars_abs_z_lt_1": int(len(stars)), "sqrt_2_N_data": float(temper),
               "median_boot_error_kms": {q: float(np.nanmedian(np.where(mask, np.nan, err[(s, k)])))
                                         for q, k, s in QUANTS}}
    sub = (f"0.5 kpc square pixels, |z| < {Z_CUT:g} kpc, bar frame (bar along x, dashed); Ω_b = {TRUTH['Omega_abs']}, "
           f"φ_b = {TRUTH['phi_b_deg']:g}°; grey = < {N_MIN} stars; χ² of these 2D maps, tempered by "
           f"√(2·4·{n_pix}) = {temper:.1f}")
    for name, title, alt in (("test1_nbody_maps_length", "Bar LENGTH", ("S 1.25 (longer)", "S 0.8 (shorter)")),
                             ("test1_nbody_maps_strength", "Bar STRENGTH", ("f 1.4 (stronger)", "f 0.7 (weaker)"))):
        rows = [("data (folded)", "value", "data", None), ("ground truth", "value", "ground truth", None)]
        rows += [(a, "value", a, None) for a in alt]
        rows += [("data error (bootstrap)", "error", None, None)]
        rows += [("ground truth − data", "diff", "ground truth", "data")]
        rows += [(f"{a} − data", "diff", a, "data") for a in alt]
        rows += [(f"{a} − ground truth", "diff", a, "ground truth") for a in alt]
        figure(name, f"{title}: data vs models\n{sub}", rows, data, err, lib, mask, temper, numbers)
    with open(OUT_CSV.parent / "SF_square_maps_numbers.json", "w") as fh:
        json.dump(numbers, fh, indent=2)
    print(f"{len(stars)} stars with |z| < {Z_CUT}; {n_pix} pixels with N >= {N_MIN} (folded); "
          f"median N per pixel {np.median(data['N'][~mask]):.0f}; sqrt(2 N_data) = {temper:.1f}")
    print("median bootstrap error per pixel [km/s]:", {k: round(v, 1) for k, v in numbers["median_boot_error_kms"].items()})
    for name in ("test1_nbody_maps_length", "test1_nbody_maps_strength"):
        print(name)
        for lab, v in numbers[name].items():
            print(f"  {lab:32s} " + "  ".join(f"{k} " + " ".join(f"{n} {x:.2f}" for n, x in d.items()) for k, d in v.items()))
        for lab in [l for l in numbers[name] if l.endswith("− data") and not l.startswith("ground")]:
            print(f"  total over 4 maps, {lab:24s} Delta chi2 {sum(d['delta_chi2'] for d in numbers[name][lab].values()):+7.1f}"
                  f"  tempered Delta logL {sum(d['dlogL_tempered'] for d in numbers[name][lab].values()):+6.2f}")


if __name__ == "__main__":
    main()
