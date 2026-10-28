"""Data vs model: the m = 0, 2, 4 Fourier coefficients of <v_R>, <v_phi>, sigma_R, sigma_phi per 0.3 kpc ring, i.e. the
numbers the per-star Fourier likelihood compares (test1_nbody_qmc_focused_fourier_stars.py, proposed scheme).

Data: the 50k mock, rotated by the true bar angle (27 deg), with 1-sigma bootstrap errors (diagonal of the joint
covariance). Models (orbit libraries from the focused global search, coefficients from fourier_stars_NNN.npz):
  #0    the truth (also the best of the 257 evaluated points),
  #235  and #15, the two best others, on the S - bar_frac ridge (longer/weaker and shorter/stronger),
  #17   a rejected model: S, bar_frac, phi_b close to the truth but |Omega_b| 33.4 instead of 26.0.
Grey line: the same fit to all 1.49 M N-body star particles (what a perfect model would give).
Dispersions are shown in km/s: m = 0 as sigma_0 = sqrt(c_0) and the other terms as c_m / (2 sigma_0), the first-order
Fourier coefficients of sigma(phi) (the fit itself is of sigma^2(phi)). The data are plotted at 27 deg; each model's own
trial angle differs by <= 1.2 deg, and the per-panel chi^2 / N (diagonal errors, R < 6 kpc) uses the data at the
model's own angle, as the likelihood does. Shaded: R > 6 kpc, outside the proposed scheme.
Figures: figures/test1_nbody_qmc_focused_fourier_coefficients_data_vs_models_absz_{0to1,1to2}.png
Numbers: results/test1_nbody_qmc_focused/fourier_coefficients_data_vs_models.json
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import ROOT, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import test1_nbody_qmc_focused_fourier_stars as fs
from test1_nbody_qmc_focused_run import OUT_DIR, DESIGN_CSV
from check_fourier_star_fit_vs_true_maps import load, POOL

MODELS = [(0, "truth model (#0)", ACCENT, "-", 1.8), (235, "#235", BLUE_RAMP[4], "--", 1.3),
          (15, "#15", BLUE_RAMP[6], "-.", 1.3), (17, "#17 (rejected)", INK, ":", 1.4)]
FIELDS = [("<v_R>", 0, "mean"), ("<v_phi>", 1, "mean"), ("sigma_R", 3, "disp"), ("sigma_phi", 4, "disp")]
TERMS = ["m = 0", "cos 2φ", "sin 2φ", "cos 4φ", "sin 4φ"]
LAYERS = [(0, "absz_0to1", "|z| < 1 kpc"), (1, "absz_1to2", "1 < |z| < 2 kpc")]
RC = 0.5 * (fs.RING_EDGES[1:] + fs.RING_EDGES[:-1])
NK = fs.NK


def as_kms(c, err=None, kind="mean"):
    """Coefficients (NR, 5) of a mean field as they are; of a variance field converted to sigma terms in km/s."""
    if kind == "mean":
        return c, err
    s0 = np.sqrt(np.maximum(c[:, 0], 0))
    out = np.column_stack([s0, c[:, 1:] / (2 * s0[:, None])])
    if err is None:
        return out, None
    return out, np.column_stack([err[:, 0] / (2 * s0), err[:, 1:] / (2 * s0[:, None])])


def main():
    set_style()
    design = pd.read_csv(DESIGN_CSV).set_index("idx")
    evals = pd.read_csv(ROOT / "results" / "test1_nbody_qmc_focused" / "evals.csv").set_index("idx")
    col = "logLs_fourier_stars_R0to6_vR_m24_vphi_sigR_sigphi_m024"
    npz = {i: np.load(OUT_DIR / f"fourier_stars_{i:03d}.npz") for i, *_ in MODELS}
    shape = (fs.NR, fs.NL, 6, NK)
    d0 = npz[0]["coef_data"].reshape(shape)
    e0 = np.sqrt(np.diagonal(npz[0]["cov_data_joint"], axis1=2, axis2=3)).reshape(shape)
    valid0 = npz[0]["valid"]
    ok = npz[0]["n_stars"].reshape(fs.NB) >= fs.N_MIN_BLOCK
    bp, Xp, vp, _ = fs.prepare(load(POOL, sel_R=8.0))
    c_nbody = fs.fit(bp, Xp, vp, np.ones(len(bp)), ok).reshape(shape)

    res = dict(models={str(i): dict(label=lab, **{k: float(design.loc[i, k]) for k in ("Omega_abs", "phi_b_deg", "S",
                                                                                         "bar_frac")},
                                    dlogL_vs_truth=round(float(evals.loc[i, col] - evals.loc[0, col]), 3),
                                    chi2_per_N_scheme=round(float(evals.loc[i, col.replace("logLs_", "chi2_")]
                                                                  / evals.loc[i, col.replace("logLs_", "ndata_")]), 3))
                       for i, lab, *_ in MODELS}, panels={})
    # how much the data coefficients move between 27 deg and each model's angle, in units of the data error
    shift = {}
    for i, *_ in MODELS:
        dd = (npz[i]["coef_data"].reshape(shape) - d0) / e0
        shift[str(i)] = round(float(np.nanmax(np.abs(dd[:fs.N_RING_DEFAULT, :2][..., [0, 1, 3, 4], :]))), 3)
    res["max_data_shift_vs_27deg_in_errors"] = shift

    for layer, tag, label in LAYERS:
        fig, axes = plt.subplots(len(FIELDS), len(TERMS), figsize=(17.5, 12.5),
                                 gridspec_kw=dict(hspace=0.42, wspace=0.32))
        use = valid0[:, layer]
        r6 = use & (np.arange(fs.NR) < fs.N_RING_DEFAULT)
        res["panels"][tag] = {}
        for r, (name, f, kind) in enumerate(FIELDS):
            dat, err = as_kms(d0[:, layer, f], e0[:, layer, f], kind)
            nb, _ = as_kms(c_nbody[:, layer, f], kind=kind)
            for k, term in enumerate(TERMS):
                ax = axes[r, k]
                ax.axvspan(6.0, 8.0, color=INK_MUTED, alpha=0.12, lw=0)
                ax.axhline(0, color=INK_MUTED, lw=0.6)
                ax.plot(RC[use], nb[use, k], color="#b9b9b4", lw=2.6, solid_capstyle="round",
                        label="all 1.49 M N-body particles" if (r, k) == (0, 0) else None, zorder=1)
                notes = []
                for i, lab, colr, ls, lw in MODELS:
                    lib, _ = as_kms(npz[i]["coef_lib"].reshape(shape)[:, layer, f], kind=kind)
                    ax.plot(RC[use], lib[use, k], color=colr, ls=ls, lw=lw, zorder=3,
                            label=f"{lab}: Ω {design.loc[i, 'Omega_abs']:.1f}, φ {design.loc[i, 'phi_b_deg']:.1f}°, "
                                  f"S {design.loc[i, 'S']:.2f}, f {design.loc[i, 'bar_frac']:.2f}"
                            if (r, k) == (0, 0) else None)
                    # chi2 / N of this panel against the data at the model's own angle (diagonal errors, R < 6)
                    di, ei = as_kms(npz[i]["coef_data"].reshape(shape)[:, layer, f],
                                    np.sqrt(np.diagonal(npz[i]["cov_data_joint"], axis1=2, axis2=3)).reshape(shape)[:, layer, f],
                                    kind)
                    pull = ((di[:, k] - lib[:, k]) / ei[:, k])[r6]
                    chi = float(np.mean(pull**2))
                    res["panels"][tag].setdefault(f"{name} {term}", {})[str(i)] = round(chi, 3)
                    notes.append(f"{chi:.2f}")
                ax.errorbar(RC[use], dat[use, k], yerr=err[use, k], fmt="o", ms=3.2, color=INK, ecolor=INK_2,
                            elinewidth=0.9, capsize=0, zorder=4, label="data (50k mock)" if (r, k) == (0, 0) else None)
                ax.set_xlim(0, 8)
                ax.tick_params(labelsize=7.5)
                if r == 0:
                    ax.set_title(term, fontsize=10)
                if k == 0:
                    unit = "km/s" if kind == "mean" else "σ [km/s]"
                    ax.set_ylabel(f"{name}\n{unit}", fontsize=9.5)
                if r == len(FIELDS) - 1:
                    ax.set_xlabel("R [kpc]", fontsize=9)
                excluded = name == "<v_R>" and k == 0
                ax.text(0.02, 0.97, ("not in the likelihood; " if excluded else "")
                        + "χ²/N  #0 " + " / ".join(notes[:1]) + ",  #235 " + notes[1] + ",  #15 " + notes[2]
                        + ",  #17 " + notes[3], transform=ax.transAxes, fontsize=6.3, va="top", color=INK_2,
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
        fig.legend(loc="upper center", ncol=3, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.965))
        fig.suptitle(f"Fourier coefficients per 0.3 kpc ring, data vs models, bar frame, {label}. χ²/N per panel: "
                     f"diagonal errors, R < 6 kpc (shaded: outside the proposed scheme). σ terms: m = 0 is σ₀, others "
                     f"c_m / 2σ₀", fontsize=10.5, y=0.995)
        savefig(fig, f"test1_nbody_qmc_focused_fourier_coefficients_data_vs_models_{tag}")
    out = ROOT / "results" / "test1_nbody_qmc_focused" / "fourier_coefficients_data_vs_models.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "panels"}, indent=1))


if __name__ == "__main__":
    main()
