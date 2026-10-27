"""Analyse the N-body Test 1 pilot: Delta log L (offset - truth) for each likelihood, compared with the
noise floor (truth evaluated at a different integrator accuracy), and its radial decomposition.

Likelihoods: per-star KDE; per-star local Gaussian N(v | mu, C); two robust variants computed from the
saved mu, C (Student-t with nu = 5 and the same scale matrix; Gaussian mixed with a 2% component of
9x the covariance); moment chi^2 (means, all R / R >= 1 kpc; means + dispersions, R >= 1 kpc) as -chi2/2.
"""
import json

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from scipy.special import gammaln

from common import set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from test1_nbody_pilot import OUT_DIR, MOCK

NU, EPS, INFL = 5.0, 0.02, 9.0
R_BINS = ((0, 1), (1, 2), (2, 4), (4, 6), (6, 8))


def per_star_variants(v, mu, C):
    """Student-t and outlier-mixture log densities from local mean and covariance."""
    r = v - mu
    Cinv = np.linalg.inv(C)
    maha = np.einsum("ni,nij,nj->n", r, Cinv, r)
    _, logdet = np.linalg.slogdet(C)
    d = 3
    # Student-t with covariance C (scale matrix C (nu-2)/nu)
    scale_logdet = logdet + d * np.log((NU - 2) / NU)
    t = (gammaln((NU + d) / 2) - gammaln(NU / 2) - 0.5 * d * np.log(NU * np.pi) - 0.5 * scale_logdet
         - 0.5 * (NU + d) * np.log1p(maha * NU / (NU - 2) / NU))
    g = -0.5 * maha - 0.5 * (d * np.log(2 * np.pi) + logdet)
    gb = -0.5 * maha / INFL - 0.5 * (d * np.log(2 * np.pi) + logdet + d * np.log(INFL))
    mix = np.logaddexp(np.log(1 - EPS) + g, np.log(EPS) + gb)
    return t, mix


def main():
    set_style()
    df = pd.read_csv(OUT_DIR / "pilot.csv")
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    R = np.hypot(stars[:, 0], stars[:, 1])
    per = {}
    for i, row in df.iterrows():
        a = np.load(OUT_DIR / f"pilot_{int(row.i):02d}.npz")
        # star velocities in the bar frame of this evaluation (rotation by the trial bar angle)
        ang = np.radians(row.phi_b_deg)
        c, s = np.cos(ang), np.sin(ang)
        v = np.column_stack([c * stars[:, 3] - s * stars[:, 4], s * stars[:, 3] + c * stars[:, 4], stars[:, 5]])
        t, mix = per_star_variants(v, a["mu"].astype(np.float64), a["C"].astype(np.float64))
        per[int(row.i)] = dict(kde=a["logp_kde"].astype(np.float64), gauss=a["logp_gauss"].astype(np.float64),
                               student_t=t, gauss_mix=mix)
    i_truth = int(df.loc[df.label == "truth", "i"].iloc[0])
    names = {"kde": "per-star KDE", "gauss": "per-star Gaussian", "student_t": "per-star Student-t (ν=5)",
             "gauss_mix": "per-star Gaussian + 2% outliers", "chi2_means": "moments: means (all R)",
             "chi2_means_R1": "moments: means (R ≥ 1)", "chi2_meansdisp_R1": "moments: means+disp (R ≥ 1)"}
    table, radial = {}, {}
    for _, row in df.iterrows():
        i = int(row.i)
        if i == i_truth:
            continue
        d = {k: float(np.sum(per[i][k] - per[i_truth][k])) for k in ("kde", "gauss", "student_t", "gauss_mix")}
        for k in ("chi2_means", "chi2_means_R1", "chi2_meansdisp_R1"):
            d[k] = float(-0.5 * (row[k] - df.loc[df.i == i_truth, k].iloc[0]))
        table[row.label] = d
        radial[row.label] = {k: {f"{a}-{b}": float(np.sum((per[i][k] - per[i_truth][k])[(R >= a) & (R < b)]))
                                 for a, b in R_BINS} for k in ("kde", "gauss", "gauss_mix")}
    noise_label = [l for l in table if "noise" in l][0]
    print(f"Delta log L = offset - truth (noise floor: '{noise_label}')")
    print(f"{'offset':24s}" + "".join(f"{names[k][:22]:>24s}" for k in names))
    for lab, d in table.items():
        print(f"{lab:24s}" + "".join(f"{d[k]:24.1f}" for k in names))
    print("\nradial decomposition (per-star likelihoods):")
    for lab, rr in radial.items():
        for k, v in rr.items():
            print(f"  {lab:24s} {k:9s} " + "  ".join(f"R {b}: {x:+8.1f}" for b, x in v.items()))
    with open(OUT_DIR / "pilot_summary.json", "w") as f:
        json.dump(dict(table=table, radial=radial), f, indent=2)

    # figure: Delta log L in units of the noise floor, per likelihood and offset
    offs = [l for l in table if l != noise_label]
    fig, ax = plt.subplots(figsize=(11, 4.6))
    cols = [BLUE_RAMP[5], BLUE_RAMP[3], BLUE_RAMP[1], "#1baf7a", ACCENT, "#e87ba4", INK_MUTED]
    for j, (k, col) in enumerate(zip(names, cols)):
        nz = max(abs(table[noise_label][k]), 1.0)
        ax.plot(np.arange(len(offs)) + (j - 3) * 0.1, [table[o][k] / nz for o in offs], "o", color=col, ms=6,
                label=f"{names[k]} (noise {table[noise_label][k]:+.1f})")
    ax.axhline(0, color=INK, lw=0.8)
    ax.axhspan(-1, 1, color="#f0efec", lw=0, zorder=0)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(range(len(offs)), offs)
    ax.set_ylabel(r"$\Delta\log L$ (offset − truth) / |noise|")
    ax.set_title("N-body Test 1 pilot: every likelihood should be negative (truth preferred); "
                 "grey band = noise floor", color=INK_2, fontsize=10)
    ax.legend(fontsize=7.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    savefig(fig, "test1_nbody_pilot")


if __name__ == "__main__":
    main()
