"""How much hotter is the ground-truth orbit library than the data? Ring-averaged dispersions (the m = 0 Fourier term
used in test1_fourier_compression.py: error-weighted mean of the 15-deg sector values with N >= N_MIN) of the N-body
mock and of the libraries at the true Omega_b, phi_b (truth, S 1.25, S 0.8, f 1.4, f 0.7), per R ring and |z| layer.
Data dispersions are unbiased (N / (N - 1)); errors from the bootstrap.
Figure: figures/test1_nbody_ring_avg_dispersion_residuals_vs_R.png (model minus data, km/s, with the data error band).
Numbers: results/test1_nbody/ring_avg_dispersion_residuals.json
"""
import json

import numpy as np
import h5py
import matplotlib.pyplot as plt

from common import set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import steady_state_likelihood as ssl
from test1_nbody_moments import boot_errors, N_MIN, SEED
from test1_nbody_qmc import MOCK, TRUTH, OUT_CSV
from test1_fourier_compression import lib

NR, NP, NZ = len(ssl.CELL_R) - 1, len(ssl.CELL_PHI) - 1, len(ssl.CELL_Z) - 1
RC = 0.5 * (ssl.CELL_R[1:] + ssl.CELL_R[:-1])
LAYERS = ["|z| 0–0.5", "|z| 0.5–1", "|z| 1–2", "|z| > 2"]
N_PLOT = 3                     # the |z| > 2 layer has no sector with >= N_MIN stars
QUANTS = [("σ_R", 0), ("σ_φ", 1), ("σ_z", 2)]
MODELS = [("ground truth", "pilot_00.npz", INK, "-"), ("S 1.25 (longer)", "pilot_09.npz", ACCENT, "-"),
          ("S 0.8 (shorter)", "pilot_08.npz", ACCENT, "--"), ("f 1.4 (stronger)", "pilot_07.npz", BLUE_RAMP[5], "-"),
          ("f 0.7 (weaker)", "pilot_06.npz", BLUE_RAMP[5], "--")]


def ring_average(val, err, good):
    """Error-weighted mean over the usable sectors of each (ring, layer): value and its error, (NR, NZ)."""
    w = np.where(good, 1.0 / err**2, 0.0)
    W = w.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nansum(w * np.nan_to_num(val), 1) / W, 1.0 / np.sqrt(W)


def main():
    set_style()
    with h5py.File(MOCK, "r") as f:
        stars = np.column_stack([f[k][:].astype(np.float64) for k in ("x", "y", "z", "vx", "vy", "vz")])
    N, m1, s, em, es = boot_errors(ssl.rotate_z(stars, np.radians(TRUTH["phi_b_deg"])), np.random.default_rng(SEED))
    libs = {lab: lib(fn) for lab, fn, _, _ in MODELS}
    good = ((N >= N_MIN) & np.all(np.isfinite(libs["ground truth"][0]), 1) & np.all(em > 0, 1)
            & np.all(es > 0, 1)).reshape(NR, NP, NZ)
    out = {}
    fig, axes = plt.subplots(N_PLOT, 3, figsize=(15, 10.5), sharex=True, gridspec_kw=dict(hspace=0.36, wspace=0.25))
    for j, (qn, q) in enumerate(QUANTS):
        e = es[:, q].reshape(NR, NP, NZ)
        d_avg, d_err = ring_average(s[:, q].reshape(NR, NP, NZ), e, good)
        o = out.setdefault(qn, {"data": d_avg.tolist(), "data_err": d_err.tolist()})
        for lab, fn, c, ls in MODELS:
            l_avg, _ = ring_average(libs[lab][1][:, q].reshape(NR, NP, NZ), e, good)
            o[lab] = (l_avg - d_avg).tolist()
            for k in range(N_PLOT):
                ax = axes[k, j]
                ax.plot(RC, l_avg[:, k] - d_avg[:, k], ls, color=c, lw=1.8 if lab == "ground truth" else 1.2,
                        marker="o" if lab == "ground truth" else None, ms=3.5, label=lab)
        for k in range(N_PLOT):
            ax = axes[k, j]
            ax.fill_between(RC, -d_err[:, k], d_err[:, k], color=INK_MUTED, alpha=0.25, lw=0, label="±1σ data error")
            ax.axhline(0, color=INK_MUTED, lw=0.8)
            r = np.array(o["ground truth"])[:, k]
            ok = np.isfinite(r)
            w = 1 / d_err[ok, k]**2
            ax.set_title(f"{qn} ring average, {LAYERS[k]}: model − data\nground truth hotter by "
                         f"{np.sum(w * r[ok]) / w.sum():+.2f} ± {1 / np.sqrt(w.sum()):.2f} km/s (weighted mean)",
                         color=INK_2, fontsize=9)
            ax.set_ylabel("km/s")
            ax.set_ylim(-9, 12)
            if k == N_PLOT - 1:
                ax.set_xlabel("R [kpc]")
    axes[0, 0].legend(fontsize=7.5, loc="lower right")
    fig.suptitle("Ring-averaged velocity dispersions (m = 0 term): orbit libraries at the true Ω_b, φ_b minus the N-body "
                 "data\n(error-weighted mean of the 15° sectors with ≥ 20 stars; data dispersions unbiased, N/(N−1); "
                 "band = ±1σ bootstrap error of the data ring average)", y=0.995, fontsize=11)
    savefig(fig, "test1_nbody_ring_avg_dispersion_residuals_vs_R")

    # ---- numbers: truth library minus data, pooled over rings, per field and layer ----
    summary = {}
    for qn, _ in QUANTS:
        o = out[qn]
        d, de = np.array(o["data"]), np.array(o["data_err"])
        for lab, *_ in MODELS:
            r = np.array(o[lab])
            ok = np.isfinite(r) & np.isfinite(de)
            w = 1 / de[ok]**2
            summary.setdefault(qn, {})[lab] = dict(
                weighted_mean_kms=float(np.sum(w * r[ok]) / w.sum()), err_kms=float(1 / np.sqrt(w.sum())),
                median_kms=float(np.median(r[ok])), median_rel_pct=float(100 * np.median(r[ok] / d[ok])),
                frac_positive=float(np.mean(r[ok] > 0)), median_ring_err_kms=float(np.median(de[ok])),
                n_rings=int(ok.sum()))
    with open(OUT_CSV.parent / "ring_avg_dispersion_residuals.json", "w") as fh:
        json.dump(dict(per_ring=out, summary=summary, R=RC.tolist(), layers=LAYERS), fh, indent=1)
    for qn, _ in QUANTS:
        print(qn)
        for lab, x in summary[qn].items():
            print(f"  {lab:18s} model − data: weighted mean {x['weighted_mean_kms']:+.2f} ± {x['err_kms']:.2f} km/s, "
                  f"median {x['median_kms']:+.2f} km/s ({x['median_rel_pct']:+.2f}%), "
                  f"{100 * x['frac_positive']:.0f}% of {x['n_rings']} rings hotter; median ring error {x['median_ring_err_kms']:.2f}")


if __name__ == "__main__":
    main()
