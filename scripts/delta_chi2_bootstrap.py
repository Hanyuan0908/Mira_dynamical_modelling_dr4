"""Confidence regions from the scatter of chi2 DIFFERENCES (Wang, Thomas, Lipka et al. 2026; Thomas & Lipka 2022),
instead of the tempering logL = -chi2 / (2 sqrt(2 N)).

Wang et al. draw mock data from the best model plus noise at the observed errors, refit, and find that the AIC of one
model scatters by ~sqrt(2 N) but the DIFFERENCE between two models on the same noise realisation scatters much less
(their sigma_dAIC ~ 11 against ~42), because the noise is common to both. For two predictions m1, m2 the difference
has std 2 Delta with Delta^2 = (m1 - m2)^T C^-1 (m1 - m2) (PROGRESS 2026-11-01).

Here, using only the saved library coefficients of every design point (no integration):
  1. whiten every design point's prediction with the data covariance (scheme B terms, Hartlap), reproducing the
     pipeline's chi2 (checked against evals.csv);
  2. parametric bootstrap: N_MC mock data sets = best prediction + Gaussian noise, chi2 of every design point
     against each mock. sigma(dchi2) between the best and each other point is measured and compared with 2 Delta.
     The threshold for a confidence region is the 68% / 95% quantile of dchi2(best model vs grid minimum) over the
     mocks: the likelihood-ratio threshold for this grid, look-elsewhere included;
  3. optionally adding LIBRARY noise. Our predictions are Monte Carlo estimates: libraries at |Omega_b| +-0.25
     km/s/kpc already differ by D ~ 2.3-3.1 tempered (scan_omega_at_fixed_bar.py), i.e. dchi2 ~ 170-230 between
     essentially identical potentials. Modelled as isotropic noise in whitened space with E|eta_i - eta_j|^2 =
     D_floor (an assumption: real library noise is neither isotropic nor independent between points), drawn afresh
     per mock for every design point.
The confidence region on the real data is the set of design points with observed dchi2 below the threshold.

Usage: python delta_chi2_bootstrap.py
Output: results/delta_chi2_bootstrap.csv, printed summary
"""
import json

import numpy as np
import pandas as pd

from common import ROOT, RESULTS_DIR, BAR_ANGLE_DEG
import test1_nbody_qmc_focused_fourier_stars as fs
from test2b_qmc_distance_error_run import load_design
import diagnose_library_sensitivity as dl

N_MC = 500
SEED = 2126
B = dl.B
FLOOR_TEMPERED = {"OSARG data, Hunter+24": 2.4, "t = 7 mock, own family": 2.7}   # D at +-0.25 km/s/kpc (scan)
TRUTH = {"t = 7 mock, own family": 25.20}


def whitened(run, sdir, var, mock, mask):
    meta = json.load(open(RESULTS_DIR / run / "run_meta.json"))
    design, _ = load_design(meta["galaxy"], meta.get("design_csv") and ROOT / meta["design_csv"])
    design = design.set_index("idx")
    ev = pd.read_csv(RESULTS_DIR / run / "evals.csv").set_index("idx")
    ev = ev[ev.kind == "qmc"]
    (cd0, cov0, _, ok), _ = dl.data_fit(mock, mask)
    nr, spec = fs.SCHEMES[B]
    idx = [f * fs.NK + k for f, terms in spec.items() for k in terms]
    libs = {}
    for i in ev.index:
        cl = np.load(f"{sdir}/{var}_{i:03d}.npz")["coef_lib"].reshape(fs.NB, -1)
        T = fs.rotation(BAR_ANGLE_DEG - float(design.loc[i, "phi_b_deg"]))
        libs[i] = cl @ T.T
    blocks = [r * fs.NL + z for r in range(nr) for z in range(fs.NL)
              if ok[r * fs.NL + z] and all(np.all(np.isfinite(libs[i][r * fs.NL + z][idx])) for i in libs)]
    hart = np.sqrt((fs.N_BOOT - len(idx) - 2) / (fs.N_BOOT - 1))
    W = {b: hart * np.linalg.inv(np.linalg.cholesky(cov0[b][np.ix_(idx, idx)])) for b in blocks}
    vec = lambda c: np.concatenate([W[b] @ c[b][idx] for b in blocks])
    y_d = vec(cd0.reshape(fs.NB, -1))
    Y = np.array([vec(libs[i]) for i in ev.index])
    chi2 = ((Y - y_d) ** 2).sum(1)
    saved = ev[f"chi2_{B}_{var}"].values
    return design, ev, Y, y_d, chi2, saved, len(blocks) * len(idx)


def main():
    rows = []
    rng = np.random.default_rng(SEED)
    for case in (dl.CASES[0], dl.CASES[2]):
        lab, run, sdir, var, mock, mask = case[:6]
        design, ev, Y, y_d, chi2, saved, ndim = whitened(run, sdir, var, mock, mask)
        n_pts = len(Y)
        om = design.loc[ev.index, "Omega_abs"].values
        ib = int(np.argmin(chi2))
        print(f"== {lab}: {n_pts} design points, {ndim} whitened dimensions; chi2 reconstruction vs evals.csv: "
              f"max |diff| {np.max(np.abs(chi2 - saved)):.2f} (median chi2 {np.median(chi2):.0f}); best idx "
              f"{ev.index[ib]} (Omega {om[ib]:.2f}), chi2 {chi2[ib]:.1f}")
        delta = np.sqrt(((Y - Y[ib]) ** 2).sum(1))                     # model separation from the best
        floor = FLOOR_TEMPERED[lab] * 2 * np.sqrt(2 * ndim)            # raw dchi2 between near-identical libraries
        for with_lib in (False, True):
            sig_eta = np.sqrt(floor / (2 * ndim)) if with_lib else 0.0
            d_best, d_nb = [], []
            near = np.argsort(delta)[1:21]
            for _ in range(N_MC):
                y_m = Y[ib] + rng.standard_normal(ndim)
                Ym = Y + sig_eta * rng.standard_normal(Y.shape) if with_lib else Y
                c = ((Ym - y_m) ** 2).sum(1)
                d_best.append(c[ib] - c.min())
                d_nb.append(c[near] - c[ib])
            d_best, d_nb = np.array(d_best), np.array(d_nb)
            q68, q95 = np.percentile(d_best, [68, 95])
            obs = chi2 - chi2.min()
            reg68, reg95 = obs <= q68, obs <= q95
            tag = "data + library noise" if with_lib else "data noise only (Wang et al.)"
            ratio = np.median(d_nb.std(0) / (2 * delta[near]))
            print(f"   [{tag}] sigma(dchi2) of the 20 nearest points / (2 Delta): median {ratio:.2f} "
                  f"(2 Delta {np.median(2 * delta[near]):.0f}); threshold dchi2: 68% {q68:.1f}, 95% {q95:.1f} "
                  f"(sqrt(2N) = {np.sqrt(2 * ndim):.1f})")
            print(f"      confidence region on the data: 68% {reg68.sum()} points, Omega {om[reg68].min():.2f}-"
                  f"{om[reg68].max():.2f}; 95% {reg95.sum()} points, Omega {om[reg95].min():.2f}-{om[reg95].max():.2f}"
                  + (f"; truth {TRUTH[lab]} inside 95%: {om[reg95].min() <= TRUTH[lab] <= om[reg95].max()}"
                     if lab in TRUTH else ""))
            rows.append(dict(case=lab, noise=tag, q68=q68, q95=q95, n68=int(reg68.sum()), n95=int(reg95.sum()),
                             om68_lo=om[reg68].min(), om68_hi=om[reg68].max(), om95_lo=om[reg95].min(),
                             om95_hi=om[reg95].max(), sigma_over_2delta=ratio))
        # what the current tempering implies, for comparison: points within dlogLs 0.5 / 2
        t = (chi2 - chi2.min()) / (2 * np.sqrt(2 * ndim))
        print(f"   [current tempering] dlogLs <= 0.5: {np.sum(t <= 0.5)} points, Omega {om[t <= 0.5].min():.2f}-"
              f"{om[t <= 0.5].max():.2f}; <= 2: {np.sum(t <= 2)} points, Omega {om[t <= 2].min():.2f}-{om[t <= 2].max():.2f}")
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "delta_chi2_bootstrap.csv", index=False)


if __name__ == "__main__":
    main()
