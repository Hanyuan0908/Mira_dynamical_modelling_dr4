"""Gaussian-process emulator of the Test 1 steady-state log-likelihood and the posterior it implies.

Training set: the N_QMC Sobol points of test1_qmc.py; the truth evaluation is kept out and only used
to validate the emulator. Inputs are mapped to the unit cube of test1_qmc.BOX
(|Omega_b|, phi_b, ln S, ln bar_frac). Model: quadratic mean (least squares) + GP on the residuals
with a Matern-5/2 kernel (one length scale per parameter) and a white-noise term.
The posterior is sampled with emcee on the emulator mean, with flat priors in the box and logL
taken at face value (formal likelihood).
Usage: python test1_gp_emulator.py [--run au18|nbody] [--column logL|logL_gauss|chi2_means|...]
       [--rmin R] [--rmax R]   (chi2 columns are converted to log L = -chi2/2; radial cuts use the
       saved per-star log p of the KDE (logL) or local-Gaussian (logL_gauss) likelihood)
"""
import argparse
import json

import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import emcee
import corner
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from common import set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from test1_qmc import OUT_CSV, OUT_DIR, LOGP_DIR, MOCK, BOX, TRUTH
import test1_nbody_qmc as _nq

NAMES = [r"$|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]", r"$\phi_b$ [deg]", r"$S$", r"bar_frac"]
KEYS = ["Omega_abs", "phi_b_deg", "S", "bar_frac"]
RUNS = {
    "au18": dict(csv=OUT_CSV, box=BOX, truth=TRUTH, out=OUT_DIR, mock=MOCK, prefix="test1",
                 perstar=lambda i, col: np.load(LOGP_DIR / f"{i:03d}.npy")),
    "nbody": dict(csv=_nq.OUT_CSV, box=_nq.BOX, truth=_nq.TRUTH, out=_nq.OUT_CSV.parent, mock=_nq.MOCK,
                  prefix="test1_nbody",
                  perstar=lambda i, col: np.load(_nq.ARR_DIR / f"{i:03d}.npz")[
                      "logp_gauss" if col == "logL_gauss" else "logp_kde"]),
}
RUN = RUNS["au18"]
LO = np.array([v[0] for v in BOX.values()])
HI = np.array([v[1] for v in BOX.values()])


def use_run(name):
    """Select the run whose design box, truth and files the module functions use."""
    global RUN, LO, HI
    RUN = RUNS[name]
    LO = np.array([v[0] for v in RUN["box"].values()])
    HI = np.array([v[1] for v in RUN["box"].values()])
    return RUN


def load_table(column="logL", rmin=None, rmax=None):
    """Design table with the chosen likelihood in column 'logL', and a tag for output names."""
    df = pd.read_csv(RUN["csv"]).sort_values("idx").reset_index(drop=True)
    tag = "" if column == "logL" else f"_{column}"
    if rmin is not None or rmax is not None:
        with h5py.File(RUN["mock"], "r") as f:
            R = np.hypot(f["x"][:], f["y"][:])
        sel = (R >= (rmin or 0)) & (R < (rmax or np.inf))
        df["logL"] = [float(np.sum(RUN["perstar"](i, column)[sel])) for i in df["idx"]]
        tag += f"_R{rmin or 0:g}-{rmax or 8:g}"
    elif column.startswith("chi2"):
        df["logL"] = -0.5 * df[column]
    else:
        df["logL"] = df[column]
    return df, tag


def to_unit(df):
    p = np.column_stack([df["Omega_abs"], df["phi_b_deg"], np.log(df["S"]), np.log(df["bar_frac"])])
    return (p - LO) / (HI - LO)


def from_unit(u):
    p = LO + np.asarray(u) * (HI - LO)
    return np.column_stack([p[:, 0], p[:, 1], np.exp(p[:, 2]), np.exp(p[:, 3])])


def quad_features(u):
    u = np.atleast_2d(u)
    cols = [np.ones(len(u))] + [u[:, i] for i in range(4)]
    cols += [u[:, i] * u[:, j] for i in range(4) for j in range(i, 4)]
    return np.column_stack(cols)


class Emulator:
    def __init__(self, u, y):
        self.ymax = y.max()
        y = y - self.ymax
        self.beta, *_ = np.linalg.lstsq(quad_features(u), y, rcond=None)
        r = y - quad_features(u) @ self.beta
        kern = (ConstantKernel(np.var(r), (1e-3 * np.var(r) + 1e-6, 1e3 * np.var(r) + 1)) *
                Matern(length_scale=[0.3] * 4, length_scale_bounds=(0.02, 10), nu=2.5) +
                WhiteKernel(100.0, (1e-2, 1e5)))
        self.gp = GaussianProcessRegressor(kern, normalize_y=False, n_restarts_optimizer=15,
                                           random_state=0).fit(u, r)

    def predict(self, u, return_std=False):
        u = np.atleast_2d(u)
        m, s = self.gp.predict(u, return_std=True)
        mean = quad_features(u) @ self.beta + m + self.ymax
        return (mean, s) if return_std else mean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="au18", choices=list(RUNS))
    ap.add_argument("--column", default="logL")
    ap.add_argument("--rmin", type=float, default=None, help="use only stars with R >= rmin [kpc]")
    ap.add_argument("--rmax", type=float, default=None, help="use only stars with R < rmax [kpc]")
    args = ap.parse_args()
    set_style()
    run = use_run(args.run)
    TRUTH, OUT_DIR, pre = run["truth"], run["out"], run["prefix"]
    df, tag = load_table(args.column, args.rmin, args.rmax)
    tr = df[df.kind == "qmc"]
    truth_row = df[df.kind == "truth"]
    u, y = to_unit(tr), tr["logL"].values
    emu = Emulator(u, y)
    res = dict(run=args.run, column=args.column, n_train=len(tr), kernel=str(emu.gp.kernel_),
               radial_cut=[args.rmin, args.rmax])

    # leave-one-out accuracy
    loo = np.array([Emulator(np.delete(u, i, 0), np.delete(y, i)).predict(u[i:i + 1])[0] for i in range(len(u))])
    res["loo_rms"] = float(np.sqrt(np.mean((loo - y)**2)))
    near = y > y.max() - 2000
    res["loo_rms_within_2000_of_max"] = float(np.sqrt(np.mean((loo[near] - y[near])**2)))
    u_truth = to_unit(pd.DataFrame([TRUTH]))
    pt, st = emu.predict(u_truth, return_std=True)
    if len(truth_row):
        res["truth_logL_evaluated"] = float(truth_row["logL"].iloc[0])
    res["truth_logL_emulated"] = [float(pt[0]), float(st[0])]

    # maximum of the emulator
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    p_best = from_unit(best.x[None])[0]
    res["emulator_max"] = dict(zip(KEYS, map(float, p_best)))
    res["emulator_max_logL"] = float(-best.fun)
    res["best_evaluated"] = {k: float(tr.loc[tr.logL.idxmax(), k]) for k in KEYS + ["logL"]}

    # posterior on the emulator (formal likelihood, flat priors in the box)
    def logpost(q):
        return -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0]
    nw = 32
    p0 = np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((nw, 4)), 0, 1)
    sampler = emcee.EnsembleSampler(nw, 4, logpost)
    sampler.run_mcmc(p0, 3000, progress=False)
    chain = sampler.get_chain(discard=1000, thin=5, flat=True)
    post = from_unit(chain)
    q = np.percentile(post, [16, 50, 84], axis=0)
    res["posterior"] = {k: dict(p16=float(a), p50=float(b), p84=float(c)) for k, a, b, c in zip(KEYS, *q)}
    res["truth"] = TRUTH
    res["truth_offset_in_sigma"] = {k: float((TRUTH[k] - b) / (0.5 * (c - a)))
                                    for k, a, b, c in zip(KEYS, *q)}
    with open(OUT_DIR / f"emulator_results{tag}.json", "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=1))

    # figure 1: 1D slices through the truth and through the maximum, with the design points
    tv = np.array([TRUTH[k] for k in KEYS])
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.6), gridspec_kw=dict(hspace=0.35, wspace=0.3))
    grid = np.linspace(0, 1, 200)
    for j in range(4):
        for row, centre, lab in ((0, u_truth[0], "slice through truth"), (1, best.x, "slice through emulator max")):
            ax = axes[row, j]
            Q = np.repeat(centre[None], len(grid), 0)
            Q[:, j] = grid
            m, s = emu.predict(Q, return_std=True)
            xv = from_unit(Q)[:, j]
            ax.fill_between(xv, m - 2 * s - y.max(), m + 2 * s - y.max(), color=BLUE_RAMP[0], lw=0)
            ax.plot(xv, m - y.max(), color=BLUE_RAMP[4], lw=1.8)
            ax.axvline(tv[j], color=ACCENT, lw=1.2, ls="--")
            ax.axvline(p_best[j], color=INK, lw=0.8, ls=":")
            ax.set_ylim(-3000, 300)
            ax.set_xlabel(NAMES[j])
            if j == 0:
                ax.set_ylabel(r"$\log L - \max_{\rm design}$" + f"\n({lab})")
    fig.suptitle("Test 1 GP emulator: orange dashed = truth, black dotted = emulator maximum; band = ±2σ",
                 y=0.98, fontsize=11)
    savefig(fig, f"{pre}_emulator_slices{tag}")

    # figure 2: posterior corner plot
    fig = corner.corner(post, labels=NAMES, truths=tv, truth_color=ACCENT, color=BLUE_RAMP[5],
                        quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".3f",
                        title_kwargs=dict(fontsize=9), label_kwargs=dict(fontsize=10))
    fig.suptitle("Test 1 posterior from the GP emulator (formal likelihood); orange = truth", y=1.02)
    savefig(fig, f"{pre}_emulator_corner{tag}")

    # figure 3: leave-one-out
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.plot(y - y.max(), loo - y.max(), "o", color=BLUE_RAMP[4], ms=4)
    lim = [min(y.min(), loo.min()) - y.max(), 200]
    ax.plot(lim, lim, color=INK_MUTED, lw=0.8)
    ax.set_xlabel(r"evaluated $\log L$ − max")
    ax.set_ylabel(r"leave-one-out emulator $\log L$ − max")
    ax.set_title(f"LOO rms {res['loo_rms']:.0f} (within 2000 of max: {res['loo_rms_within_2000_of_max']:.0f})",
                 color=INK_2, fontsize=9.5)
    savefig(fig, f"{pre}_emulator_loo{tag}")


if __name__ == "__main__":
    main()
