"""OSARG posterior with the χ² scale set by the measured model jitter instead of √(2N).

logL = −χ² / (2 σ), so a χ² difference equal to the jitter scale σ costs 0.5 in logL. σ is measured at the best fit
by chi2_jitter_small_perturbations.py --random (tiny perturbations, same stars and seeds):
  σ_1  std of χ² − χ²(best run): the jitter of one model against a fixed reference;
  σ_2  std of the χ² difference between two perturbed runs (≈ √2 σ_1): two models that both jitter, as any two
       design points do.
The emulator (test1_gp_emulator.Emulator: quadratic mean + Matérn-5/2 GP + white noise) is refitted for each scale
on the 384 evaluated χ² (pipeline covariance, sfratio, K = 10, scheme B); posterior by emcee as in the pipeline.
The old √(2N) scale is included for comparison.

Usage: python posterior_with_jitter_scale.py [--run osarg_zhang2024_sormani22 --case osarg_s22]
Output: results/<run>/posterior_jitter_scale.json, figures/<run>_posterior_vs_chi2_scale.png
"""
import argparse
import json

import numpy as np
import pandas as pd
import emcee
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge

RUN = "osarg_zhang2024_mw_hunter24"
B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"


def posterior(df, scale):
    tr = df[df.kind == "qmc"].copy()
    chi2 = tr[f"chi2_{B}_sfratio"].values
    y = -(chi2 - chi2.min()) / (2 * scale)
    u = ge.to_unit(tr)
    emu = ge.Emulator(u, y)
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    s = emcee.EnsembleSampler(32, 4, lambda q: -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0])
    s.random_state = np.random.RandomState(2).get_state()
    s.run_mcmc(np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((32, 4)), 0, 1), 3000, progress=False)
    post = ge.from_unit(s.get_chain(discard=1000, thin=5, flat=True))
    return post, ge.from_unit(best.x[None])[0], str(emu.gp.kernel_)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=RUN)
    ap.add_argument("--case", default="osarg_h24", help="the chi2_jitter_<case>_random.csv to take sigma from")
    args = ap.parse_args()
    set_style()
    ge.use_run(args.run)
    df = pd.read_csv(RESULTS_DIR / args.run / "evals.csv")
    jit = pd.read_csv(RESULTS_DIR / f"chi2_jitter_{args.case}_random.csv")
    x = jit.loc[jit.param == "random", "C_data"].values
    x0 = float(jit.loc[jit.param == "repeat", "C_data"].iloc[0])
    sig1 = float(np.std(x - x0, ddof=1))
    sig2 = float(np.std((x[:, None] - x[None, :])[np.triu_indices(len(x), 1)], ddof=1))
    nd = int(df[f"ndata_{B}_sfratio"].mode().iloc[0])
    scales = {f"√(2N) = {np.sqrt(2 * nd):.1f} (old)": np.sqrt(2 * nd), f"σ_1 = {sig1:.1f}": sig1,
              f"σ_2 = {sig2:.1f}": sig2}
    res = {}
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2))
    colours = [INK_2, ACCENT, BLUE_RAMP[5]]
    for (lab, sc), col in zip(scales.items(), colours):
        post, mx, kern = posterior(df, sc)
        q = np.percentile(post, [16, 50, 84], axis=0)
        res[lab] = dict(scale=sc, max=dict(zip(ge.KEYS, map(float, mx))), kernel=kern,
                        posterior={k: dict(p16=float(a), p50=float(b), p84=float(c)) for k, a, b, c in zip(ge.KEYS, *q)})
        print(f"[{lab}] " + ", ".join(f"{k} {b:.2f} ({a:.2f}-{c:.2f})" for k, a, b, c in zip(ge.KEYS, *q)), flush=True)
        for j, ax in enumerate(axes):
            ax.hist(post[:, j], bins=40, histtype="step", density=True, color=col, lw=1.8, label=lab if j == 0 else None)
    for j, ax in enumerate(axes):
        ax.set_xlabel(ge.NAMES[j])
        ax.set_yticks([])
    axes[0].legend(frameon=False, fontsize=8.5, loc="upper right")
    fig.suptitle(f"{args.run}: posterior with the χ² scale from √(2N) (old) and from the measured model jitter",
                 color=INK, y=1.02, fontsize=11)
    savefig(fig, f"{args.run}_posterior_vs_chi2_scale")
    json.dump(res, open(RESULTS_DIR / args.run / "posterior_jitter_scale.json", "w"), indent=2)


if __name__ == "__main__":
    main()
