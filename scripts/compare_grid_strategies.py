"""Which Fourier grid gives the best pattern-speed constraint, judged at each grid's OWN temperature?

For every grid of run_grid_tests.sh (and the default):
  - T = sigma_2 of that grid (chi2_jitter_per_grid.py): the chi2 jitter between two near-identical models;
  - the posterior at T (emulator refitted on -dchi2 / (2 T), sampled as the pipeline), and the pipeline's own
    sqrt(2N)-tempered posterior for reference (emulator_results_<column>.json);
  - an emulator-free check: the best evaluated chi2 per 2.5 km/s/kpc bin of |Omega_b|, above the overall best, in units of
    that grid's sigma_2. This is how many "jitter units" the data separate pattern speeds by on that grid.
With --truth (a mock), bias and whether the truth lies inside 16-84% are reported.

Usage: python compare_grid_strategies.py --run osarg_zhang2024_sormani22 --case osarg_s22
       python compare_grid_strategies.py --run mw_hunter24_mock_from_osarg_ics --case h24mock --truth 32.5 25 1 1
Output: results/<run>/grid_strategies.csv, figures/<run>_pattern_speed_posterior_by_grid.png
"""
import argparse
import json

import numpy as np
import pandas as pd
import emcee
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
GRIDS = [("default", "", "0.3 kpc rings, |z| 0-1, 1-2, 2-3"), ("_r0.5_z2", "_r0.5_z2", "0.5 kpc, |z| 0-1, 1-2"),
         ("_r0.5_z2_dz2", "_r0.5_z2_dz2", "0.5 kpc, |z| < 2"), ("_r1_z2_dz2", "_r1_z2_dz2", "1 kpc, |z| < 2"),
         ("_r0.5_z3_dz3", "_r0.5_z3_dz3", "0.5 kpc, |z| < 3"), ("_r0.3_z2_dz2", "_r0.3_z2_dz2", "0.3 kpc, |z| < 2")]


def posterior(u, y):
    emu = ge.Emulator(u, y)
    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -emu.predict(q)[0], x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    s = emcee.EnsembleSampler(32, 4, lambda q: -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0])
    s.random_state = np.random.RandomState(2).get_state()
    s.run_mcmc(np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((32, 4)), 0, 1), 3000, progress=False)
    return ge.from_unit(s.get_chain(discard=1000, thin=5, flat=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--truth", type=float, nargs=4, default=None, metavar=("OMEGA", "PHI", "S", "F"))
    args = ap.parse_args()
    set_style()
    ge.use_run(args.run)
    ev = pd.read_csv(RESULTS_DIR / args.run / "evals.csv")
    tr = ev[ev.kind == "qmc"].reset_index(drop=True)
    u = ge.to_unit(tr)
    om = tr.Omega_abs.values
    jit = pd.read_csv(RESULTS_DIR / f"chi2_jitter_per_grid_{args.case}.csv")
    rows, chains = [], {}
    edges = np.arange(15, 45.01, 2.5)
    for label, tag, desc in GRIDS:
        col = f"chi2_{B}_sfratio{tag}"
        if col not in tr:
            print(f"missing column {col}")
            continue
        chi2 = tr[col].values
        nd = int(tr[f"ndata_{B}_sfratio{tag}"].mode().iloc[0])
        j = jit[jit.grid == label]
        x0, x = j.chi2.iloc[0], j.chi2.values[1:]
        s2 = float(np.std((x[:, None] - x[None, :])[np.triu_indices(len(x), 1)], ddof=1))
        post = posterior(u, -(chi2 - chi2.min()) / (2 * s2))
        chains[label] = post
        q = np.percentile(post, [16, 50, 84], axis=0)
        pj = RESULTS_DIR / args.run / f"emulator_results_logLs_{B}_sfratio{tag}.json"
        pipe = json.load(open(pj)) if pj.exists() else None
        prof = [(chi2[(om >= a) & (om < b)].min() - chi2.min()) / s2 for a, b in zip(edges[:-1], edges[1:])]
        r = dict(grid=desc, tag=label, N_data=nd, best_chi2_per_N=chi2.min() / nd, sigma_2=s2, sqrt_2N=np.sqrt(2 * nd),
                 chi2_span_over_sigma2=(np.percentile(chi2, 90) - chi2.min()) / s2,
                 Om_p16=q[0, 0], Om_p50=q[1, 0], Om_p84=q[2, 0], Om_width_pct=100 * (q[2, 0] - q[0, 0]) / q[1, 0],
                 phi_p50=q[1, 1], S_p50=q[1, 2], S_p16=q[0, 2], S_p84=q[2, 2], f_p50=q[1, 3], f_p16=q[0, 3], f_p84=q[2, 3],
                 loo_rms_tempered=pipe["loo_rms"] if pipe else np.nan)
        if pipe:
            o = pipe["posterior"]["Omega_abs"]
            r.update(pipe_Om=f"{o['p50']:.2f} ({o['p16']:.2f}-{o['p84']:.2f})", pipe_width_pct=100 * (o["p84"] - o["p16"]) / o["p50"])
        if args.truth:
            r.update(Om_bias_pct=100 * (q[1, 0] / args.truth[0] - 1), Om_inside=bool(q[0, 0] <= args.truth[0] <= q[2, 0]))
        rows.append(r)
        print(f"{desc:30s} N {nd:4d}  chi2/N {chi2.min() / nd:.3f}  sigma_2 {s2:5.1f} (sqrt(2N) {np.sqrt(2 * nd):5.1f})  "
              f"|Omega_b| at T=sigma_2: {q[1, 0]:.2f} ({q[0, 0]:.2f}-{q[2, 0]:.2f}) width {r['Om_width_pct']:.1f}%"
              + (f" bias {r['Om_bias_pct']:+.1f}% inside {r['Om_inside']}" if args.truth else "")
              + (f"; sqrt(2N): {r['pipe_Om']} {r['pipe_width_pct']:.1f}%" if pipe else "")
              + f"; S {q[1, 2]:.2f} ({q[0, 2]:.2f}-{q[2, 2]:.2f}), f {q[1, 3]:.2f} ({q[0, 3]:.2f}-{q[2, 3]:.2f}); LOO {r['loo_rms_tempered']:.2f}",
              flush=True)
        print(f"{'':30s} best chi2 per 2.5 bin above the best, in sigma_2, Omega 15 -> 45: "
              + " ".join(f"{v:4.1f}" for v in prof), flush=True)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / args.run / "grid_strategies.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.5, 5))
    cols = [INK, BLUE_RAMP[2], BLUE_RAMP[4], BLUE_RAMP[6], ACCENT, INK_MUTED]
    for (label, tag, desc), c in zip(GRIDS, cols):
        if label in chains:
            ax.hist(chains[label][:, 0], bins=np.arange(15, 45.01, 1.0), histtype="step", density=True, lw=2 if label == "default" else 1.5,
                    color=c, label=f"{desc} (T = {[r for r in rows if r['tag'] == label][0]['sigma_2']:.0f})")
    if args.truth:
        ax.axvline(args.truth[0], color=INK, ls="--", lw=1, label=f"truth {args.truth[0]:g}")
    ax.set(xlabel=r"$|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$]", ylabel="posterior density", yticks=[])
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title(f"{args.run}: pattern-speed posterior per Fourier grid, each at its own T = σ₂ (jitter)", fontsize=10, color=INK)
    savefig(fig, f"{args.run}_pattern_speed_posterior_by_grid")


if __name__ == "__main__":
    main()
