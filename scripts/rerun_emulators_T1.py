"""Re-run the emulator posterior of a Test 2 search at T = 1 (logL = −Δχ²/2, no tempering), next to the tempered one.

Same emulator and sampler as test1_gp_emulator (quadratic mean + Matérn-5/2 ARD GP + white noise, trained on the
kind == "qmc" rows; emcee 32 walkers, 3000 steps, 1000 burn-in, fixed seeds), on the run's raw scheme B sfratio χ²
(chi2_<B>_sfratio in evals.csv). For each run it saves:
  results/<run>/emulator_results_T1.json   posterior (16/50/84%) for T = 1 and for the tempered scale √(2 N_data),
                                            emulator maximum, kernel and leave-one-out rms (top 32 points)
  results/<run>/posterior_chains_T1.npz     both chains, for drawing potentials (bar strength profiles)
With --phi-prior MEAN SIGMA a Gaussian prior on φ_b is added to the log-posterior and the chain is re-sampled (not
reweighted: at T = 1 the flat-prior chain covers too little of φ_b for importance weights to reach the prior's peak).
Outputs then carry the suffix _phiprior<MEAN>pm<SIGMA>, and the leave-one-out check is skipped (it does not depend
on the prior).
The truth (finite-difference pattern speed for the N-body, own measurement for Auriga) is attached for reference.

Usage: python rerun_emulators_T1.py --runs test2_t7_sf_12p5pct [...]
"""
import argparse
import json

import numpy as np
import pandas as pd
import emcee
from scipy.optimize import minimize

from common import RESULTS_DIR
import test1_gp_emulator as ge

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
LOO_TOP = 32


def fit_and_sample(u, y, phi_prior=None):
    emu = ge.Emulator(u, y)

    def logprior(q):
        if phi_prior is None:
            return 0.0
        phi = ge.from_unit(q[None])[0][1]
        return -0.5 * ((phi - phi_prior[0]) / phi_prior[1]) ** 2

    def logpost(q):
        return -np.inf if np.any(q < 0) | np.any(q > 1) else emu.predict(q)[0] + logprior(q)

    best = None
    for x0 in np.vstack([u[np.argsort(-y)[:8]], np.random.default_rng(0).random((24, 4))]):
        r = minimize(lambda q: -(emu.predict(q)[0] + logprior(q)), x0, bounds=[(0, 1)] * 4, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    s = emcee.EnsembleSampler(32, 4, logpost)
    s.random_state = np.random.RandomState(2).get_state()
    s.run_mcmc(np.clip(best.x + 1e-3 * np.random.default_rng(1).standard_normal((32, 4)), 0, 1), 3000, progress=False)
    loo = None
    if phi_prior is None:
        sel = np.argsort(-y)[:LOO_TOP]
        pred = np.array([ge.Emulator(np.delete(u, i, 0), np.delete(y, i)).predict(u[i:i + 1])[0] for i in sel])
        loo = float(np.sqrt(np.mean((pred - y[sel]) ** 2)))
    return ge.from_unit(s.get_chain(discard=1000, thin=5, flat=True)), ge.from_unit(best.x[None])[0], str(emu.gp.kernel_), loo


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--phi-prior", type=float, nargs=2, metavar=("MEAN", "SIGMA"), default=None)
    args = ap.parse_args()
    suffix = f"_phiprior{args.phi_prior[0]:g}pm{args.phi_prior[1]:g}" if args.phi_prior else ""
    for run in args.runs:
        ge.use_run(run)
        ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
        tr = ev[ev.kind == "qmc"]
        chi2 = tr[f"chi2_{B}_sfratio"].values
        nd = int(tr[f"ndata_{B}_sfratio"].mode().iloc[0])
        u = ge.to_unit(tr)
        out, chains = dict(run=run, column=f"chi2_{B}_sfratio", n_train=len(tr), ndata=nd,
                           truth=ge.RUN["truth"], phi_prior=args.phi_prior), {}
        for tag, scale in (("T1", 1.0), ("tempered", np.sqrt(2 * nd))):
            post, mx, kern, loo = fit_and_sample(u, -(chi2 - chi2.min()) / (2 * scale), args.phi_prior)
            q = np.percentile(post, [16, 50, 84], axis=0)
            out[tag] = dict(scale=float(scale), kernel=kern, loo_rms=loo, emulator_max=dict(zip(ge.KEYS, map(float, mx))),
                            posterior={k: dict(p16=float(a), p50=float(b), p84=float(c))
                                       for k, a, b, c in zip(ge.KEYS, *q)})
            chains[tag] = post.astype(np.float32)
            truth = ge.RUN["truth"]
            print(f"{run} [{tag}] Omega {q[1, 0]:.2f} ({q[0, 0]:.2f}-{q[2, 0]:.2f}), phi {q[1, 1]:.1f}, S {q[1, 2]:.2f}, "
                  f"f {q[1, 3]:.2f}" + (f"; LOO {loo:.2f}" if loo is not None else "")
                  + (f"; truth Omega {truth['Omega_abs']:.2f}" if truth is not None else ""), flush=True)
        json.dump(out, open(RESULTS_DIR / run / f"emulator_results_T1{suffix}.json", "w"), indent=2)
        np.savez_compressed(RESULTS_DIR / run / f"posterior_chains_T1{suffix}.npz", **chains)


if __name__ == "__main__":
    main()
