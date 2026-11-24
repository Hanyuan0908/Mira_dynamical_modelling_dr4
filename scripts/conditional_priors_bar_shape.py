"""Which prior constrains the bar's length and strength: one on the bar angle (as in the mocks) or one on the pattern
speed? (user, 2026-11-24)

In the mocks' own families bar_frac is degenerate with phi_b, so a phi_b prior constrains it; in the Milky Way family
bar_frac and S are degenerate with |Omega_b| instead (posterior_degeneracies.py). Here each run's emulator chain (scheme
B, sqrt(2N) tempering, as before) is importance-weighted with
  flat | phi_b ~ N(phi_0, 2 deg) | |Omega_b| ~ N(Omega_0, 2 km/s/kpc) | both,
and the medians and 16-84% widths of S and bar_frac are compared (truth S = bar_frac = 1 for the mocks).
phi_0 and Omega_0 are the truth for a mock; for OSARG phi_0 = 25 and Omega_0 = 32 or 37.5, two illustrative values.

Usage: python conditional_priors_bar_shape.py
Output: results/conditional_priors_bar_shape.csv and a printed table
"""
from multiprocessing import Pool

import numpy as np
import pandas as pd

from common import RESULTS_DIR
import test1_gp_emulator as ge
import posterior_degeneracies as pdg

CASES = [("t7, own family", "test2_t7_sf_12p5pct", "sfratio", 27.0, 25.20),
         ("t7, OSARG coverage + N", "test2_t7_sf_12p5pct_osarg_cover", "sfratio_osarg_cover", 27.0, 25.20),
         ("au18, own family", "test2_au18_sf_12p5pct", "sfratio", 27.0, 36.27),
         ("steady-state mock, Hunter+24", "mw_hunter24_mock_from_osarg_ics", "sfratio", 25.0, 32.5),
         ("OSARG, Sormani+22 (Omega_0 32)", "osarg_zhang2024_sormani22", "sfratio", 25.0, 32.0),
         ("OSARG, Sormani+22 (Omega_0 37.5)", "osarg_zhang2024_sormani22", "sfratio", 25.0, 37.5)]


def chain(run, var):
    import emcee
    from scipy.optimize import minimize
    ge.use_run(run)
    ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
    tr = ev[ev.kind == "qmc"].reset_index(drop=True)
    y = tr[f"logLs_{pdg.B}_{var}"].values
    y = y - y.max()
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
    return ge.from_unit(s.get_chain(discard=1000, thin=5, flat=True))


def one(case):
    lab, run, var, phi0, om0 = case
    P = chain(run, var)
    wphi = np.exp(-0.5 * ((P[:, 1] - phi0) / 2.0) ** 2)
    wom = np.exp(-0.5 * ((P[:, 0] - om0) / 2.0) ** 2)
    rows = []
    for tag, w in (("flat", np.ones(len(P))), ("phi prior", wphi), ("Omega prior", wom), ("both", wphi * wom)):
        r = dict(case=lab, prior=tag, ess=float(w.sum() ** 2 / np.sum(w ** 2) / len(w)))
        for k, i in (("Omega", 0), ("phi", 1), ("S", 2), ("f", 3)):
            q = [pdg.wq(P[:, i], w, x) for x in (0.16, 0.5, 0.84)]
            r[f"{k}_p50"], r[f"{k}_p16"], r[f"{k}_p84"] = q[1], q[0], q[2]
        rows.append(r)
    return rows


def main():
    with Pool(3) as pool:
        rows = [r for rr in pool.map(one, CASES) for r in rr]
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / "conditional_priors_bar_shape.csv", index=False)
    for lab, s in t.groupby("case", sort=False):
        print(f"== {lab}")
        for _, r in s.iterrows():
            print(f"   {r.prior:12s} S {r.S_p50:.2f} ({r.S_p16:.2f}-{r.S_p84:.2f}) width {r.S_p84 - r.S_p16:.2f} | "
                  f"bar_frac {r.f_p50:.2f} ({r.f_p16:.2f}-{r.f_p84:.2f}) width {r.f_p84 - r.f_p16:.2f} | "
                  f"Omega {r.Omega_p50:.1f} | phi {r.phi_p50:.1f} | ESS {r.ess:.2f}")


if __name__ == "__main__":
    main()
