"""Is the weak pattern-speed constraint mainly the unconstrained bar strength? (user, 2026-11-25)

If |Omega_b| is degenerate with the bar's strength (and length) in the Milky Way family, constraining the bar should
tighten |Omega_b|. Posteriors at the adopted temperature T1 (sigma_1 of the grid, chi2_jitter_per_grid.py), emulator
and sampler as the pipeline, then importance-weighted with
  none | bar_frac ~ N(1, 0.1) | S ~ N(1, 0.1) | both | both + phi_b ~ N(phi_0, 2 deg)
Runs: OSARG in Sormani+22 and the steady-state mock in Hunter+24 (truth 32.5, 25, 1, 1), on the default grid and on
0.3 kpc rings with one |z| < 2 layer. The prior centre 1 is the published bar (Sormani+22 / Hunter+24 at S = bar_frac
= 1) and, for the mock, the truth.

Usage: python bar_strength_prior_effect.py
Output: results/bar_strength_prior_effect.csv and a printed table
"""
import numpy as np
import pandas as pd

from common import RESULTS_DIR
import test1_gp_emulator as ge
from compare_grid_strategies import posterior
from posterior_degeneracies import wq

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
CASES = [("OSARG, Sormani+22", "osarg_zhang2024_sormani22", "osarg_s22", None),
         ("steady-state mock, Hunter+24", "mw_hunter24_mock_from_osarg_ics", "h24mock", 32.5)]
GRIDS = [("default", "", "0.3 kpc, 3 layers"), ("_r0.3_z2_dz2", "_r0.3_z2_dz2", "0.3 kpc, |z| < 2")]


def main():
    rows = []
    for lab, run, case, om_true in CASES:
        ge.use_run(run)
        ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
        tr = ev[ev.kind == "qmc"].reset_index(drop=True)
        jit = pd.read_csv(RESULTS_DIR / f"chi2_jitter_per_grid_{case}.csv")
        for glabel, tag, gdesc in GRIDS:
            j = jit[jit.grid == glabel]
            T1 = float(np.std(j.chi2.values[1:] - j.chi2.iloc[0], ddof=1))
            chi2 = tr[f"chi2_{B}_sfratio{tag}"].values
            P = posterior(ge.to_unit(tr), -(chi2 - chi2.min()) / (2 * T1))
            g = lambda x, c, s: np.exp(-0.5 * ((x - c) / s) ** 2)
            priors = {"none": np.ones(len(P)), "bar_frac 1 +- 0.1": g(P[:, 3], 1, 0.1), "S 1 +- 0.1": g(P[:, 2], 1, 0.1),
                      "both": g(P[:, 3], 1, 0.1) * g(P[:, 2], 1, 0.1),
                      "both + phi_b 25 +- 2": g(P[:, 3], 1, 0.1) * g(P[:, 2], 1, 0.1) * g(P[:, 1], 25, 2)}
            print(f"== {lab}, {gdesc}, T1 = {T1:.1f}")
            for pr, w in priors.items():
                q = [wq(P[:, 0], w, x) for x in (0.16, 0.5, 0.84)]
                ess = w.sum() ** 2 / np.sum(w ** 2) / len(w)
                r = dict(case=lab, grid=gdesc, T1=T1, prior=pr, ess=ess, Om_p16=q[0], Om_p50=q[1], Om_p84=q[2],
                         Om_width_pct=100 * (q[2] - q[0]) / q[1])
                for k, i in (("S", 2), ("f", 3), ("phi", 1)):
                    qq = [wq(P[:, i], w, x) for x in (0.16, 0.5, 0.84)]
                    r[f"{k}_p50"], r[f"{k}_p16"], r[f"{k}_p84"] = qq[1], qq[0], qq[2]
                if om_true:
                    r["Om_bias_pct"] = 100 * (q[1] / om_true - 1)
                rows.append(r)
                print(f"   {pr:22s} |Omega_b| {q[1]:.2f} ({q[0]:.2f}-{q[2]:.2f}) width {r['Om_width_pct']:.1f}%"
                      + (f" bias {r['Om_bias_pct']:+.1f}%" if om_true else "")
                      + f" | S {r['S_p50']:.2f} ({r['S_p16']:.2f}-{r['S_p84']:.2f}) | bar_frac {r['f_p50']:.2f} "
                      f"({r['f_p16']:.2f}-{r['f_p84']:.2f}) | ESS {ess:.2f}", flush=True)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "bar_strength_prior_effect.csv", index=False)


if __name__ == "__main__":
    main()
