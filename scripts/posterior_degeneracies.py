"""Posterior degeneracies of every run, measured the same way: are OSARG's different from the mocks', and is that the
data or the potential family? (user, 2026-11-24)

For each run the pipeline's emulator is refitted on its tempered scheme B column (sqrt(2N), as all earlier
posteriors) and sampled as the pipeline. From the chain:
  - the correlation matrix of (|Omega_b|, phi_b, S, bar_frac);
  - fractional 16-84% widths of |Omega_b|, S and bar_frac;
  - the same widths after a phi_b prior N(phi_0, 2 deg) applied by importance weighting (phi_0 = the truth for a mock,
    25 for OSARG and for the steady-state mock), and the effective sample size of those weights.
Cases: the 25 mock galaxies in their own families; t7 at OSARG's coverage and N; t7 in Hunter+24; the steady-state
OSARG-like mock in Hunter+24 (truth 32.5, 25, 1, 1); OSARG in Hunter+24 and in Sormani+22.
Note the boxes differ: the mocks search phi_b 15-39, S 0.6-1.65, bar_frac 0.45-2; the Milky Way runs phi_b 20-30,
S and bar_frac 0.5-1.5. Widths are quoted as fractions of the median to make them comparable.

Usage: python posterior_degeneracies.py [--workers 5]
Output: results/posterior_degeneracies.csv, figures/posterior_degeneracies_mocks_vs_osarg.png
"""
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd
import emcee
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
import test1_gp_emulator as ge

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
NBODY = ["t2", "t2.375", "t2.625", "t3", "t3.375", "t3.625", "t4", "t4.5", "t5", "t5.5", "t6.25", "t7"]
AURIGA = ["au5", "au9", "au10", "au12", "au14", "au17", "au18", "au20", "au22", "au24", "au26", "au27", "au28"]
CASES = ([(g, "N-body", f"test2_{g}_sf_12p5pct", "sfratio", 27.0) for g in NBODY]
         + [(g, "Auriga", f"test2_{g}_sf_12p5pct", "sfratio", 27.0) for g in AURIGA]
         + [("t7, OSARG coverage + N", "t7 variants", "test2_t7_sf_12p5pct_osarg_cover", "sfratio_osarg_cover", 27.0),
            ("t7 in Hunter+24", "t7 variants", "test2_t7_sf_12p5pct_mwH24_osargdesign", "sfratio", 27.0),
            ("steady-state mock, Hunter+24", "Milky Way family", "mw_hunter24_mock_from_osarg_ics", "sfratio", 25.0),
            ("OSARG, Hunter+24", "Milky Way family", "osarg_zhang2024_mw_hunter24", "sfratio", 25.0),
            ("OSARG, Sormani+22", "Milky Way family", "osarg_zhang2024_sormani22", "sfratio", 25.0)])
KEYS = ["Omega_abs", "phi_b_deg", "S", "bar_frac"]


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return np.interp(q, c, x[o])


def one(case):
    lab, group, run, var, phi0 = case
    ge.use_run(run)
    ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
    tr = ev[ev.kind == "qmc"].reset_index(drop=True)
    y = tr[f"logLs_{B}_{var}"].values
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
    P = ge.from_unit(s.get_chain(discard=1000, thin=5, flat=True))
    C = np.corrcoef(P.T)
    w = np.exp(-0.5 * ((P[:, 1] - phi0) / 2.0) ** 2)
    row = dict(case=lab, group=group, run=run, phi0=phi0, ess_prior=float(w.sum() ** 2 / np.sum(w ** 2)) / len(w))
    for i, a in enumerate(KEYS):
        for j, b in enumerate(KEYS):
            if j > i:
                row[f"r_{a}_{b}"] = float(C[i, j])
    # the same correlations with the chain cut to the Milky Way box (phi_b 20-30, S and bar_frac 0.5-1.5): are the mocks'
    # degeneracies a matter of their wider box?
    m = (P[:, 1] >= 20) & (P[:, 1] <= 30) & (P[:, 2] >= 0.5) & (P[:, 2] <= 1.5) & (P[:, 3] >= 0.5) & (P[:, 3] <= 1.5)
    Cb = np.corrcoef(P[m].T) if m.sum() > 200 else np.full((4, 4), np.nan)
    row.update(frac_in_mw_box=float(m.mean()), r_box_phi_f=float(Cb[1, 3]), r_box_phi_S=float(Cb[1, 2]),
               r_box_S_f=float(Cb[2, 3]), r_box_Om_S=float(Cb[0, 2]), r_box_Om_f=float(Cb[0, 3]))
    for k, i in (("Omega", 0), ("phi", 1), ("S", 2), ("f", 3)):
        q = np.percentile(P[:, i], [16, 50, 84])
        qw = [wq(P[:, i], w, x) for x in (0.16, 0.5, 0.84)]
        row[f"{k}_p50"], row[f"{k}_w"] = q[1], q[2] - q[0]
        row[f"{k}_p50_prior"], row[f"{k}_w_prior"] = qw[1], qw[2] - qw[0]
        row[f"{k}_fracw"] = (q[2] - q[0]) / q[1] if k != "phi" else q[2] - q[0]
        row[f"{k}_fracw_prior"] = (qw[2] - qw[0]) / qw[1] if k != "phi" else qw[2] - qw[0]
    print(f"{lab:30s} r(phi,f) {row['r_phi_b_deg_bar_frac']:+.2f} r(phi,S) {row['r_phi_b_deg_S']:+.2f} "
          f"r(Om,S) {row['r_Omega_abs_S']:+.2f} r(Om,f) {row['r_Omega_abs_bar_frac']:+.2f} r(Om,phi) {row['r_Omega_abs_phi_b_deg']:+.2f} "
          f"r(S,f) {row['r_S_bar_frac']:+.2f} | widths Om {row['Omega_fracw']:.2f} S {row['S_fracw']:.2f} f {row['f_fracw']:.2f} "
          f"-> with phi prior Om {row['Omega_fracw_prior']:.2f} S {row['S_fracw_prior']:.2f} f {row['f_fracw_prior']:.2f} "
          f"(ESS {row['ess_prior']:.2f})", flush=True)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--plot-only", action="store_true", help="redraw from results/posterior_degeneracies.csv")
    args = ap.parse_args()
    set_style()
    if args.plot_only:
        t = pd.read_csv(RESULTS_DIR / "posterior_degeneracies.csv")
    else:
        with Pool(args.workers) as pool:
            rows = pool.map(one, CASES)
        t = pd.DataFrame(rows)
        t.to_csv(RESULTS_DIR / "posterior_degeneracies.csv", index=False)
    print("\nmedians by group:")
    cols = ["r_phi_b_deg_bar_frac", "r_phi_b_deg_S", "r_S_bar_frac", "r_Omega_abs_S", "r_Omega_abs_bar_frac",
            "r_Omega_abs_phi_b_deg", "Omega_fracw", "S_fracw", "f_fracw", "S_fracw_prior", "f_fracw_prior", "Omega_fracw_prior"]
    with pd.option_context("display.width", 250, "display.precision", 2, "display.max_columns", 30):
        print(t.groupby("group", sort=False)[cols].median())
        print("within the Milky Way box (phi_b 20-30, S and bar_frac 0.5-1.5):")
        print(t.groupby("group", sort=False)[["frac_in_mw_box", "r_box_phi_f", "r_box_phi_S", "r_box_S_f", "r_box_Om_S",
                                              "r_box_Om_f"]].median())
        t["S_shrink"] = t.S_fracw_prior / t.S_fracw
        t["f_shrink"] = t.f_fracw_prior / t.f_fracw
        print(t.groupby("group", sort=False)[["S_shrink", "f_shrink"]].median())

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    gcol = {"N-body": BLUE_RAMP[5], "Auriga": BLUE_RAMP[3], "t7 variants": INK_MUTED, "Milky Way family": ACCENT}
    for g, s in t.groupby("group", sort=False):
        kw = dict(color=gcol[g], s=55 if g == "Milky Way family" else 30, label=g, edgecolor="white", linewidth=0.6)
        axes[0].scatter(s.r_phi_b_deg_bar_frac, s.r_phi_b_deg_S, **kw)
        axes[1].scatter(s.r_Omega_abs_S, s.r_Omega_abs_bar_frac, **kw)
        axes[2].scatter(s.S_fracw_prior / s.S_fracw, s.f_fracw_prior / s.f_fracw, **kw)
    for _, r in t[t.group.isin(["Milky Way family", "t7 variants"])].iterrows():
        axes[0].annotate(r.case, (r.r_phi_b_deg_bar_frac, r.r_phi_b_deg_S), fontsize=6.5, color=INK_2, xytext=(4, 2),
                         textcoords="offset points")
        if abs(r.S_fracw_prior / r.S_fracw - 1) > 0.04 or abs(r.f_fracw_prior / r.f_fracw - 1) > 0.04:
            axes[2].annotate(r.case, (r.S_fracw_prior / r.S_fracw, r.f_fracw_prior / r.f_fracw), fontsize=6.5,
                             color=INK_2, xytext=(4, 2), textcoords="offset points")
    axes[2].annotate("at (1, 1): OSARG in Hunter+24 and Sormani+22,\nthe steady-state mock, t7 in Hunter+24",
                     (1.0, 1.0), xytext=(0.66, 1.12), fontsize=7.5, color=INK_2,
                     arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8))
    axes[0].set(xlabel="corr(φ_b, bar_frac)", ylabel="corr(φ_b, S)", title="bar angle against bar strength / length")
    axes[1].set(xlabel="corr(|Ω_b|, S)", ylabel="corr(|Ω_b|, bar_frac)", title="pattern speed against bar length / strength")
    axes[2].set(xlabel="width of S with / without the φ_b prior", ylabel="width of bar_frac with / without the φ_b prior",
                title="what a ±2° φ_b prior buys")
    for ax in axes[:2]:
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.axvline(0, color=INK_MUTED, lw=0.6)
    axes[2].axhline(1, color=INK_MUTED, lw=0.6)
    axes[2].axvline(1, color=INK_MUTED, lw=0.6)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Posterior degeneracies, all runs measured the same way (emulator chains, √(2N) tempering)",
                 color=INK, fontsize=11)
    fig.tight_layout()
    savefig(fig, "posterior_degeneracies_mocks_vs_osarg")


if __name__ == "__main__":
    main()
