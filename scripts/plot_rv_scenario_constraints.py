"""Pattern-speed constraint of the t = 7 mock thinned to Gaia-RV sub-samples (run_rv_scenarios_t7.sh), against the
full mock and the real OSARG run.

For each sub-sample: the emulator posterior of |Omega_b| (scheme B, sfratio, K = 10, 50 Myr), and, straight from
the 257 evaluated points with no emulator, the profile of the tempered logL in |Omega_b| (best point per bin) and
its span. Truth: the simulation table's finite-difference pattern speed at t = 7 (25.20 km/s/kpc).
OSARG (real data, Hunter+24 family, its own 384-point box) is drawn for reference only: different potential family,
design box and data.

Usage: python plot_rv_scenario_constraints.py
Figure: figures/test2_t7_pattern_speed_constraint_vs_gaia_rv_subsample.png
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
RUN = "test2_t7_sf_12p5pct"
STEM = "nbody_mock_selection_sigmaDM0.271_N50k"
TRUTH = 25.20
SCEN = [("full mock", RUN, "sfratio"), ("G_RVS < 14", f"{RUN}_grvs14", "sfratio_grvs14"),
        ("G < 16", f"{RUN}_g16", "sfratio_g16"), ("G < 15", f"{RUN}_g15", "sfratio_g15"),
        ("G < 14", f"{RUN}_g14", "sfratio_g14"), ("OSARG coverage, OSARG N", f"{RUN}_osarg_cover", "sfratio_osarg_cover"),
        ("random, OSARG N", f"{RUN}_random_osargN", "sfratio_random_osargN")]
MASK_NAME = {"full mock": "full mock", "G_RVS < 14": "grvs14", "G < 16": "g16", "G < 15": "g15", "G < 14": "g14",
             "OSARG coverage, OSARG N": "osarg_cover", "random, OSARG N": "random_osargN"}


def profile(df, col, lo, hi, nbin=8):
    L = df[col]
    e = np.linspace(lo, hi, nbin + 1)
    return [float(L[(df.Omega_abs >= a) & (df.Omega_abs <= b)].max() - L.max()) for a, b in zip(e[:-1], e[1:])]


def main():
    set_style()
    summ = pd.read_csv(RESULTS_DIR / "rv_selection_masks" / f"{STEM}__summary.csv").set_index("scenario")
    rows = []
    for lab, run, v in SCEN:
        r = json.load(open(RESULTS_DIR / run / f"emulator_results_logLs_{B}_{v}.json"))
        p = r["posterior"]["Omega_abs"]
        ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
        ev = ev[ev.kind == "qmc"]
        col = f"logLs_{B}_{v}"
        prof = profile(ev, col, ev.Omega_abs.min(), ev.Omega_abs.max())
        s = summ.loc[MASK_NAME[lab]]
        rows.append(dict(sample=lab, n=int(s.n), n_R6=int(s.R_lt_6), n_R6_z1to2=int(s.R_lt_6_z1to2),
                         beyond_centre=int(s.beyond_centre), p16=p["p16"], p50=p["p50"], p84=p["p84"],
                         width_pct=100 * (p["p84"] - p["p16"]) / p["p50"], bias_pct=100 * (p["p50"] / TRUTH - 1),
                         loo=r["loo_rms"], span=float(ev[col].max() - ev[col].min()),
                         prof_ends=f"{prof[0]:.1f} / {prof[-1]:.1f}", min_chi2N=float(ev[f"chi2N_{B}_{v}"].min())))
    t = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(t.round(2).to_string(index=False))
    o = json.load(open(RESULTS_DIR / "osarg_zhang2024_mw_hunter24" / f"emulator_results_logLs_{B}_sfratio.json"))
    op = o["posterior"]["Omega_abs"]
    o_width = 100 * (op["p84"] - op["p16"]) / op["p50"]
    print(f"OSARG (real, Hunter+24, 384 points): {op['p50']:.2f} ({op['p16']:.2f}-{op['p84']:.2f}), "
          f"width {o_width:.1f}%, 20,161 stars at R < 6")
    t.to_csv(RESULTS_DIR / "rv_selection_masks" / f"{RUN}__constraints.csv", index=False)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw=dict(wspace=0.3))
    y = np.arange(len(t))[::-1]
    a1.axvline(TRUTH, color=ACCENT, lw=1.4, ls="--", label=f"truth {TRUTH:.2f}")
    a1.errorbar(t.p50, y, xerr=[t.p50 - t.p16, t.p84 - t.p50], fmt="o", color=BLUE_RAMP[5], ms=7, capsize=3, lw=1.6)
    a1.set_yticks(y, [f"{s}  ({n:,} at R<6)" for s, n in zip(t["sample"], t.n_R6)])
    a1.set_xlabel(r"recovered $|\Omega_b|$ [km s$^{-1}$ kpc$^{-1}$] (median, 16–84%)")
    a1.set_ylim(-0.6, len(t) - 1 + 1.1)                 # headroom above the first row for the legend
    a1.legend(loc="upper left", frameon=False)
    a1.set_title("t = 7 mock, 12.5% distance error: thinned to Gaia-RV sub-samples", color=INK, fontsize=11)
    a2.plot(t.n_R6, t.width_pct, "o", color=BLUE_RAMP[5], ms=8)
    # label placement by hand: the points sit close together (offset in points, horizontal alignment)
    place = {"full mock": ((0, 10), "center"), "G_RVS < 14": ((0, 10), "center"), "G < 16": ((0, -16), "center"),
             "G < 15": ((8, -14), "left"), "G < 14": ((0, 10), "center"),
             "OSARG coverage, OSARG N": ((-8, -16), "right"), "random, OSARG N": ((8, 8), "left")}
    for _, r in t.iterrows():
        off, ha = place[r["sample"]]
        a2.annotate(r["sample"], (r.n_R6, r.width_pct), textcoords="offset points", xytext=off, fontsize=8.5,
                    color=INK_2, ha=ha)
    a2.plot(20161, o_width, "D", color=ACCENT, ms=8)
    a2.annotate("OSARG, real data\n(Hunter+24 family, own box)", (20161, o_width), textcoords="offset points",
                xytext=(10, -4), fontsize=8.5, color=ACCENT, va="top")
    nn = np.linspace(t.n_R6.min() * 0.9, t.n_R6.max() * 1.05, 50)
    full = t.iloc[0]
    a2.plot(nn, full.width_pct * np.sqrt(full.n_R6 / nn), color=INK_MUTED, lw=1, ls=":",
            label=r"full mock scaled as $1/\sqrt{N}$")
    a2.set(xlabel="stars at R < 6 kpc", ylabel=r"16–84% width of $|\Omega_b|$ / median  [%]")
    a2.legend(loc="center right", frameon=False)
    a2.set_title("width of the pattern-speed posterior against sample size", color=INK, fontsize=11)
    savefig(fig, "test2_t7_pattern_speed_constraint_vs_gaia_rv_subsample")


if __name__ == "__main__":
    main()
