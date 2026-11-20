"""Does a potential-family mismatch flatten the likelihood? Profiles in |Omega_b| of the tempered logL (best evaluated
point per bin, no emulator) for:
  t = 7 mock in its own family      test2_t7_sf_12p5pct            (focused design, 257 points)
  t = 7 mock in Hunter+24           test2_t7_sf_12p5pct_mwH24_osargdesign (the OSARG design, 384 points)
  OSARG data in Hunter+24           osarg_zhang2024_mw_hunter24     (384 points)
all scheme B, sfratio, K = 10, 50 Myr. The last two share design, family and pipeline; only the data differ.

Usage: python plot_family_mismatch_test.py
Figure: figures/test2_t7_family_mismatch_omega_profile_vs_osarg.png
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP

B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
TRUTH_T7 = 25.20
CASES = [("t = 7 mock, own family", "test2_t7_sf_12p5pct", BLUE_RAMP[5], "o"),
         ("t = 7 mock, Hunter+24 family", "test2_t7_sf_12p5pct_mwH24_osargdesign", ACCENT, "s"),
         ("OSARG data, Hunter+24 family", "osarg_zhang2024_mw_hunter24", INK_2, "D")]


def profile(ev, col, edges):
    L = ev[col]
    return np.array([L[(ev.Omega_abs >= a) & (ev.Omega_abs < b)].max() - L.max()
                     if ((ev.Omega_abs >= a) & (ev.Omega_abs < b)).any() else np.nan
                     for a, b in zip(edges[:-1], edges[1:])])


def main():
    set_style()
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    edges = np.arange(15.0, 45.01, 1.25)
    mid = 0.5 * (edges[1:] + edges[:-1])
    for lab, run, c, mk in CASES:
        ev = pd.read_csv(RESULTS_DIR / run / "evals.csv")
        ev = ev[ev.kind == "qmc"]
        col = f"logLs_{B}_sfratio"
        prof = profile(ev, col, edges)
        r = json.load(open(RESULTS_DIR / run / f"emulator_results_{col}.json"))
        p = r["posterior"]
        i = ev[col].idxmax()
        print(f"{lab:30s} N_pts {len(ev)}; best chi2/N {ev[f'chi2N_{B}_sfratio'].min():.3f} "
              f"(median {ev[f'chi2N_{B}_sfratio'].median():.3f}); logLs span {ev[col].max() - ev[col].min():.1f}; "
              f"best point Omega {ev.Omega_abs[i]:.1f} phi {ev.phi_b_deg[i]:.1f} S {ev.S[i]:.2f} f {ev.bar_frac[i]:.2f}")
        print("   posterior: " + ", ".join(f"{k} {v['p50']:.2f} ({v['p16']:.2f}-{v['p84']:.2f})" for k, v in p.items())
              + f"; LOO {r['loo_rms']:.2f}")
        print("   Omega profile per 1.25:", " ".join(f"{x:6.1f}" for x in prof))
        ax.plot(mid, prof, "-", marker=mk, ms=5, color=c, lw=1.6, label=f"{lab}  (best χ²/N "
                f"{ev[f'chi2N_{B}_sfratio'].min():.2f})")
    ax.axvline(TRUTH_T7, color=BLUE_RAMP[3], ls="--", lw=1.2, label=f"t = 7 truth {TRUTH_T7:.2f}")
    for lev in (-0.5, -2.0, -4.5):
        ax.axhline(lev, color=INK_2, lw=0.6, ls=":")
    ax.set(xlim=(15, 45), ylim=(-30, 1), xlabel=r"$|\Omega_b|$  [km s$^{-1}$ kpc$^{-1}$]",
           ylabel=r"best evaluated $\log L_s$ in the bin $-$ overall best")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=9.5)
    ax.set_title("Pattern-speed profile: does fitting with the wrong potential family flatten it?", color=INK,
                 fontsize=11)
    savefig(fig, "test2_t7_family_mismatch_omega_profile_vs_osarg")


if __name__ == "__main__":
    main()
