"""Apple-to-apple: the real OSARG run against mocks run the same way. Do the PREDICTIONS move, and do the DATA follow?

For every case the saved library coefficients of all design points are whitened with the data covariance (scheme B,
Hartlap; delta_chi2_bootstrap.whitened, which reproduces evals.csv). For a design point i and the best point b:
  Delta_i^2   = |y_i - y_b|^2          how far apart the two PREDICTIONS are, in units of the data errors
  dchi2_i     = chi2_i - chi2_b        how much worse the DATA find point i
If the data were the best model plus noise, dchi2_i = Delta_i^2 +- 2 Delta_i (Wang et al. 2026), so the ratio
dchi2 / Delta^2 measures what fraction of the model differences the data register. Library noise inflates Delta^2 by
about its floor (~180-230), which pushes the ratio below 1 when Delta^2 is small.
Points are grouped by their fractional pattern-speed offset from the best point, |Omega_i / Omega_b - 1|, in
5-10%, 10-20% and 20-30%, faster and slower separately (other parameters free, as the design gives them).

Cases (all scheme B, sfratio, K = 10, 50 Myr):
  OSARG data in Hunter+24 (384-point OSARG design)
  OSARG-like steady-state mock in Hunter+24 (the same design; truth |Omega_b| 32.5)
  t = 7 mock in Hunter+24 (the same design; wrong family)
  t = 7 mock in its own family, OSARG coverage + N (its own 257-point design)
  t = 7 mock in its own family, all 50k stars

Usage: python compare_data_reaction_osarg_vs_mocks.py
Output: results/data_reaction_osarg_vs_mocks.csv, figures/data_reaction_osarg_vs_mocks.png
"""
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import RESULTS_DIR, set_style, savefig, INK, INK_2, INK_MUTED, ACCENT, BLUE_RAMP
from delta_chi2_bootstrap import whitened

LOCAL = "/Users/hanyuan/mira_dr4_runs"
USB = "/Volumes/USB/hanyuan/projects/mira_dynamics_dr4/Mira_dynamical_modelling_dr4/runs"
STEM = "nbody_mock_selection_sigmaDM0.271_N50k"
B = "fourier_stars_R0to6_vR_sigR_sigphi_m24_vphi_m024"
CASES = [  # label, run, score dir, variant, mock, mask, truth |Omega_b|, colour
    ("OSARG data, Hunter+24", "osarg_zhang2024_mw_hunter24", f"{LOCAL}/osarg_zhang2024_mw_hunter24", "sfratio",
     "osarg_zhang2024_observed_R8.hdf5", None, None, INK),
    ("OSARG-like steady-state mock, Hunter+24", "mw_hunter24_mock_from_osarg_ics",
     f"{LOCAL}/mw_hunter24_mock_from_osarg_ics", "sfratio", "mw_hunter24_mock_from_osarg_ics_sigmaDM0.208.hdf5", None,
     32.5, ACCENT),
    ("t = 7 mock, Hunter+24 (wrong family)", "test2_t7_sf_12p5pct_mwH24_osargdesign",
     f"{USB}/test2_t7_sf_12p5pct_mwH24_osargdesign", "sfratio", f"{STEM}.hdf5", None, 25.20, INK_MUTED),
    ("t = 7 mock, own family, OSARG coverage + N", "test2_t7_sf_12p5pct_osarg_cover", f"{USB}/test2_t7_sf_12p5pct",
     "sfratio_osarg_cover", f"{STEM}.hdf5", f"results/rv_selection_masks/{STEM}__osarg_cover.npy", 25.20, BLUE_RAMP[3]),
    ("t = 7 mock, own family, all stars", "test2_t7_sf_12p5pct", f"{USB}/test2_t7_sf_12p5pct", "sfratio",
     f"{STEM}.hdf5", None, 25.20, BLUE_RAMP[6]),
]
WINDOWS = [(0.05, 0.10), (0.10, 0.20), (0.20, 0.30)]


def emulator_posterior(run, var):
    for name in (f"emulator_results_logLs_{B}_{var}.json", f"emulator_results_logLs_{B}_sfratio.json"):
        p = RESULTS_DIR / run / name
        if p.exists():
            r = json.load(open(p))
            post = r.get("posterior", r)
            o = post["Omega_abs"]
            return o["p16"], o["p50"], o["p84"]
    return (np.nan,) * 3


def main():
    set_style()
    rows, curves = [], {}
    for lab, run, sdir, var, mock, mask, truth, col in CASES:
        design, ev, Y, y_d, chi2, saved, ndim = whitened(run, sdir, var, mock, mask)
        om = design.loc[ev.index, "Omega_abs"].values
        ib = int(np.argmin(chi2))
        d2 = ((Y - Y[ib]) ** 2).sum(1)
        dchi2 = chi2 - chi2[ib]
        frac = om / om[ib] - 1
        temp = 2 * np.sqrt(2 * ndim)
        p16, p50, p84 = emulator_posterior(run, var)
        row = dict(case=lab, n_points=len(Y), N_data=ndim, best_chi2_per_N=chi2[ib] / ndim, best_Omega=om[ib],
                   truth=truth, post_p16=p16, post_p50=p50, post_p84=p84, post_width_pct=100 * (p84 - p16) / p50,
                   chi2_reconstruction_maxdiff=float(np.max(np.abs(chi2 - saved))))
        for lo, hi in WINDOWS:
            for side, sgn in (("faster", 1), ("slower", -1)):
                m = (sgn * frac >= lo) & (sgn * frac < hi)
                key = f"{side}_{int(100 * lo)}to{int(100 * hi)}"
                row[f"n_{key}"] = int(m.sum())
                row[f"Delta2_med_tempered_{key}"] = float(np.median(d2[m]) / temp) if m.any() else np.nan
                row[f"dchi2_med_tempered_{key}"] = float(np.median(dchi2[m]) / temp) if m.any() else np.nan
                row[f"ratio_med_{key}"] = float(np.median(dchi2[m] / d2[m])) if m.any() else np.nan
        rows.append(row)
        curves[lab] = (frac, d2 / temp, dchi2 / temp, col)
        print(f"== {lab}: {len(Y)} points, N {ndim}, chi2 reconstruction max |diff| {row['chi2_reconstruction_maxdiff']:.2f}; "
              f"best chi2/N {row['best_chi2_per_N']:.3f} at Omega {om[ib]:.2f}; posterior {p50:.2f} ({p16:.2f}-{p84:.2f}), "
              f"width {row['post_width_pct']:.1f}%" + (f"; truth {truth}" if truth else ""))
        for lo, hi in WINDOWS:
            for side in ("slower", "faster"):
                key = f"{side}_{int(100 * lo)}to{int(100 * hi)}"
                print(f"   {side:6s} {int(100 * lo):2d}-{int(100 * hi):2d}%: n {row[f'n_{key}']:3d}; model separation Delta^2 "
                      f"{row[f'Delta2_med_tempered_{key}']:6.2f}, data dchi2 {row[f'dchi2_med_tempered_{key}']:6.2f} "
                      f"(tempered medians); dchi2/Delta^2 {row[f'ratio_med_{key}']:.2f}", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS_DIR / "data_reaction_osarg_vs_mocks.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
    bins = np.arange(-0.5, 0.51, 0.05)
    for lab, (frac, d2, dc, col) in curves.items():
        mid = 0.5 * (bins[1:] + bins[:-1])
        for ax, val in ((axes[0], d2), (axes[1], dc)):
            med = [np.median(val[(frac >= a) & (frac < b)]) if np.sum((frac >= a) & (frac < b)) >= 3 else np.nan
                   for a, b in zip(bins[:-1], bins[1:])]
            ax.plot(100 * mid, med, "o-", color=col, lw=1.8 if "OSARG" in lab else 1.2, ms=4, label=lab)
    axes[0].set(xlabel=r"$|\Omega_b| / |\Omega_b|_{\rm best} - 1$ [%]", ylabel="median over design points, tempered",
                title="how far the PREDICTIONS move from the best point (Δ²)")
    axes[1].set(xlabel=r"$|\Omega_b| / |\Omega_b|_{\rm best} - 1$ [%]", title="how much worse the DATA find them (Δχ²)")
    for ax in axes:
        ax.set_ylim(0, 36)
        ax.axhspan(2.0, 2.9, color=INK_MUTED, alpha=0.3, lw=0)
    axes[0].text(-48, 3.2, "library noise floor between near-identical potentials (2.0–2.9)", fontsize=8, color=INK_2)
    axes[0].legend(frameon=False, fontsize=8, loc="upper center")
    fig.suptitle("Apple-to-apple: real OSARG against mocks run through the same pipeline (all other parameters free "
                 "as the design gives them)", color=INK, fontsize=11)
    fig.tight_layout()
    savefig(fig, "data_reaction_osarg_vs_mocks")


if __name__ == "__main__":
    main()
