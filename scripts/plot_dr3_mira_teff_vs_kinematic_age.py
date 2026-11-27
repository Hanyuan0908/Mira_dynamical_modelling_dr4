"""Teff against kinematic age for the DR3 Miras, from both Teff methods, over the full Teff range (user, 2026-11-27).

Data: data/gaia_dr3_lpv_mira.csv (55,074 O-rich Miras).
  - age: the catalogue's `age_kine` column (a function of the pulsation period; it correlates 1.000 with `age_both`).
  - Teff, method 1: `rv_template_teff`, the Teff of the RVS synthetic template used for the radial velocity (only
    stars with a DR3 RV, 3,284; discrete grid, 100 K steps to 4,000 K then 250 K; floor 3,100 K, below which DR3
    removed combined RVs, Katz et al. 2023).
  - Teff, method 2: `teff_gspphot`, GSP-Phot from the BP/RP spectra (14,996 stars).
Panels: one per method, a 2D histogram of Teff against age with the running median and 16-84% band in 1 Gyr age
bins, over the full Teff range up to 8,000 K. A third panel compares the two medians, and a fourth gives the
fraction of each age bin with a template Teff (i.e. with a DR3 RV).
Printed: per age bin, N and the median (16-84%) Teff of each method; the stars above 8,000 K (off the plot).
Usage: python plot_dr3_mira_teff_vs_kinematic_age.py
  -> figures/dr3_mira_teff_vs_kinematic_age_rv_template_and_gspphot.png,
     results/dr3_mira_teff_vs_kinematic_age.csv
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from common import ROOT, RESULTS_DIR, set_style, savefig, INK, INK_2, ACCENT, BLUE_RAMP

METHODS = [("rv_template_teff", "RVS template Teff (stars with a DR3 RV)", BLUE_RAMP[5]),
           ("teff_gspphot", "GSP-Phot Teff (BP/RP spectra)", ACCENT)]
AGE_EDGES = np.arange(1.0, 11.01, 1.0)
TEFF_MAX = 8000.0


def main():
    set_style()
    d = pd.read_csv(ROOT / "data" / "gaia_dr3_lpv_mira.csv", usecols=["age_kine", "rv_template_teff", "teff_gspphot"])
    ac = 0.5 * (AGE_EDGES[:-1] + AGE_EDGES[1:])
    rows = [dict(age_lo=a, age_hi=b, n_all=int(((d.age_kine >= a) & (d.age_kine < b)).sum()))
            for a, b in zip(AGE_EDGES[:-1], AGE_EDGES[1:])]
    fig, ax = plt.subplots(2, 2, figsize=(13, 10), layout="constrained")
    for k, (col, lab, colr) in enumerate(METHODS):
        s = d[["age_kine", col]].dropna()
        a = ax[0, k]
        tb = (np.concatenate([np.arange(3050, 4051, 100), np.arange(4125, TEFF_MAX + 1, 250)]) if col == "rv_template_teff"
              else np.arange(2900, TEFF_MAX + 1, 50))
        h = a.hist2d(s.age_kine, s[col].clip(upper=TEFF_MAX), bins=[np.arange(1, 11.01, 0.25), tb], cmap="Greys",
                     cmin=1, norm="log")
        fig.colorbar(h[3], ax=a, label="stars per cell")
        med, lo, hi = [], [], []
        for r, (x0, x1) in zip(rows, zip(AGE_EDGES[:-1], AGE_EDGES[1:])):
            t = s[col][(s.age_kine >= x0) & (s.age_kine < x1)]
            q = np.percentile(t, [16, 50, 84]) if len(t) >= 10 else [np.nan] * 3
            r[f"n_{col}"], r[f"p16_{col}"], r[f"median_{col}"], r[f"p84_{col}"] = len(t), *q
            lo.append(q[0]); med.append(q[1]); hi.append(q[2])
        a.plot(ac, med, "o-", color=colr, lw=2, ms=5, label="median per 1 Gyr bin")
        a.fill_between(ac, lo, hi, color=colr, alpha=0.18, lw=0, label="16–84%")
        if col == "rv_template_teff":
            a.axhline(3100, color=INK_2, ls=":", lw=1)
            a.text(1.1, 3030, "DR3 floor 3,100 K (cooler RVs removed)", fontsize=8, color=INK_2, va="top")
        rho, p = spearmanr(s.age_kine, s[col])
        nhi = int((s[col] > TEFF_MAX).sum())
        a.set(xlabel="kinematic age, age_kine [Gyr]", ylabel=f"{col} [K]", ylim=(2850, TEFF_MAX + 100), xlim=(1, 11),
              title=f"{lab}\nN = {len(s):,}; Spearman ρ = {rho:+.2f}; {nhi} above {TEFF_MAX:,.0f} K not shown")
        a.legend(frameon=False, fontsize=8, loc="upper right")
        ax[1, 0].plot(ac, med, "o-", color=colr, lw=2, ms=5, label=lab)
        ax[1, 0].fill_between(ac, lo, hi, color=colr, alpha=0.15, lw=0)
        print(f"{col}: N {len(s):,}, Spearman rho {rho:+.3f} (p {p:.1e}); {nhi} stars above {TEFF_MAX:.0f} K: "
              + (", ".join(f"{v:.0f}" for v in np.sort(s[col][s[col] > TEFF_MAX].values)[:12]) if nhi else "none"))
    ax[1, 0].axhline(3100, color=INK_2, ls=":", lw=1)
    ax[1, 0].set(xlabel="kinematic age, age_kine [Gyr]", ylabel="Teff [K], median and 16–84%", xlim=(1, 11),
                 title="The two methods compared (medians per 1 Gyr bin)")
    ax[1, 0].legend(frameon=False, fontsize=8)
    df = pd.DataFrame(rows)
    frac = df["n_rv_template_teff"] / df["n_all"]
    ax[1, 1].bar(ac, 100 * frac, width=0.8, color=BLUE_RAMP[4])
    for x, f, n in zip(ac, frac, df["n_all"]):
        ax[1, 1].text(x, 100 * f + 0.3, f"{100 * f:.1f}%\nof\n{n:,}", ha="center", fontsize=6.5, color=INK)
    ax[1, 1].set(xlabel="kinematic age, age_kine [Gyr]", ylabel="stars with a template Teff (a DR3 RV) [%]", xlim=(1, 11),
                 title="Who has an RVS template Teff: the RV selection by age")
    for a in ax.flat:
        a.tick_params(top=False, right=False)
    fig.suptitle("DR3 Miras: effective temperature against kinematic age, from the RVS template and from GSP-Phot "
                 "(full Teff range to 8,000 K)", fontsize=11, x=0.01, ha="left")
    savefig(fig, "dr3_mira_teff_vs_kinematic_age_rv_template_and_gspphot")
    df.to_csv(RESULTS_DIR / "dr3_mira_teff_vs_kinematic_age.csv", index=False)
    pd.set_option("display.width", 200)
    print(df.round(0).to_string(index=False))


if __name__ == "__main__":
    main()
